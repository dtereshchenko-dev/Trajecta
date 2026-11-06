# Developer Guide

## Project structure

```text
Trajecta/
├── OrbitalCore/  Orbital dynamics and orbit determination
├── TrajectaUI/   Graphical interface
├── TrajectaCLI/  Command-line interface
├── StrewnField/  Atmospheric entry and strewn-field modeling
└── data/         Runtime data and observation examples
```

Core fitting clients include `<orbital_core/...>` and link `OrbitalCore`.
Atmospheric clients use `<strewn_core/...>` and link `StrewnCore`; the optional
`StrewnOrbit` adapter links both libraries. Headers in `src/` are internal,
not a supported third-party interface. These are same-toolchain C++ APIs,
not a stable C ABI across arbitrary compilers.

## Quick start (C++ API)

This example loads MPC or ADES observations, determines an initial orbit and
fits it. It uses planetary perturbations and weighted observations; asteroid
perturbations, debiasing and sigma clipping remain disabled. Settings here are
explicit examples, not a guarantee of the UI defaults.

Link your C++20 target to `OrbitalCore`. When adding a target within this
project's CMake build, use `target_link_libraries(your_target PRIVATE OrbitalCore)`
and `orbital_core_deploy_data(your_target)` to deploy the configured runtime
data. The `OrbitalCore` DLL and its dependencies must also be available at
runtime. The example takes the runtime-data directory as its second argument;
the deployment helper puts these files beside the executable. See
[Runtime data files](#runtime-data-files).

```cpp
#include <orbital_core/DE440Source.h>
#include <orbital_core/FitPipeline.h>
#include <orbital_core/ObservationSet.h>
#include <orbital_core/Types.h>
#include <exception>
#include <filesystem>
#include <iostream>

namespace oc = orbital_core;

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: orbit_fit <observations> <runtime-data-dir>\n";
        return 1;
    }

    try {
        getExecutableDir() = std::filesystem::absolute(argv[2]);
        DE440Source de440;
        if (!de440.isAvailable()) {
            std::cerr << "DE440 not available in the runtime-data directory\n";
            return 1;
        }

        // These are process-global, non-owning pointers.
        setDefaultEphemerisSource(&de440);
        struct ResetSources {
            ~ResetSources() { setDefaultEphemerisSource(nullptr); }
        } resetSources;

        auto observations = oc::loadObservations(argv[1]);
        oc::prepareObservations(observations);

        oc::FitConfig config;
        config.useWeights = true;
        config.nbody.enabled = true;
        config.nbody.useSb441 = false;

        oc::FitRunControl control;
        control.progress = [](const std::string& message) {
            std::clog << message << '\n';
        };
        auto result = oc::FitPipeline{config}.run({observations}, control);
        if (!result.fit.converged || !result.fit.orbit) {
            std::cerr << "Fit did not converge\n";
            return 1;
        }

        std::cout << "a = " << result.fit.orbit->getSemiMajorAxis() << " AU\n";
        std::cout << "Optical RMS = " << result.fit.cost << " arcsec\n";
        if (result.fit.uncertainty.hasElementCovariance) {
            std::cout << "sigma(a) = "
                      << result.fit.uncertainty.derivedSigmas.sigma_a << " AU\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Orbit fit failed: " << e.what() << '\n';
        return 1;
    }
}
```

Run the program with paths to an observation file and to the directory holding
runtime data. The application must set `getExecutableDir()` before loading
sources; it is not automatically initialized by the library. `OBSCODE.dat` is needed
for observation loading; radar data also require `RADCODE.dat`. Weighting uses
the configured profile, whose rules must be deployed separately. To enable
debiasing, set `ObservationPreparationOptions::debias` before preparation and
provide the matching bias data. See [Runtime data files](#runtime-data-files)
for fallbacks and optional data.

The default ephemeris pointers are global and non-owning: keep sources alive
for every calculation that uses them, and do not replace or clear them while
calculations are running. To use SB441, keep a `SpkReader` alive alongside DE440,
check `isAvailable()`, register it with `setDefaultAsteroidSource`, enable
`config.nbody.useSb441`, and clear the pointer before destroying the reader.

If an initial orbit and prepared observations already exist, `FitSession`
performs the final fit directly:

```cpp
oc::FitSession session{config};
auto result = session.optimize({
    .initialOrbit = initialOrbit,
    .observations = observations
}, control);
```

This fragment expects an `Orbit_p` and a `std::vector<ObsNode>`; the latter
must remain alive and unchanged for the synchronous call. `FitRequest` also
accepts a borrowed selection mask and optional warm-start non-gravitational
model. Set `FitRunControl::cancel` to a `std::stop_token` for cooperative
cancellation. Progress callbacks are not retained after the call.

`uncertainty.elementCovariance` is a 6×6 matrix for `(a, e, i, Ω, ω, ν)` at
the output-orbit epoch; check `hasElementCovariance` before using it. The
fit-epoch Cartesian and full-parameter covariances are separate fields; use
the full matrix when non-gravitational parameters and their cross terms matter.

## Impact monitoring

Close Approaches is a separate workspace view. Risk Monitor offers covariance
MC, observation MC, MOV and atmospheric MC. The two orbital MC modes propagate
clones across a selected interval and record physical boundary crossings:
covariance MC samples the fitted parameter covariance, while observation MC
perturbs observations and refits each clone. Atmospheric MC uses Earth entry
states from an orbital MC archive.

MOV uses systematic ranging with grid or cobweb sampling and weighted nodes,
not equally probable MC clones. Its current implementation is limited to
ground-based optical observations, gravity-only RA15, and observation arcs and
prediction horizons of at most 30 days. Cobweb requires a suitable converged
nominal, covariance and curvature SNR; the admissible-region grid can work
without a nominal. See `hazard/ManifoldOfVariations.h` for the current contract.

The interval MC API is declared in `hazard/IntervalMonteCarlo.h`. Both
covariance and observation sampling use the same physical boundary scan.
Covariance sampling includes the full fitted state/non-grav covariance.

The experimental LOV implementation traces returns and refines virtual impactors;
it is still under development and is not an MC sample. Its configuration
includes the coordinate system, coordinate scaling, and uniform-sigma or
probability-based sampling. Use the full-parameter overload when fitted
non-gravitational parameters contribute to uncertainty.

For an accepted fit, use the fit-epoch state and its full covariance:

```cpp
#include <orbital_core/hazard/IntervalMonteCarlo.h>

namespace hz = orbital_core::hazard;

hz::IntervalMcConfig config;
config.method = hz::McMethod::Covariance; // or Observations
config.target = hz::Planet::Earth;        // Earth, Moon or Mars
config.startMjdTdb = fit.fit.fitEpochOrbit->epochMjd();
config.endMjdTdb = config.startMjdTdb + 365.25;
config.samples = 10000;
auto mc = hz::runIntervalMonteCarlo(fit, config);
```

Run hazard analysis only with a valid fit-epoch orbit and available covariance.
State and covariance must share their epoch; do not pair the fit-epoch covariance
with an orbit transported to a different display epoch.

The returned clones carry individual status and boundary-entry states, and
events group crossings by time. For saved solutions, the covariance-only
overload accepts an `OrbitSolution` and explicit force settings. The archive
writer records the calculation and supports atmospheric handoff.

## Atmospheric simulation

`StrewnCore` accepts fixed entry states or an ensemble of atmospheric entries.
It handles parameter MC, fragmentation cascades, ablation and dark flight over
prepared wind/terrain grids. Jobs can acquire and cache GFS weather and terrain,
or use explicit offline inputs. `StrewnCLI` and the Qt worker use the same library.

Results include replay inputs, impacts, breakup events and per-realization
parameters. The optional map output contains count/mass density, HPD contours,
GeoTIFF/GeoJSON and a standalone `maps/field.html` report with a coordinate grid,
trajectory animation, altitude profile, wind plots and print views.
Atmospheric outcomes are conditional on the supplied entries, not a new estimate
of the asteroid's orbital impact probability. See the
[StrewnField guide](StrewnField/README.md) for configuration, inputs and output formats.

## Dependencies

- **C++20** compiler; the current Windows setup uses MSVC 2022 x64
- **CMake** 3.21+
- **Ninja** for the supplied presets
- **Qt** 6.9.3 in the Windows presets, including Quick, QuickControls2, Quick3D,
  StateMachine, Charts, Widgets, Network and Test. The root build includes the UI;
  the numerical libraries do not depend on Qt.
- [**vcpkg**](https://vcpkg.io/) for the libraries below:
  - Eigen 3.3+ within the 3.x series; the Linux preset points to Eigen 3.4.1
  - nlohmann_json
  - libcurl >= 7.85 (StrewnField)
  - GoogleTest (tests only)
  - GDAL >= 3.4 with GRIB, GeoTIFF and JPEG2000 support (StrewnField):
    `vcpkg install "gdal[core,openjpeg,png]:x64-windows"`; see
    [StrewnField build instructions](StrewnField/README.md#build-and-run).
- **OpenMP** support from the compiler toolchain (OrbitalCore)

The ERFA subset and HTML viewer libraries are bundled with their license notices.

## Runtime data files

Large ephemeris and catalog-bias files are ignored by Git. The current CMake
deployment copies the required files below, so they must exist under `data/`
when building targets that deploy runtime data. At runtime, their necessity
depends on the selected calculation and observation types.

| File | Size | Source |
|---|---|---|
| `data/linux_p1550p2650.440` | 97.5 MiB | [JPL DE440 binary](https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/de440/linux_p1550p2650.440) |
| `data/debias/efcc18.dat` | 35.4 MiB | [JPL 2018 debias archive](https://ssd.jpl.nasa.gov/ftp/ssd/debias/debias_2018.tgz); extract `bias.dat` as `efcc18.dat` |
| `data/debias/fcct14.dat` | 25.9 MiB | [JPL 2014 debias archive](https://ssd.jpl.nasa.gov/ftp/ssd/debias/debias_2014.tgz); extract `bias.dat` as `fcct14.dat` |

SB441-N16 is optional and is not deployed by the build. The application can
download [JPL SB441-N16](https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp)
(615.8 MiB) in Data sources. It checks the full SHA-256 and installs the file
atomically under the user's local application data; restart to activate it.
A manually installed `sb441-n16.bsp` beside the executable is still supported.
A previously installed file is not removed by a build.

The repository includes `OBSCODE.dat`, `RADCODE.dat`, `weights/vfcc17.rules`,
the test snapshot and the EOP series in `data/earth_orientation/`.
See the [EOP source notes](data/earth_orientation/README.md) for the bundled
series and their interpretation. Observation examples and EN/RU documentation
are under `data/examples/` and `data/docs/`.

`OBSCODE.dat` is required to load observations. `RADCODE.dat` is needed for
radar observations, not optical-only fits; it is the project's converted radar
catalogue, not the original JPL API response. DE440 is needed for planetary
ephemerides and N-body calculations. SB441 is needed only when massive-asteroid
perturbations are enabled. VFCC17 rules and the matching debias table are used
when their respective options are enabled. The current fitter uses EFCC18 for
both weight profiles; FCCT14 is present for compatibility and inspection.
At runtime the new paths take precedence, with `debias_2018/bias.dat` and
`debias_fcct14/xcat.bias_fcct14` accepted as legacy fallbacks. If VFCC17 rules cannot be loaded,
the current implementation falls back to 1 arcsecond; a missing debias table
leaves observations uncorrected.

`finals2000A.all` can be retrieved from [USNO Bulletin A](https://maia.usno.navy.mil/ser7/finals2000A.all).
It provides the UT1 data used by the optical clock. C04 and its consistency
file support the Earth-pole model; historical `ET-UT.dat` provides another
time-data range. Missing files or dates outside their coverage can cause
fallbacks; see the EOP notes for details. These files have different roles and
must not be substituted by simply renaming another IERS product.

The application can download SB441-N16; other source files are not yet downloaded
or updated by the application. Files already loaded by a running process may
require a restart after replacement.

For the debias tables currently used by the project, SHA-256 after
extraction is `3ad243953731218565cbb25937e87803165ac28e4c948f06e52c6f8b9e1af887`
for FCCT14 and `ef6a50830bb83d8b1161e7acb391eed988ea9708e9fbcccc0c834da6109476e3`
for the 2018 table. Both were checked against their linked JPL archives.
The independently checked JPL SB441-N16 SHA-256 is
`919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90`.

The CMake `orbital_core_deploy_data()` function copies runtime files next to
executables that use it. Override `ORBITAL_CORE_DATA_DIR` when storing them
outside this layout. Strewn weather/terrain caches are acquired per job and
are separate from these orbital datasets.

## Versioning

The Trajecta product version is defined in `cmake/TrajectaVersion.cmake`.
Rebuild after changing it: the UI, HTTP User-Agent, Windows executable metadata
and generated `trajecta-version.txt` use that value. The installer deployment
script checks the generated version against `Trajecta.exe` and fills both
installer XML templates. QML module and data-schema versions are independent.

## Build

From an MSVC x64 developer PowerShell, in the repository root:

```powershell
$env:VSLANG = '1049'
chcp.com 65001 | Out-Null
[Console]::InputEncoding = [Text.UTF8Encoding]::new($false)
[Console]::OutputEncoding = [Text.UTF8Encoding]::new($false)
cmake --preset Release
cmake --build --preset Release
```

`CMakePresets.json` also defines `Debug`, `Profile` and `linux-*` presets.
Windows paths currently target Qt 6.9.3, MSVC 2022 and `C:/vcpkg`;
override machine-specific values with CMake arguments or an ignored
`CMakeUserPresets.json`. Windows Release outputs are under
`build/Desktop_Qt_6_9_3_MSVC2022_64bit-Release/`, in each module's directory.

For a Qt-free standalone atmospheric build, see
[StrewnField build instructions](StrewnField/README.md#build-and-run).

### Tests

```powershell
$env:PATH = 'C:\Qt\6.9.3\msvc2022_64\bin;' + $env:PATH
$env:QT_QPA_PLATFORM = 'offscreen'
ctest --preset Release
ctest --preset Release -R '^Strewn\.'
```

CTest runs Core/Strewn GoogleTest cases and Qt test executables. Coverage includes
orbit fitting, covariance and coordinate transformations, LOV/MC behavior,
fragmentation, weather/raster handling, numerical archives and replay.
Strewn HTML/JavaScript tests live under `StrewnField/tests/`; browser-based tests
have their own runtime requirements. A successful build alone does not run tests
or establish cross-platform numerical equivalence.

### Continuous integration

The Windows CI workflow builds the project and runs CTest on pushes, pull
requests, and manual dispatches. It caches vcpkg binary packages and obtains
DE440 and the two debias tables from JPL or a cache, verifying their SHA-256
hashes each time. An initial run must populate the caches; GitHub may evict
them, and cache access can differ between branches and pull requests.
SB441-N16 is not downloaded, so tests requiring it are skipped. The workflow
does not build or publish an installer and does not establish numerical
equivalence across platforms.

