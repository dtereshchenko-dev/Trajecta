# Trajecta CLI (`TrajectaCLI`)

Command-line orbit fitting and impact analysis using the same Core `FitPipeline`
as the UI. This is the renamed FitPipelineCLI implementation, not the removed
legacy OrbitSolver/CLI11 program.

## Build

From the repository root with the configured MSVC developer shell:

```powershell
cmake --preset Release
cmake --build --preset Release --target TrajectaCLI
& ./build/Desktop_Qt_6_9_3_MSVC2022_64bit-Release-Trajecta/TrajectaCLI/TrajectaCLI.exe --help
```

The C++20 executable links to `OrbitalCore`, with the same public compiler flags
and Eigen ABI settings. On Windows, its build refreshes `OrbitalCore.dll` beside
the executable even when relinking is unnecessary. The shared deployment helper
copies kernels, observatory data, weights and debias tables. The CLI reads PCK,
C04 and `ET-UT.dat` from the same user directory as the desktop
application: `%LOCALAPPDATA%/Trajecta/Trajecta/earth_orientation` on Windows or
`${XDG_DATA_HOME:-$HOME/.local/share}/Trajecta/Trajecta/earth_orientation` on
Linux. The installer places missing files there; for a build without the installer,
copy the three PCK files, two C04 files and `ET-UT.dat` from `data/earth_orientation`.

For a CLI-only Linux build without Qt, configure this subdirectory directly:

```bash
cmake -S TrajectaCLI -B build/cli-release -DCMAKE_BUILD_TYPE=Release -DORBITAL_CORE_BUILD_TESTS=OFF
cmake --build build/cli-release --target TrajectaCLI
./build/cli-release/TrajectaCLI --help
```

Provide the platform's Core dependencies and runtime data first. Configure the
toolchain through CMake options/presets, not a hard-coded path in this module.
Override the data source with `-DORBITAL_CORE_DATA_DIR=/path/to/data` if needed.

## Usage

```text
TrajectaCLI <observations.xml|observations.psv|observations.txt> [options]
```

```powershell
TrajectaCLI observations.xml --residuals residuals.txt
TrajectaCLI observations.xml --lov-impact --survey-integrator ra15 --survey-horizon 30
TrajectaCLI observations.xml --mc 10000 --mc-method covariance --mc-start-mjd 61000 --mc-end-mjd 62000 --survey-integrator ra15
TrajectaCLI observations.xml --mc 1000 --mc-method observations --mc-start-mjd 61000 --mc-end-mjd 62000
TrajectaCLI observations.xml --relativity sun-outer --no-earth-j3 --no-earth-j4 --earth-pole fixed-j2000
```

The force model defaults to full EIH, Earth J2/J3/J4 and a time-varying Earth
gravity pole. Use `--relativity full|sun-outer|off`, `--earth-j2` /
`--no-earth-j2` (likewise for J3 and J4), and `--earth-pole
time-varying|fixed-j2000` to change it. Sun-outer retains planetary terms in
the inner EIH sums. The Earth-pole setting applies to the Earth zonal gravity
model, not the position of an observer. Solar J2 remains enabled as usual.
These settings are used in the fit and subsequent propagation, including close
approaches, LOV, MC and `--state-at`. They do not affect a two-body fit when
`--no-nbody` is selected. CLI options do not read the GUI's `Trajecta.ini`.

Use `--help` for weighting, clipping, non-gravitational fitting, close approaches,
LOV and diagnostic options. `--survey-integrator` selects survey propagation;
the nominal fit retains its existing RA15 policy. OEF/RWO export uses Core's
shared `Exporter` and reports the output directory.

Both MC methods use the same interval engine as Risk Monitor. Times are MJD TDB;
the start must be at or after the fit epoch. If omitted, the start is the fit
epoch and the end is the fit epoch plus `--survey-horizon` years. The target is
selected by `--survey-planet` (Earth, Moon or Mars). The old `--obs-mc` option
remains an alias for `--mc`; by default it samples covariance, not observations.
The archive is now the interval-MC format; `--mc-output` writes Earth-fixed
entries when the target is Earth. Nominal encounter filters do not affect MC.
