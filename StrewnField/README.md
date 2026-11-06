# StrewnField C++ engine

The atmospheric solver is a shared C++20 library (`StrewnCore`) with a thin
standalone executable (`StrewnCLI`). Both clients use the same numerical path.
The Qt application calls the library in a worker thread; no process protocol
or server is required. Orbital fitting and clone propagation remain in OrbitalCore.

## Implementation status

Implemented: fixed atmospheric states and weighted entry samples, parameter MC,
cascade fragmentation, ablation, dark flight, adaptive integration, supplied 3D
atmosphere and terrain grids, explicit idealized environment, cancellation,
ordered parallel execution, replayable numerical output, and (in the full repository
build) loading a saved work solution and propagating its existing covariance MC
directly into the atmospheric simulation with one command.
Automatic native GFS/DEM acquisition, height-datum conversion, per-asteroid raw
caches and an explicit offline fallback policy are connected to that command.
Optional count/mass JS stopping and configurable representatives per tail bin
are implemented; fixed-N jobs and one representative remain the defaults.

### Fragment mass closure

New jobs use `schema_version: 2`. The ram-pressure convention remains
`rho_v2`; material strengths and their sampling ranges are unchanged.
`half_rho_v2` remains an explicit alternative requiring consistently defined strengths.

The power-law generator clips the last predicted fragment to the remaining
available mass when that remainder is at least `physics.min_mass_kg`. This
fragment participates in grouping, strength scaling and all later flight/breakup
events. It is not left behind in the mass ledger. A remainder below the cutoff,
or a tail whose next power-law mass is below the cutoff, stays in `residualKg`;
the unresolved tail is not combined into an artificial large fragment.

With M=1000 kg, the single-breakup CC example (dust=.6, largest=.15, beta=.55)
has 117 fragments totaling 400 kg, including a final 0.031322109851249 kg
remainder. The iron example (.2, .45, .35) has 6 totaling 800 kg, including
2.407719106529015 kg. These follow the residual-clipping rule in Moscati et al.
section 2.3, not the 116/5 counts labelled in Figure 2. Counts include multiplicity,
not just the grouped trajectories. This verifies the mass-generation algorithm,
not a physical reconstruction of a particular meteorite fall.

Schema-1 jobs retain the old omitted-remainder behavior for archived replay.
The resolved `physics.retain_mass_remainder` boolean is always saved explicitly;
it overrides either schema's default. New API/UI jobs default to true. Existing
editable UI configs retain their stored behavior; to rerun one with corrected
closure, set this field to true. Immutable archives are never migrated.
Entry, environment and result document schemas remain version 1.

Not implemented: time-varying weather, trajectory-history
sampling, HDF5, and restart checkpoints. Unsupported job modes fail explicitly.

### Orbital MC handoff and UI

Earth covariance MC in Qt has a configurable stop altitude, default 100 km.
Qt and the orbital CLI archive the actual nominal, full covariance, parameter
order and scaling provenance, sampling/propagation settings, and completed
entry states under `<fit>/MC/run-*/`. Epochs are explicitly MJD TDB for the
nominal and MJD UTC for atmospheric states. External ephemeris data are not
embedded. A partial/cancelled run, a failed clone, another body/altitude or
missing/duplicate entry cannot be reported as ready merely because a rounded
probability displays 100%.

The atmospheric panel is available for complete Earth MC with every requested
clone crossing the 100 km boundary and a fitted H for size estimation. It uses
those saved ECEF states without repeating orbital MC. This means 100% entry in
that finite sample, not certainty of meteorites reaching the ground. Nothing
starts automatically: the user edits settings and presses Run atmospheric MC.
The orbital CLI only archives its MC; it does not launch the atmospheric job.

Valid edits are atomically persisted to `MC/run-*/atmosphere.json`, directly
runnable with StrewnCLI. Each launch also saves an immutable configuration and
uses a new `<fit>/StrewnField/runs/run-*/` directory. Weather/terrain caches are
shared in `<fit>/StrewnField/`. Network failure is an explicit error; offline
mode uses the cache and does not silently substitute zero wind/flat terrain.

The per-launch copy in `MC/run-*` is pending-launch provenance: Qt removes it
only after successfully saving `StrewnField/runs/run-*/launch.json`. If preparation
fails, the copy is retained. `atmosphere.json` remains the editable UI configuration;
`launch.json` records the original launch, while `job.json` replays the archived
entry states and environment without downloads.

### Saved files

Open `StrewnField/runs/run-*/maps/field.html` for the self-contained report.
Animation uses a fixed selected speed (1x by default, also 5x, 20x and 100x),
with no automatic acceleration after fragmentation.

| Location | Purpose |
| --- | --- |
| `orbit_*_oef.txt`, `weighted_residuals_rwo.txt` | Orbit/covariance exports and observation residuals |
| `MC/run-*/request.json`, `result.json`, `entries.json` | Orbital MC input, complete results and atmospheric handoff states |
| `MC/run-*/manifest.json`, `handoff.json` | MC status, readiness and handoff contract |
| `StrewnField/runs/run-*/entries.json`, `environment.json`, `job.json` | Self-contained atmospheric replay inputs; keep even when they overlap the orbital archive |
| `impacts.csv`, `runs.jsonl`, `breakups.jsonl` in each atmospheric run | Full-precision ground impacts, per-run parameters/mass balances and breakup events |
| `summary.json`, `manifest.json`, `launch.json` in each atmospheric run | Ensemble stopping/status and original launch provenance |
| `maps/density.tif`, `density.json`, `contours.geojson`, `artifacts.json` | GIS raster/contours, normalization/kernel metadata and rendering provenance |
| `StrewnField/weather/`, `terrain/`, `imagery/` | Shared download caches with integrity/source sidecars, not simulation results |

The caches are not needed to view an existing HTML report or replay a completed
run from its archived environment. Removing them prevents offline preparation of
new runs until the sources are downloaded again. Do not remove cache sidecars
individually or treat HTML preview samples as a replacement for numerical archives.

`Realization timeout (s)` persists as `monte_carlo.run_timeout_s` (default 300,
positive, at most 86400). It measures wall-clock time for one atmospheric
realization, not simulated flight time or the whole MC job. The UI disables the
legacy particle-count limit, including when loading an older editable configuration.
Immutable launch/replay archives are not migrated. Standalone CLI jobs still
honor an explicitly positive `physics.max_particles`; zero or omission disables it.

The default material is OC (ordinary chondrite, strong branch). Selecting OC,
CC (carbonaceous chondrite) or Iron applies the supplied material-table ranges:

| Parameter | OC | CC | Iron |
| --- | --- | --- | --- |
| Density, kg/m3 | 2500-3500 | 2100-3100 | 7000-7800 |
| Albedo min/mode/max | .075/.208/.33 | .02/.06/.18 | .07/.18/.54 |
| Strength, MPa | .9-5 | .25-.7 | 5-50 |
| Dust fraction | 0-.9 | .3-.95 | 0-.25 |
| Largest fragment fraction | .05-.4 | .01-.25 | .4-.9 |
| Beta | .4-.8 | .5-.95 | .1-.4 |

All three use wind/drag multipliers .85-1.15 and ablation 7e-8-9e-8 kg/J.
Albedo is triangular; equal min/mode/max fixes it. The selected material albedo
replaces the initial orbital-risk albedo for atmospheric size sampling; fitted
H is unchanged. The OC weak interval .04-.12 MPa is retained but its mixture
probability remains zero, as agreed; the table supplies no mixing weight.
Manual material changes show Custom, and all numeric ranges survive reload.
Selecting a preset leaves entries, H, integrator, grouping, run controls and
maps unchanged. New UI jobs use alpha .2, three representatives per tail bin,
15 workers, at most 1000 runs, and JS .02 with three consecutive checks.
The convergence extent/cell are derived from entry coverage and stored; an
out-of-grid or incomplete ensemble cannot declare convergence. Display grids
can coarsen automatically to fit memory without changing the KDE kernel.

### Rendering ownership

Production rendering belongs to `StrewnCore`, not CLI or Qt clients.
An optional `maps` object enables weighted count/mass KDE, HPD contours,
A self-contained interactive HTML report, GeoTIFF and GeoJSON after simulation.
`StrewnJob` returns `JobResult { simulation, artifacts, directory }`.
Rendering failure leaves the numerical status and archives intact; its separate
status/error is returned and recorded in the run manifest. Qt is not required
at runtime. Windows is verified; Linux has not been built in this stage.

The public methods are **adaptive** (default) and **fixed** only. Adaptive uses
a weighted covariance pilot with effective-sample-size bandwidth and Abramson
local factors. Scott is not a selectable mode.
The JS stopping grid remains independent of presentation smoothing.

## Boundaries

| Layer | Responsibility |
| --- | --- |
| `api/strewn_core` | Typed requests, owning results, units, cancellation |
| `src/Environment.cpp` | Immutable prepared atmosphere and terrain, geodesy |
| `src/TrajectoryIntegrator.cpp` | Single-fragment flight and event localization |
| `src/FragmentationModel.cpp` | Mass groups, strengths, transverse velocities |
| `src/Simulator.cpp` | Cascades, parameter draws, bounded parallel scheduler |
| `src/Convergence.cpp` | Fixed-grid weighted count/mass KDE and JS stopping |
| `src/FieldDensity.cpp` | Pure count/mass density estimation, Gaussian FFT, HPD thresholds |
| `output/` | Archive import, contours, numerical exports, imagery, interactive HTML |
| `app/JobConfig.cpp` | Versioned JSON configuration |
| `app/JobFiles.cpp` | Environment/state import and result serialization |
| `app/StrewnJob.cpp` | Prepare inputs, execute, finalize the output directory |
| `../Common/Http/HttpClient.cpp` | Shared bounded HTTP(S) transfers, TLS verification, timeout/cancellation |
| `app/EnvironmentLoader.cpp` | Download and validate prepared environment documents, per-asteroid cache |
| `app/AutoEnvironment.cpp` | Resolve coverage, direct Copernicus tiles, geoid, explicit fallback |
| `app/GfsSource.cpp` | Archive/forecast selection, GRIB inventory and pressure-level profiles |
| `app/RasterSource.cpp` | GDAL decoding and geographic raster sampling |
| `app/SourceCache.cpp` | Raw source snapshots, URL/range identity, integrity and atomic publication |
| `orbit/OrbitJob.cpp` | Optional `StrewnOrbit` adapter: saved solution, OrbitalCore MC, atmospheric job |
| `cli/main.cpp` | Command line, progress, Ctrl+C |

The numerical files do not depend on JSON, Qt, network access, or file formats.
There are no virtual dispatches or heap allocations in the integration RHS.
The seven-component state and RK stages have fixed stack storage. Large tables
are shared read-only. Independent realizations use persistent workers; at most
twice the worker count of unconsumed realizations is in flight.

## Build and run

With the repository's configured MSVC developer shell:

```powershell
cmake --preset Release
cmake --build --preset Release --target StrewnCLI StrewnCoreTests --parallel 15
ctest --preset Release -R '^Strewn\.'
& ./build/Desktop_Qt_6_9_3_MSVC2022_64bit-Release/StrewnField/StrewnCLI.exe ./StrewnField/examples/cx1-offline.json
```

The module also builds independently, without Qt or OrbitalCore:

```powershell
cmake -S StrewnField -B build/StrewnField-Release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build/StrewnField-Release --parallel 15
ctest --test-dir build/StrewnField-Release --output-on-failure
```

Dependencies: C++20 compiler, threads, nlohmann_json, libcurl >= 7.85,
GDAL >= 3.4 with GRIB, GeoTIFF and JPEG2000 support, and Eigen3;
GoogleTest when
`STREWN_BUILD_TESTS=ON`. DLL and executable must be deployed together. This is
a same-toolchain C++ interface, not a stable C ABI across arbitrary compilers.
MSVC builds use precise floating point, without fast-math.
Install the raster dependency into the vcpkg instance used by the CMake toolchain:

```powershell
C:/vcpkg/vcpkg.exe install "gdal[core,openjpeg,png]:x64-windows"
```

GRIB and GeoTIFF are included in the core build; OpenJPEG and PNG enable the
corresponding GRIB2 compression formats. The default GIS/database feature bundle
is unnecessary. If vcpkg requests rebuilding existing dependencies, inspect its
installation plan before repeating the command with `--recurse`.
The repository's Windows presets use `C:/vcpkg/installed/x64-windows`, not a
dependency directory under `build/`. GDAL is discovered through its CMake config
package. Set `GDAL_DIR` or `CMAKE_PREFIX_PATH` only when intentionally using a
separate dependency prefix. No native plotting package is required.
Windows builds deploy raster DLLs beside StrewnCore, GDAL's resource tables under
`share/gdal`, and PROJ's database under `share/proj` in the same directory.
StrewnCore locates these resources relative to its own module, independently of
the working directory. Explicit `GDAL_DATA`, `PROJ_DATA` and `PROJ_LIB` settings
retain precedence.
Linux uses the installed GDAL/PROJ resources and system CA certificates.

The example is an offline smoke calculation using one CX1 entry state, D=0.72 m,
density=3300 kg/m3, strong bodies only, zero wind and zero ground height. It is
not a reconstruction with the measured CX1 environment. Output directories must
be new or empty. The generated `job.json` reruns into a `replay` subdirectory.
Exit codes: 0 complete, 3 truncated/failed realizations, 4 artifact export failed, 130 cancelled,
1 job/configuration/IO error, 2 command-line usage error.

### Saved work solution, one command

Map `cell_m` is a requested resolution (25 m by default). `auto_cell_size: true`
increases it by powers of two if grid or FFT memory limits would be exceeded;
KDE bandwidths, fragment weights and raw impacts are unchanged. Set it to false
for a strict grid/error policy. `density.json`, image captions and
`artifacts.json.resolved_options` record the actual cell size, and the artifact
report warns when it changed. `artifacts.json.options` preserves the request.
The resulting per-cell probabilities naturally depend on cell area.

The full-repository build also accepts `examples/cx1-work.json`. Its local fixture
is the existing CX1 work directory under `diagnostics/entry_mc_20261010/`; it is
not a newly fitted or externally downloaded orbit. From the repository root:

```powershell
& ./build/Desktop_Qt_6_9_3_MSVC2022_64bit-Release/StrewnField/StrewnCLI.exe ./StrewnField/examples/cx1-work.json
```

Set `input.path` to another existing `work/<date>/<object>` directory, relative
to the config file or absolute. Without `output_directory`, a unique run directory
is allocated inside that asteroid's `StrewnField/runs`. An explicit output path
remains supported for existing jobs. No
preliminary CLI call or manual conversion of entry states is required. The
example uses 200 orbital clones and 200 atmospheric realizations on 15 workers,
with the same explicit idealized environment as the fixed-state smoke test.
It is an integration check, not a meteorite-field reconstruction.

For the same work input with real NOAA GFS and Copernicus DEM:

```powershell
& ./build/Desktop_Qt_6_9_3_MSVC2022_64bit-Release/StrewnField/StrewnCLI.exe ./StrewnField/examples/cx1-weather.json
```

This preset retains 200 entries, 200 strong-body atmospheric realizations and 15
workers, with the wind multiplier drawn from 0.85-1.15. Size, density, material
priors and fragmentation policy remain explicit; acquiring real weather does
not validate the resulting meteorite-field model against an observed event.

Default layout (the object directory name is used as-is):

```text
work/2026-10-04/2023 VD3/
  orbit_cartesian_oef.txt
  StrewnField/
    weather/<source-key>.bin + <source-key>.json
    terrain/<source-key>.bin + <source-key>.json
    runs/run-<unix-milliseconds>-<sequence>/
      request.json
      orbital_mc.json
      entries.json
      environment.json
      job.json
      manifest.json
      summary.json
      runs.jsonl
      impacts.csv
      breakups.jsonl
```

The shared cache is created only when remote sources are requested. Each successful
run archives the actual prepared environment, so replay does not depend on cache
contents or network availability. No existing run is overwritten. Future generated
maps/HTML belong in the same run directory; their generator is not implemented yet.

`input.kind: "orbit_solution"` requires `orbit_cartesian_oef.txt`, including CAR
state, covariance and dynamical MJD epoch in our ECLM J2000 format. Non-grav
solutions additionally require matching `fit_summary.txt` and `covariance.txt`
with full fitted-parameter covariance and cross terms. Old exports without those
sidecars are rejected for non-grav, not silently treated as gravitational.
No fit is repeated and no extra covariance multiplier is applied.

`input.clones`, `input.seed` and `input.altitude_km` configure the existing orbital
MC (defaults: 1000, 42, 100 km). `input.propagation` requires `horizon_days`,
`planet_mask`, and `sb441`: the old saved solution does not retain enough force
metadata to infer these. Match the original fit's force settings. Optional
settings are `integrator` (`ra15` or `ias15`) and `step_days`.
The orbital MC checks the full interval, without a nominal encounter filter.
DE440/SB441 and the ordinary OrbitalCore runtime data are deployed alongside the
executable. A requested but unavailable ephemeris source is an error.

`monte_carlo.runs` independently sets the atmospheric sample count. Atmosphere
runs select from the successfully crossing orbital samples with replacement,
then draw the configured physical parameters. They are not a second orbital
sampler or a one-to-one clone continuation. Size remains explicit in the
atmospheric configuration; OEF `MAG`'s second number is slope G, not albedo.
Non-grav parameters act only during orbital propagation; the atmospheric solver
receives ECEF position/velocity, MJD UTC, IDs and weights.

No entries, invalid covariance or failed orbital clones stop the job explicitly.
Physical misses are reported
separately and atmospheric results are conditional on the selected entries.
Cancellation is cooperative around the nominal encounter scan and between MC
clones; an ongoing nominal scan is not interrupted internally.

## Library calls

```cpp
#include <strewn_core/StrewnJob.h>

auto spec = strewn::loadJob(configPath);
std::stop_source cancellation;
strewn::RunControl control;
control.cancel = cancellation.get_token();
control.progress = onProgress;
auto result = strewn::StrewnJob{}.run(spec, control);
// result.simulation and result.artifacts have independent status.
```

`Simulator::run` accepts entry states and an optional per-run consumer without
filesystem access. `runOne` accepts fixed physical parameters. The environment
must outlive the simulator and every call. Calls block until workers join.
Progress and result callbacks run serially on the calling thread in run-ID
order. A GUI worker must forward progress to the GUI thread using a queued
signal. Result references are borrowed for the callback duration only.
Callback exceptions stop and join all workers before propagating.

Every realization has a seed derived from the master seed and its ID. Entry
selection, physical parameters, and kicks use separate RNG streams. Worker
count and scheduling do not change results on the same build. NumPy's RNG
sequence and bit-identical results across different math libraries are not
promised. Cancellation preserves a completed ordered prefix; this is not a
resumable integration checkpoint.

## Maps and Re-rendering

Add this optional object to a job; omitted/null retains numerical-only behavior:

```json
"maps": {
  "title": "2023 CX1",
  "method": "adaptive",
  "cell_m": 25,
  "fixed_sigma_m": 60,
  "adaptive_alpha": 0.5,
  "adaptive_limits": [0.2, 5],
  "adaptive_classes": 16,
  "basemap": "esri_imagery",
  "offline": false,
  "image_width": 1600
}
```

`fixed_sigma_m` is the Gaussian standard deviation, not a grid-cell size or
diameter. For `fixed`, it sets the isotropic kernel. For a degenerate adaptive
sample covariance, it is an explicitly reported fallback. Adaptive count and
mass maps have separate weighted bandwidths. Every weight includes multiplicity:
`run_weight * multiplicity`, additionally times surviving mass for the mass map.
Each grid sums to one. These are conditional number/mass distributions, not the
probability that at least one meteorite lands in a cell or Earth-impact probability.

Optional `center_deg: [latitude, longitude]` and `half_size_km: [north, east]`
freeze coverage for comparisons; absent center derives coverage from points and
kernel margins. Out-of-grid impacts or a 99% contour touching the grid edge fail
export explicitly. Grid size is bounded to four million cells; increase cell size
or choose a suitable explicit extent for unusually wide fields. Polar/dateline
fields are not supported. Display bounds are common to both maps, with physical
east/north aspect ratio and the actual kernel and wind multiplier labelled.

Results live in `<run>/maps`: `field.html`, `density.tif`, `density.json`,
`contours.geojson`, `artifacts.json`. PNG generation and Blend2D are removed.
GeoTIFF has two Float64 WGS84 bands, count then mass, north-first; public C++ grids
are south-first. Metadata stores axes, normalization, bandwidth and HPD levels.
GeoJSON line coordinates are longitude/latitude, with quantity and HPD attributes.
GDAL contour input is scaled by the lowest HPD threshold to keep its absolute
interpolation tolerance from dominating small per-cell probabilities.

The HTML contains its display data and vendored Leaflet/Plotly scripts, with no runtime
chart CDN or local server. Its map draws the numeric density grid, geographic HPD contours and
selectable impacts rather than a rendered map picture. Embedded imagery tiles
remain available offline even if the imagery cache is removed. With Esri imagery
and `offline=false`, a Leaflet tile layer loads surrounding areas and higher zoom
levels from Esri. Embedded tiles remain underneath as a fallback on network failure.
`offline=true` and `basemap=none` never enable that network layer.

Ensemble statistics include weighted histograms, mass spectra, JS convergence,
mass budgets, medians and P05/P95. Counts include multiplicity; per-realization
statistics retain zero-survivor runs. Breakup height excludes runs with no
breakup. Equal-weight percentiles use linear interpolation (type 7); unequal
weights interpolate between weighted CDF midpoints.

Any archived realization can be selected. Its preview tracks, complete breakup
table, exact parameters and scaled wind profile are included. The library
replays saved atmospheric states in parallel and validates impact IDs, masses,
locations, times, multiplicities, mass ledgers and breakup events before exposing
traces. Failed replays are explicitly unavailable, without invalidating the KDE.
Display traces have at most 96 samples per displayed fragment; exact birth/end
samples are retained separately. The display budget per realization is
clamp(32000 / completed runs, 32, 500) tracks and impact points. Terminal mass
ranks select tracks with complete ancestry; point previews use mass ranks within
each run. These are labelled visualization subsets, not weighted estimators.
All-run metrics and impact histograms, KDE and HPD use the entire archive, never
the display subsets. The browser does not create a million Leaflet markers or
parse every realization's breakup table at startup. Larger per-run payloads are
Deflate-compressed JSON and decoded only when selected; Float32 previews and
grids additionally use byte-shuffle compression. This uses GDAL and the browser's
native DecompressionStream API, with no runtime plotting dependency. Scientific
archives and propagation limits remain unchanged. Old plain-JSON reports remain
supported by the viewer loader.

The interactive map can overlay the selected realization's ground-projected
trajectory and animate its fragments, with play/pause, time scrubbing, speed,
restart and trajectory fitting. Up to 120 terminal branches are sampled by mass
rank, retaining their ancestry; at most 120 representatives are visible at once.
This display subset is not a new Monte Carlo sample and does not change any
ensemble statistics. Birth, breakup, mass-cutoff and ground events determine
when markers appear or disappear; landed markers remain at the archived endpoint.
Latitude, longitude and air-relative speed are exported with the trace columns.
Playback interpolates compact display samples, not a new integration in JavaScript.
Height and breakup markers have event details. The dark-flight marker follows a
displayed branch's configured air-speed threshold, not an observed end of light
or a simulated flare. Old reports without coordinates require report regeneration.
Animation is interactive-only and pauses when printing or leaving the map.
Enabling Trajectory also shows a synchronized side projection (stacked below
the map on narrow screens). Both views share the exact same particle states,
time, selection and playback controls. Its horizontal coordinate is signed
great-circle along-track distance from entry, using the initial ECEF velocity's
horizontal heading. Lateral displacement is projected out; it is not unsigned
ground range or cumulative path length. Older payloads without an entry velocity
use the first nonzero parent-trace displacement; exactly vertical entries use
an explicitly labelled northward projection. The vertical coordinate retains
WGS84 height, including nonzero terrain endpoints. Axis ranges are fixed for the
selected run and scaled independently. Side-view particles and breakup/dark-flight
markers can be inspected; selecting an event seeks both views to that time.

Wind direction uses a polar plot: north at the top, clockwise bearings TO, radial
WGS84 height in kilometres, and colour for horizontal wind speed in m/s. Calm
samples have no direction; the display does not invent a bearing for them.

The print button supports the current view or the ensemble report plus the
selected realization. Browser printing can save PDF; no files are generated per
realization. Recovery observations are not inferred or embedded automatically.
Printing prepares separate fixed-size SVG charts and a composed map in memory;
it does not resize the live map. The printed map uses embedded imagery, including
offline, and covers the full ensemble rather than the current zoomed viewport.

For a saved complete, nontruncated archive, put just the map object in `maps.json`:

```powershell
StrewnCLI --render path/to/run path/to/new-maps maps.json
```

Or call `renderSavedRun(run, output, loadMapOptions(config), stopToken)` directly.
The output must be new/empty; the source archive is never changed. No propagation
of orbital states is repeated. Atmospheric trace replay is required for old
archives and adds report-generation time. Imagery cache defaults to the new output directory for re-rendering;
jobs use `<asteroid>/StrewnField/imagery` or `<run>/inputs/imagery`.

`basemap=none` avoids network access. Esri tiles use existing libcurl TLS/cache
validation and are embedded as original tile bytes, with no native image decoder.
`offline=true` uses cached tiles only. Unavailable imagery produces a labelled
untextured map and warning; TLS/corrupt-data errors fail export, never silently
disable verification. An empty impact ensemble has status `empty` and no map.
An empty ensemble still receives an HTML statistics report. Legacy `font_file`
and `image_width` settings are accepted but unused; fonts are browser-owned.

## Configuration conventions

Schema version is 1. Paths are relative to the job file. Unknown job keys and
ambiguous input/size models are rejected. A parameter range is either a scalar
or `[min, max]`. Continuous ranges are independent uniform draws per realization.
H/albedo sizing uses a triangular `[min, mode, max]` albedo; diameter and density
determine mass together. There is no implicit albedo when mass or diameter is
supplied directly.

`fixed_state` requires `epoch_mjd_utc` (modified Julian date, UTC) plus either ECEF position/velocity or
geodetic latitude, longitude, WGS84 ellipsoidal height, and ENU velocity. Positions
are metres, velocities m/s. ECEF velocity is the coordinate derivative in the
rotating frame, not an inertial velocity with rotated axes. ENU is east, north,
up, ground-relative. These fields are sufficient for an atmospheric start; no
orbital elements or albedo are required if initial mass is provided.

`state_samples` reads `entries.json` with `schema_version: 1`, `frame: "ECEF"`
and an `entries` array of the same state records. Legacy `r_itrs`, `v_itrs`,
`jd_utc` records and the earlier `epoch_jd_utc` field are also supported on
import only; exactly one epoch field is required. Core states and
all new JSON/CSV archives retain MJD without a round trip through JD. Selection probability is
proportional to nonnegative entry weight. The supplied weights are used once,
not multiplied into sampled impacts a second time.

The OrbitalCore interval MC accepts `IntervalMcConfig::altitudeKm`. Zero and
positive values stop at geodetic height above the reference ellipsoid (WGS84 for
Earth). Entries contain actual crossing points and MJD TDB epochs; the adapter
converts the states to rotating-frame position/velocity and MJD UTC. At positive height the
reported probability describes entry through that boundary, not ground impact.
KDE must use those crossing coordinates without projection from ground points.

The orbital CLI exposes `--mc 200 --mc-altitude-km 100 --mc-output entries.json`.
The Earth states in this file can be used directly as `state_samples` by
StrewnCLI. The file includes run counts and crossing IDs; selection remains
conditional on successful Earth entries in the selected interval.
The frame uses OrbitalCore's IAU1976/1980 UT1/GAST convention without polar motion,
not a new high-precision ITRF realization. The one-command `orbit_solution` adapter
uses the same MC engine.

Defaults are a strong-only preset, not a universal population prior:
`weak_probability=0`, strong strength 0.9-5 MPa, density 2500-3500 kg/m3,
ablation 7e-8 to 9e-8 kg/J, alpha=0.2, independent wind and drag multipliers
0.85-1.15. The initial size default is 1 kg. Weak bodies can be enabled explicitly;
their default strength range is 0.04-0.12 MPa. Small-fragment retention is a
separate parameter: minimum tracked mass is 0.001 kg.

The default dynamic-pressure convention is rho*v_relative^2. Half-rho-v^2 is
available explicitly. First/later breakups retain 30/10 individual fragments
and 4/2 logarithmic tail bins per decade. There is no particle-count cap by
default (`max_particles: 0`). A positive legacy cap counts propagated
representatives, not physical multiplicity; hitting it produces `particle_limit`.
The entire cascade instead shares a cooperative steady-clock deadline controlled
by `monte_carlo.run_timeout_s`. Checks run between particles, inside integration
(including rejected steps), and during mass enumeration/splitting. Enumeration
no longer stops at ten million children. Queuing before a worker starts, network
preparation, archive writing and map rendering do not consume this deadline.
Clock reads are amortized in numerical loops; this is not a hard real-time or
memory limit. Blocking client observers and allocator/OS calls cannot be forcibly
interrupted. Arithmetic validation and numerical integration safeguards remain.

A wall-clock expiration produces `wall_timeout`, retains unresolved mass, counts
as failed (and separately in `summary.wall_timeouts`), and blocks JS convergence
and normalized map publication. No failed realization is silently discarded or
replaced. Unexpired results retain the same random draws and numerical steps;
timeout outcomes themselves can depend on CPU contention and worker count.

`brick_area_only` and `cone_area_only` are experimental fixed-area models sharing
the spherical Cd(M) law. They are not validated attitude-dependent brick/cone
aerodynamics. The luminous entry parent remains spherical; the selected shape
is used for descendants and for a direct dark-flight start.

## Environment

### Automatic native sources

```json
{
  "mode": "automatic",
  "margin_km": 25,
  "terrain_step_deg": 0.02,
  "timeout_seconds": 120,
  "max_download_bytes": 268435456,
  "offline": false,
  "on_unavailable": "error"
}
```

No URLs, weather files, event date or coordinates need to be repeated in this
mode. Coverage is derived from the supplied entries and their straight-line
WGS84 ground intersections, padded by `margin_km`. The projection is only a
download-envelope estimate, never a replacement for numerical impact propagation.
Optional `bounds_deg: [south, north, west, east]` overrides it; leaving downloaded
coverage during flight is an error. Polar/dateline corridors are currently rejected.
Limits: 512 native weather columns, 200000 terrain samples and 100 DEM tiles.
The default 0.02-degree terrain sampling is approximately 2.2 km north-south,
not the native 90 m DEM resolution. Set a smaller step explicitly for rough terrain.

Weather comes **directly from NOAA**, not Open-Meteo:

- [GFS public archive](https://registry.opendata.aws/noaa-gfs-bdp-pds/): 0.25-degree
  GRIB2, selecting HGT/TMP/UGRD/VGRD pressure levels from 1000 to 1 hPa. HTTP byte
  ranges fetch only relevant messages; GDAL decodes their native packing.
- Older events use the 0.5-degree GFS analysis archive at the migrated
  [NCEI object store](https://www.ncei.noaa.gov/oa/prod-model/index.html#global-forecast-system/access/historical/analysis/).
  Historical availability depends on the archive. A supplied `gfs_file` can use
  an already acquired complete GRIB, with its metadata/time validated.
- Choose the nearest three-hour valid time. For recent/future events, probe up
  to eight issued six-hour cycles and reject leads beyond 384 h. Optional
  `cycle_mjd_utc` pins a particular cycle; it is exclusive with `gfs_file`.
  Missing archives/forecasts are not replaced with unrelated days.

The weather snapshot is static, not time-interpolated. Entry times spanning more
than one hour are rejected. The snapshot is within 90 minutes of the midpoint
entry time; its exact offset, issue time and lead are recorded as MJD UTC.
GDAL's possible Celsius conversion is explicitly handled; winds remain east/north
in m/s. Vertical wind is zero, not derived from latitude or the north component.
Density uses the dry ideal-gas relation. Geopotential heights convert to geometric
heights with R=6356766 m and an EGM2008 offset to the WGS84 ellipsoid (an approximate
vertical datum realization). Below the lowest pressure level values are held
constant; above the highest level, US76 is used with a 10 km transition in density,
temperature ratio and wind. The prepared grid extends from -1 to 140 km at 50 m.

Terrain comes **directly from the public Copernicus GLO-90 GeoTIFF tiles**, not
the Open-Meteo elevation API. The immutable 2021 release's tile inventory identifies
ocean-only tiles, for which orthometric height is zero as specified by the
[dataset publisher](https://copernicus-dem-90m.s3.amazonaws.com/readme.html).
A failed land-tile download never becomes zero elevation. Bilinear native-tile
sampling holds the last pixel constant over a missing shared edge, at most one
native pixel. DEM heights are a surface model including vegetation/buildings,
not a bare-earth DTM. [EGM2008](https://github.com/OSGeo/PROJ-data/blob/master/us_nga/us_nga_README.txt)
undulation is downloaded from PROJ/NGA and added: `h_WGS84 = H_EGM2008 + N`.

Raw GRIB messages, the geoid, inventory and DEM tiles are retained under the
asteroid's `weather/terrain` directories. Sidecars record URL, byte range, final
URL, download time, byte count and FNV-1a corruption checks (not cryptographic
authentication). An offline run uses these same records; it never accesses the
network. Archive keys include immutable cycle/lead URLs. Forecasts resolve the
latest available cycle anew unless explicitly pinned; for exact replay use the
generated `job.json` with its frozen environment. Busy caches fail explicitly.

`on_unavailable: "idealized"` explicitly permits the complete US76/zero-wind/flat
surface fallback when the network/source or offline cache is unavailable.
The snapshot records `idealized_fallback` and the reason, and the manifest marks
the run as idealized. Invalid data, wrong times, coverage errors, TLS verification
errors and corrupt cache entries are not downgraded to this fallback.

### Prepared environments

`idealized` explicitly selects the reference US76-based table, zero wind and
constant WGS84 ellipsoidal ground height. `file` requires a prepared schema-1
environment. A missing or invalid file fails; it never silently becomes zero
wind. The `remote` mode downloads two **prepared schema-1 JSON documents**, not
raw GFS GRIB or DEM raster products:

```json
{
  "mode": "remote",
  "weather_url": "https://your-data-host/weather-at-event-time.json",
  "terrain_url": "https://your-data-host/terrain-wgs84.json",
  "timeout_seconds": 60,
  "max_download_bytes": 67108864
}
```

Both documents require `schema_version: 1` and `height_datum: "WGS84"`.
The weather document contains `atmosphere`; the terrain document contains
`terrain`, using the structures described below. These URLs are illustrative,
not operational provider endpoints. Use `automatic` for the native providers above.

The same linked libcurl backend is used on Windows and Linux, with no QtNetwork,
WinHTTP branch, or external curl process. HTTP(S) only, verified TLS, no HTTPS to
HTTP redirect downgrade, bounded response size, timeout and cancellation. The
transport follows the [libcurl threading guidance](https://curl.se/libcurl/c/threadsafe.html)
and [protocol restrictions](https://curl.se/libcurl/c/CURLOPT_PROTOCOLS_STR.html).
On Windows, libcurl uses the native CA store; Linux uses the installed TLS
backend's CA configuration (deploy the distribution's CA certificates).
Certificate and hostname checks are never disabled.

Validated documents are cached by exact URL in `StrewnField/weather` and
`StrewnField/terrain`. URL identity is checked inside the cache record; its filename
hash is a lookup key, not a cryptographic integrity claim. Records retain the
provided source metadata, requested/final URL and download time. Files publish
through a temporary file and rename, guarded by a per-key directory lock. A
busy/stale lock produces an explicit error; no automatic stale-lock removal.
Use immutable URLs identifying source/cycle/valid time: cache freshness for live
forecast endpoints is not inferred. An existing validated record is reused without
network; malformed data and HTTP errors never silently become zero wind. Generic
fixed-state jobs keep their remote inputs under `<output_directory>/inputs`.

The environment snapshot specifies `height_datum: "WGS84"`. Atmosphere axes are
regular `latitudes_deg`, `longitudes_deg`, `first_height_m`, `height_step_m`,
`height_count`. Flattened samples run latitude, longitude, height (height fastest),
each `[density_kg_m3, temperature_k, east_mps, north_mps, up_mps]`. Terrain has its
own regular horizontal axes and `heights_wgs84_m`, longitude fastest. Horizontal
singleton axes describe a uniform column/surface. Interpolation is trilinear for
air and bilinear for terrain. Below the atmosphere table, the first layer is
held constant. Leaving a finite horizontal grid or the atmosphere top fails.
Periodic longitude is not supported. Supplied grids must already have WGS84
heights; only automatic native-source preparation performs geoid conversion.

The wind is an actual ENU flow vector, not a meteorological FROM bearing. Up is
independent of north; no `W = V*sin(latitude)` term exists. A realization's wind
multiplier scales all supplied components coherently. The environment is static
over a realization; source/cycle/valid-time provenance can be retained in the
environment JSON, which is copied into the result snapshot.

## Results and interpretation

| File | Contents |
| --- | --- |
| `job.json` | Resolved defaults and a replay job using local snapshots |
| `entries.json` | Exact atmospheric states and input weights |
| `environment.json` | Exact prepared environment and supplied provenance |
| `manifest.json` | Engine/RNG identity, conditioning, limitations, final status |
| `summary.json` | Counts, elapsed time, stop reason, JS history and blocking diagnostics |
| `runs.jsonl` | Sampled parameters, per-run status/error, mass budget, step counts |
| `impacts.csv` | Full ECEF state, geodetic position, mass, multiplicity, time, IDs |
| `breakups.jsonl` | Parent states, strengths, generation and group counts |
| `orbital_mc.json` | For work input: imported state/full COV, orbital settings, counts, events and crossing clone IDs |
| `request.json` | For work input: original request for provenance, with its original path base |

Replay atmospheric results using `job.json`, not the copied `request.json`:
the former references local frozen entries/environment and does not need the
original work directory. `request.json` preserves relative paths as submitted.

CSV uses 17 significant digits and locale-independent decimal points. JSON
round-trips binary64. Completed records are flushed after each realization;
the manifest distinguishes complete, partial, cancelled, and failed jobs. The
archive is for replay and analysis, not crash recovery mid-fragment.

Mass balance separates landed, ablated, dust, mass-spectrum residual, below-cutoff,
and unresolved mass. Initial mass must equal their sum. Multiplicity counts
represented fragments, not probability. `run_weight=1` is an unnormalized MC
weight. All results are conditioned on the supplied atmospheric entries; this
module does not calculate the probability of an asteroid hitting Earth.

For a normalized number-density map use multiplicity; for a mass-density map use
mass times multiplicity, with the chosen run weighting and normalization stated.
Neither is the probability of at least one fragment landing in a cell. Incomplete
realizations require explicit treatment; silently discarding them can bias maps.
If maps were requested for a partial ensemble, `manifest.json.artifacts` records
`status: "skipped"` and an explicit reason, also returned to the UI/CLI. Particle
limits require a new run with `physics.max_particles: 0` (or a larger explicit
cap), the same seed and input snapshots. Wall-clock timeouts name
`monte_carlo.run_timeout_s` in the report. These computational limits are separate
from the maximum number of Monte Carlo realizations.
Omitting `monte_carlo.convergence` (or setting it to null) retains fixed-N execution.

### JS stopping and directional representatives

For CX1, the tested configuration is:

```json
{
  "physics": {"representatives_per_bin": 3},
  "monte_carlo": {
    "runs": 1000,
    "workers": 15,
    "convergence": {
      "threshold": 0.02,
      "batch_size": 50,
      "consecutive_checks": 3,
      "center_deg": [49.83416935789474, 0.7873627789473683],
      "half_size_km": [14, 20],
      "cell_m": 25,
      "kernel_sigma_m": 60
    }
  }
}
```

Merge these fields into a complete job; this fragment does not define inputs or
environment. `center_deg` is required (latitude, longitude), `half_size_km` is
north/east. Use an appropriate fixed region for other events. The grid origin
does not modify flight states.

The monitor compares consecutive cumulative maps at 50, 100, 150, ... runs.
Both base-2 Jensen-Shannon divergences must be strictly below the threshold
for three consecutive comparisons. This is divergence, not its square root.
The first nonempty batch establishes the reference; empty maps cannot converge.
Count weights are multiplicities, mass weights are mass times multiplicity.
The fixed 25 m grid and isotropic 60 m Gaussian are independent of the adaptive
bandwidths used for final plots. The local coordinate radius is 6367 km.
Gaussian convolution is separable and truncated at eight sigma; no FFT library
or extra process is needed.

Only the ordered run-ID prefix is consumed and archived, independent of worker
count. On convergence, outstanding speculative work is cooperatively stopped
and excluded. `runs` is the upper bound, not an instruction to complete 1000
despite convergence. `stop_reason=js_converged` is distinct from cancellation
and `run_limit`. JS values and stable-check counts are saved in `summary.json`;
a null JS denotes a reference or unevaluable batch. Any incomplete realization
or out-of-grid impact blocks convergence, with explicit summary diagnostics.
The computation then continues to the run cap unless cancelled. Unlike the old
research code, incomplete or spatially clipped ensembles cannot pass the stop gate.

`representatives_per_bin` is 1 by default, range 1..1024. A tail interval with
N fragments gets min(K,N) independently drawn kick directions. All retain its
mean mass; multiplicity is divided equally, with the rounding remainder in the
last representative. Individual large children, bin boundaries, retained mass
and dust/residual policy are unchanged. Fractional multiplicities are statistical
weights. Momentum correction is still applied to the complete child set.
More representatives change kick consumption and subsequent trajectories, but
not the independent entry-selection or physical-parameter random streams.

JS stability measures sampling stability on this grid, not grouping accuracy,
physical validity, or agreement between different grouping resolutions. Adaptive
KDE bandwidth also changes with effective sample size when representatives change;
compare raw spatial moments and maps with a common fixed kernel as well.
