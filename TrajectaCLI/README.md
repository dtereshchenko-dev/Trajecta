# Trajecta CLI (`TrajectaCLI`)

Command-line orbit fitting and impact analysis using the same Core `FitPipeline`
as the UI. This is the renamed FitPipelineCLI implementation, not the removed
legacy OrbitSolver/CLI11 program. Its existing options and numerical defaults
are unchanged.

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
copies kernels, observatory data, weights, debias tables and Earth-orientation
data, including `earth_orientation/finals2000A.all`.

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
```

Use `--help` for weighting, clipping, non-gravitational fitting, close approaches,
LOV and diagnostic options. `--survey-integrator` selects survey propagation;
the nominal fit retains its existing RA15 policy. OEF/RWO export uses Core's
shared `Exporter` and reports the output directory.
