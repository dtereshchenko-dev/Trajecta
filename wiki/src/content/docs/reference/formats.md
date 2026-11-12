---
title: Files and outputs
description: Observation input, solution files and calculation archives.
---

## Inputs

Trajecta reads observation files in MPC 80-column and ADES XML/PSV formats.
Open the file from the application or drag it onto the workspace. The
measurements, station codes, catalogue information and supplied uncertainties
available to the fit depend on the file. See [Observations](../../guide/observations/)
for the supported records and how to review them.

A calculation also uses external data such as planetary ephemerides and,
when enabled, catalogue corrections or asteroid perturbations. These datasets
are separate from the observation file; their status is shown in
[Data sources](../../data-sources/).

## Where results are saved

After a converged orbit fit, Trajecta creates a folder for the object under a
dated `work` directory. The folder button in the top bar opens the current
solution folder. Calculation pages also provide a folder button for the
individual run.

On Windows, the exporter first tries `work` next to the executable. If it cannot
write there, it uses `%LOCALAPPDATA%\Trajecta\work`. An installed application
in a protected directory commonly uses the latter; a development build may use
its own `work` directory. On Linux, results normally go to
`$XDG_DATA_HOME/Trajecta/work`, or `~/.local/share/Trajecta/work` when
`XDG_DATA_HOME` is unset. The folder button shows the actual location.

The date in the folder name is the local date when the result was saved; it is
not the orbit's epoch. Re-solving the same object on the same date overwrites
the automatic solution files in that object folder. After a successful export,
obsolete automatic solution files are removed; separate calculation archives
and other files are kept. Save a copy before re-solving if you need the earlier
solution files. Individual calculation archives use separate `run-…` folders.
Those archives continue to describe their original runs, not the latest fit.
An explicitly exported `covariance.txt` is also kept and may refer to an
earlier solution.

## Orbit solution

The current solution and residual exports use **OrbFit-like OEF and RWO**
formats. They were introduced to make cross-validation with OrbFit easier;
they should not be assumed to support every OrbFit format variant. These
exports are planned to move to JSON in a future release.

The standard fit export includes these files when the corresponding solution
data are available:

| File | Contents |
| --- | --- |
| `orbit_cartesian_oef.txt` | Cartesian state in OEF format, with available covariance and fitted parameters. |
| `orbit_equinoctial_oef.txt` | Equinoctial orbital elements and available uncertainty information. |
| `orbit_keplerian_oef.txt` | Keplerian orbital elements and available derived quantities. |
| `orbit_cometary_oef.txt` | Cometary elements for a solution using the cometary representation. |
| `weighted_residuals_rwo.txt` | Observation records and weighted residual information in RWO format. |

The element files produced depend on the solution representation. A
`fit_summary.txt` is also written when a nominal Earth entry or its calculation
error is available. The core exporter can additionally write a
`covariance.txt` with formal and scaled matrices, but that file is not part of
the standard automatic fit export.

OEF files include an epoch and coordinate representation; use those fields
when comparing solutions. The covariance matrix view in the application is at
the fit epoch, which can differ from the displayed element epoch. See
[Uncertainty](../uncertainty/) for the distinction.

## Ephemerides and close approaches

Each completed ephemeris calculation is saved in
`Ephemerides/run-…` inside the solution folder:

- `ephemeris.txt` is the readable table with calculation settings and units.
- `result.json` contains the numerical result and its metadata.

The text table follows the selected display columns and also includes the
available uncertainty columns. The folder button on the ephemeris page opens
the saved run.

Completed close-approach searches are saved under `CloseApproaches/run-…`.
`encounters.txt` is the readable encounter catalogue; `result.json` records
the numerical result and configuration. The results table can use different
display units, so check the file headers when using the exported values.

## Risk calculations

Each orbital MC run is stored in `MC/run-…`, and each MOV run in `MOV/run-…`.
The run folder contains a `result.json` numerical archive and a `report.html`
for viewing the results. The orbital MC archive also includes a copy of the
source solution information used for that run. Open the report or folder from
the corresponding results view.

An atmospheric run is stored under `Atmosphere/run-…`. Its `entries.json`,
`source.json` and `launch.json` record the input states, their provenance and
the launch configuration. The `results` folder contains the simulation output,
including a summary and the individual-run, impact and breakup data. When maps
are produced, their files and HTML view are under `results/maps`.

For the fields in the orbital and atmospheric handoff archives, see the
[Monte Carlo archive reference](../mc-archive/).

## Reusing results

The core has an OEF reader for supported orbital elements, covariance and
non-gravitational data. Atmospheric MC can also use a suitable saved Earth
entry archive from orbital MC. Keep the archive and its accompanying files
together when moving them.

These readers do not restore a complete editable project with its observations,
settings and previous runs. Full project restoration is a separate feature.
