---
title: Data Sources
description: Runtime ephemerides, observatory catalogues, and corrections.
---

Select **Data sources** in the sidebar to see the datasets available to
Trajecta. Each entry shows its status and, where available, its file name,
record count and date coverage. The data needed for a calculation depend on
the observations and settings you choose.

![Data sources dialog showing connected ephemerides and observatory catalogues](../../assets/data-sources/source-status.png)

The screenshots on this page show an earlier layout; the current dialog also
lists the three Earth PCK files.

## Ephemerides and observatories

| Source | Used for |
| --- | --- |
| Planetary ephemeris: DE440 | Planetary positions, observer geometry and N-body calculations. |
| Massive asteroids: SB441-N16 | Perturbations from 16 massive asteroids, when enabled in the force model. |
| Optical observatories | Resolving observatory codes and locations in optical observations. |
| Radar observatories | Station information for radar observations. Optical observations do not require this catalogue. |

DE440 and the observatory catalogues are included with the installer. SB441-N16
is optional and can be downloaded from **Data sources**. Enable its
perturbations in the calculation settings when you want to use it.

## Earth orientation and observation corrections

| Source | Used for |
| --- | --- |
| IERS C04: orientation and UT1 | UT1 conversion within observed coverage and a fallback for Earth orientation. |
| PCK: historical orientation, recent orientation and long-term prediction | Earth gravity pole, stations, topocentric ephemerides and Earth-fixed crossing coordinates. The predictive file is less accurate in the long term. |
| Historical TT-UT1 | Historical time conversion. |
| Astrometric weights: VFCC17 | Observation weights when the VFCC17 scheme is selected. |
| Astrometric debias: FCCT14 and EFCC18 | Catalogue corrections when the corresponding debias model is selected. |

The time-varying Earth gravity pole uses PCK when it covers the epoch. Outside
PCK coverage, the gravity model falls back to C04; beyond its observed records,
the pole is extrapolated rather than measured. The predictive PCK also contains
predictions, not measured future orientations. These choices can affect near-Earth
trajectory calculations. Station and crossing coordinates use PCK independently
of the gravity-pole setting.

UT1 conversions require observed C04 coverage or the historical TT-UT1 table.
The bundled C04 series is used from 19 January 1972 through 15 January 2023;
earlier records in the file are not loaded. PCK provides orientation for
later epochs, but not an independent UT1 time conversion. If neither PCK nor
C04 is available for the gravity pole, the model uses a nominal pole.
Missing VFCC17 rules fall back to a weight of one arcsecond. Check the selected
weighting and debias settings as well as the source status before a fit.

## Check and update sources

The warning icon beside **Data sources** draws attention to a source that
needs checking, a pending restart or an available data update. Open the
dialog to see the current status.

![Data sources warning in the sidebar](../../assets/data-sources/sidebar-warning.png)

- **Connected** means the source is loaded by the application. A file can also
  be **Verified on disk**, **Missing** or **Invalid**. Some entries report that
  a previously loaded file is no longer on disk.
- Use **Check again** to refresh source status and check the observatory
  catalogues and Earth PCK files for updates. PCK updates are also checked
  when the dialog opens. An **Update** button appears beside each available
  update; files are downloaded only when selected.
- PCK updates are checked separately for the historical, recent and predictive
  files. A downloaded file is checked for ITRF93 coverage and prepared for
  installation when Trajecta next starts. The existing file is kept if the
  download or installation fails.
- Use **Download SB441** if you need massive-asteroid perturbations and the
  ephemeris is not available.
- Restart Trajecta when the dialog asks you to reload changed data. Also check
  date coverage: a connected file may not cover the epoch of your calculation.

When an observatory catalogue has an update, its **Update** button appears in
the corresponding row:

![An available optical observatory catalogue update](../../assets/data-sources/observatory-update.png)

Gaia FPR observations are separate from these runtime sources. For a suitable
object, use **Add Gaia FPR and re-solve** on the solution page.

## Where files are stored

On Windows and Linux, the installer places the three Earth PCK files, both C04 files and `ET-UT.dat`
in the user data directory, alongside the planetary ephemeris and observatory catalogues.
PCK updates are stored there too:

- **Windows:** `%LOCALAPPDATA%\Trajecta\Trajecta\`, in the `ephemerides`,
  `observatories` and `earth_orientation` folders.
- **Linux:** `~/.local/share/Trajecta/Trajecta/` by default. If `XDG_DATA_HOME`
  is set, it replaces `~/.local/share`.

Calculation results are stored separately, in dated folders grouped by object.
On Windows, when the application directory is not writable, Trajecta uses
`%LOCALAPPDATA%\Trajecta\work`. If the application directory is writable,
results can instead be saved in its `work` folder.

On Linux, the default is `~/.local/share/Trajecta/work`, or
`$XDG_DATA_HOME/Trajecta/work` when `XDG_DATA_HOME` is set. To open the folder
for the current solution, use the folder icon in the application's top bar.
