---
title: Monte Carlo archive reference
description: Orbital and atmospheric archives produced by the desktop Risk Monitor.
---

This page describes the files produced by the desktop application's Risk
Monitor. Orbital Monte Carlo writes `result.json`; Atmospheric Monte Carlo
reads a suitable orbital result and prepares a separate set of entry states.

## Directory layout

For an object with both calculations, the working folder contains:

```text
<object>/
  MC/run-<id>/
    result.json
    report.html
    work/                         # source-solution exports
  Atmosphere/run-<id>/
    entries.json
    source.json
    launch.json
    outcome.json                  # when saved
    environment-input.json        # when a saved environment was selected
    results/
      entries.json
      environment.json
      job.json
      manifest.json
      summary.json
      runs.jsonl
      impacts.csv
      breakups.jsonl
      maps/                       # when map output is produced
```

The `work` directory contains dated source-solution exports, including OEF
files when available.

The contents depend on the run and its completion status. For example, a
failed preparation or map export need not produce every file shown. The MC
archive is written in a `.partial` directory before it is published under its
final run name.

## Orbital `result.json`

The current exporter writes `schema: "orbital-interval-mc-4"`. Its epochs are
MJD TDB. The sampled initial states use heliocentric ecliptic J2000 coordinates,
with positions in AU and velocities in AU/day. These conventions are also
recorded in `time_scale`, `frame` and `state_units`.

### Calculation settings

| Field | Meaning |
| --- | --- |
| `schema`, `time_scale`, `frame`, `state_units` | Archive version and conventions for the orbital states. |
| `method` | `covariance` or `observations`. |
| `target`, `altitude_km` | Target-body enum value and crossing-boundary altitude. Earth is `2`; Moon is `8`. |
| `start`, `end` | Inclusive search interval, in MJD TDB. |
| `requested`, `seed`, `workers`, `integrator` | Sampling and propagation settings. |
| `earth_frame`, `body_fixed_frame`, `body_fixed_model` | Identifiers for the saved Earth and target-body frame models. |
| `forces` | Effective propagation settings: `step_days`, `planet_mask`, `sb441`, `relativity`, `sun_j2`, `earth_j2`, `earth_j3`, `earth_j4`, `earth_rotating_pole`, `target_naif_id` and `interior_mass_policy`. |
| `weight_profile`, `use_weights`, `accepted_indices` | Fit weighting settings and the indices of observations accepted for the source fit. |
| `fit_parameter_covariance` | The source fit's objective-scaled parameter covariance, serialized as rows. It includes correlations between the Cartesian state and fitted non-gravitational parameters, when present. |
| `event_grouping`, `event_max_gap_days` | The convention and maximum consecutive-crossing gap used to group events for presentation. |

Enum values such as `target`, `weight_profile` and force-policy codes are
integers, not descriptive names; their definitions are in the corresponding
OrbitalCore headers. The force fields do not by themselves identify every
external data file used in a run. The frame and model names do not record
the hashes of the PCK files used for the calculation.

The `relativity` field is `full_eih`, `sun_outer_eih` or `off`. The Earth-pole
field selects the pole for the Earth zonal terms, not the coordinate frame for
boundary crossings. Archives created before these settings were added may not
contain the corresponding fields.

### Status and results

| Field | Meaning |
| --- | --- |
| `cancelled` | Whether the overall run was cancelled. |
| `completed` | Number of complete trials. |
| `failed_fit`, `failed_propagation`, `cancelled_clones`, `unstarted` | Counts of trials with each non-complete outcome. |
| `entry_clones`, `entries_before_interval` | Counts of recorded boundary-entry trials and entries already before the requested interval. |
| `probability_denominator`, `uncertainty_interval` | Recorded conventions: complete trials only, with a conditional 95% Wilson interval. |
| `nominal` | Propagation result for the fitted nominal state, separate from the random trials. |
| `clones` | One record for each requested trial. |
| `events` | Groups of crossing trials, with time bounds, counts and probability estimates. |
| `status_codes` | Names corresponding to the numeric clone-status values. |

The nominal trajectory is not a random trial and does not contribute to the
event probability. Event grouping organizes crossings in time; it does not
change their trajectories.

### Nominal and clone records

`nominal` and each item in `clones` have the same structure:

| Field | Meaning |
| --- | --- |
| `index` | Trial index; the nominal record uses `-1`. |
| `status` | `0` unstarted, `1` complete, `2` fit failed, `3` propagation failed, `4` cancelled. |
| `epoch`, `state` | Epoch and six-component heliocentric ecliptic state at sampling, in MJD TDB and AU/AU per day. |
| `non_grav` | Model enum, `enabled`, `fit_mask`, and values `A1`, `A2`, `A3`, `DT`. |
| `entry_before_interval` | Whether entry occurred before the requested interval. |
| `entry` | Boundary crossing record, or `null` when none was saved. |
| `error` | Diagnostic text for the trial, when present. |

A crossing record contains `epoch` in MJD TDB, `relative_position_km` and
`relative_velocity_km_s` relative to the target body's centre, expressed in
the J2000 ecliptic frame. For Earth crossings, these vectors are geocentric.
The record also contains body-fixed `latitude_deg` and `longitude_deg`; for
Earth, latitude is geodetic and refers to the reference ellipsoid.

Each `events` record contains `start` and `end` in MJD TDB, `entries`,
`probability`, `probability_low` and `probability_high`. Probabilities are
conditional on successfully completed trials; the interval reflects finite
sampling rather than every source of model uncertainty.

## Atmospheric handoff

The desktop application reads the orbital `result.json` directly. It accepts
`orbital-interval-mc-4` with MJD TDB and the older `orbital-interval-mc-3`
with MJD TT, converting its epochs as needed. The source must be a complete,
successful ensemble for Earth at the 100 km boundary, with at least one valid
crossing and a compatible saved Earth-orientation model. Physical misses are
allowed; incomplete or failed trials are not. The application groups valid
crossings into events and uses the selected event for the atmospheric run.
Compatibility is checked against the saved frame and model names, not the
hashes of the underlying PCK files. Keep track of the files used if the orbital
and atmospheric calculations need to be reproduced later.

The selected crossing positions and velocities are transformed into ECEF; the
velocity conversion includes the time derivative of the rotating frame. The
application validates the height, inbound direction and saved coordinates,
then writes the resulting states to the atmospheric run's `entries.json`.
The atmosphere is conditional on these selected entries, not on the entire
set of orbital trials.

### `entries.json`

| Field | Meaning |
| --- | --- |
| `schema_version` | Entry-file schema version, currently `1`. |
| `frame` | `ECEF`. |
| `altitude_km` | Entry boundary altitude, `100`. |
| `entries` | Array of entry states supplied to the atmospheric job. |

Each entry contains `id`, `epoch_mjd_utc`, `position_ecef_m`,
`velocity_ecef_mps` and `weight`. For the current desktop orbital handoff,
each selected crossing is assigned weight `1`. A fixed-entry calculation also
writes an entry file, but its state comes from user input rather than an
orbital archive.

### `source.json` and `launch.json`

`source.json` records provenance: the source schema, conditioning, chosen
event and entry counts where applicable, frame and time conventions, and the
input values. For an orbital archive it also records the source path and
SHA-256, Earth-orientation and velocity-transform descriptions, and any
available source photometry path and hash. Fixed entry uses its own provenance
schema. Paths and hashes identify inputs; the original files and external
datasets are not embedded in `source.json`.

`launch.json` records the atmospheric job configuration: input, environment,
physics, integration, Monte Carlo and map settings, and output paths. If a
saved environment file is chosen, a copy is written as
`environment-input.json`. `outcome.json`, when written, records the launch
label, summary and error text.

## Atmospheric `results`

`results/` contains the numerical output and replay inputs for the
atmospheric calculation. `entries.json`, `environment.json` and `job.json`
capture the prepared inputs and configuration; `manifest.json` and
`summary.json` describe the run status and ensemble. `runs.jsonl`,
`impacts.csv` and `breakups.jsonl` hold per-realization and fragment data.
When available, `maps/field.html` displays the report, while the other map
files contain the numerical density, contours and metadata. See
[Atmospheric Monte Carlo](../../guide/atmospheric-mc/) for the report.

The implementation described here is in
`OrbitalCore/src/export/IntervalMcArchive.cpp`,
`TrajectaUI/src/features/atmosphere/AtmosphereArchive.cpp` and
`AtmosphereExecution.cpp`.
