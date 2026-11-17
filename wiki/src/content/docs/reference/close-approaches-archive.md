---
title: Close-approach archive reference
description: Files and numerical fields saved by a close-approach search.
---

This page describes the archive written by the desktop **Close Approaches**
page for a completed nominal search. For using the page and interpreting the
results, see [Ephemerides and close approaches](../../ephemerides-and-close-approaches/#close-approaches).

## Directory layout

Each saved search has its own identifier under the current solution folder:

```text
<object>/
  CloseApproaches/
    run-<id>/
      encounters.txt
      result.json
      report.html
```

| File | Contents |
| --- | --- |
| `encounters.txt` | Human-readable encounter catalogue. Times are displayed in TDB, distances in AU and speeds in km/s. |
| `result.json` | Numerical events, sampled distance curves, search settings and a snapshot of the source solution and force settings. |
| `report.html` | Self-contained interactive report with its own embedded data, styles and script. It works offline. |

The text catalogue rounds event times to the nearest minute and formats its
uncertainties for display. Use the numerical values in `result.json` when
that rounding is unsuitable. Display-unit selections do not change the units
stored in the JSON archive. The files are written individually; a failure
while saving can leave an incomplete run folder.

## `result.json`

The current exporter writes `schema: "close-approaches-v2"`. Epochs
in the numerical archive are **MJD TDB**, even if the interval was entered in
UTC. MJD is a numerical date representation, not a time scale on its own.

### Run fields

| Field | Meaning |
| --- | --- |
| `schema` | Archive schema identifier. |
| `object` | Object designation displayed by the application. |
| `timeScale`, `inputTimeScale` | The stored epoch scale, `TDB`, and the scale selected when entering the interval, `UTC` or `TDB`. |
| `startMjdTdb`, `endMjdTdb` | Requested search interval in MJD TDB. |
| `reportingLimitsAu` | Effective reporting distance in AU for each body. This is an output filter, not a force-model cutoff. |
| `seconds` | Elapsed calculation time in seconds before archive writing. |
| `stoppedAtSurface` | Whether the nominal search detected a surface-crossing case. Later encounters are omitted in the affected propagation direction. |
| `uncertaintyMethod` | Description of the local linear 3-sigma uncertainty method. |
| `source` | Source fit-epoch state, covariance and effective settings. |
| `events` | Array of reported nominal encounters. |
| `distanceSamples` | Array of samples used by the distance plot. |

### `source`

The saved state is heliocentric in the ecliptic J2000 frame, at the fit epoch.
The state is not necessarily at the epoch used to display orbital elements.

| Field | Meaning |
| --- | --- |
| `epochMjdTdb` | Fit epoch in MJD TDB. |
| `stateAuAuDay` | `[x, y, z, vx, vy, vz]`; position is in AU and velocity is in AU/day. |
| `fitParameterCovariance` | Objective-scaled full fit-parameter covariance as an array of matrix rows. |
| `covarianceObjectiveScale` | Objective scale associated with the fitted covariance; the saved covariance already includes this scale. |
| `planetMask`, `stepDays`, `sb441`, `relativity` | Effective body selection, initial integration step in days, massive-asteroid option and relativity mode. The initial step is not a fixed output sampling interval. |
| `sunJ2`, `earthJ2`, `earthJ3`, `earthJ4`, `earthRotatingPole` | Solar and terrestrial zonal-gravity settings and the Earth-pole option. |
| `targetNaifId`, `inactiveInteriorMassPolicy` | Target identifier and numeric interior-mass policy used by the force configuration. |
| `ngModel`, `ngEnabled`, `A1`, `A2`, `A3`, `DT`, `ngSolveMask` | Non-gravitational model, enabled state, coefficients in AU/day², time shift in days and fitted-parameter mask. |

The covariance begins with the six Cartesian state components in the order
shown above, followed by fitted non-gravitational parameters when present.
The fitted-parameter mask determines their layout; see
[Dynamics and time](../dynamics/#non-gravitational-acceleration) and
[Uncertainty](../uncertainty/). The matrix can be empty when no full fit
covariance is available; individual non-finite covariance values are `null`.

In this schema, `ngModel` is `0` for none, `1` for Marsden H2O, `2` for
Marsden CO, `3` for inverse square and `4` for Yarkovsky. `ngSolveMask` is a
bit mask: A1 is `1`, A2 is `2`, A3 is `4`, and DT is `8`. For example, `3`
indicates that A1 and A2 were fitted. The force-policy and body-mask values are
also numeric codes; their definitions are in the corresponding OrbitalCore
headers. The `relativity` value is `full_eih`, `sun_outer_eih` or `off`.

### `events`

Each element is one reported nominal encounter. Distances are measured from
the centre of the named body.

| Field | Meaning |
| --- | --- |
| `body` | Planet or Moon involved in the encounter. |
| `mjdTdb` | Epoch of the nominal closest approach, in MJD TDB. |
| `distanceAu` | Nominal closest distance in AU. |
| `minimum3SigmaAu`, `maximum3SigmaAu` | Local linearized 3-sigma distance bounds in AU, or `null` when unavailable. |
| `time3SigmaDays` | Local linearized 3-sigma time uncertainty in days, or `null` when unavailable. |
| `relativeSpeedKmS` | Nominal relative speed at the encounter, in km/s. |
| `vInfinityKmS` | Nominal hyperbolic excess speed in km/s, or `null` when not available. |
| `uncertaintyNote` | Explanation or warning associated with the uncertainty estimate. |

The distance bounds are obtained from the local 3-sigma B-plane ellipse and
nominal encounter energy. They are not simply the nominal distance plus or
minus three standard deviations, and are not guaranteed limits or impact
probabilities. Strongly nonlinear encounters can require additional sampling.
An interior osculating pericentre may be reported in a surface-crossing case;
it is not a physical trajectory through the body's interior.

### `distanceSamples`

Each row contains **19 numbers** in this order:

```text
[mjdTdb,
 Mercury distance, Venus distance, Earth distance, Mars distance,
 Jupiter distance, Saturn distance, Uranus distance, Neptune distance, Moon distance,
 Mercury Hill radius, Venus Hill radius, Earth Hill radius, Mars Hill radius,
 Jupiter Hill radius, Saturn Hill radius, Uranus Hill radius, Neptune Hill radius, Moon Hill radius]
```

All distances and radii are in AU. The Hill radii are approximations: the
Sun is the primary for planets and Earth is the primary for the Moon.
The samples describe the nominal trajectory, not its uncertainty. The
exporter samples the interval and adds points around reported encounters;
the plotted curve can still miss unreported minima. The `events` array is
calculated separately and is not extracted from the plotted curve. For a
surface-crossing case, `distanceSamples` is empty rather than drawing a
trajectory beyond the surface.

## Reproducibility

The archive records a useful snapshot, but not a complete editable project or
all provenance. In particular, it does not embed the original observations,
the full fit configuration, or the versions and hashes of the external
ephemeris and Earth-orientation files. Preserve those inputs and settings
separately when a calculation must be reproduced. The report embeds its own
copy of the JSON data; editing `result.json` does not update `report.html`.

Implementation: `TrajectaUI/src/features/encounters/CloseApproachesViewModel.cpp`
and `OrbitalCore/src/hazard/CloseApproachSearch.cpp`.
