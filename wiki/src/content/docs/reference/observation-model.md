---
title: Observation model
description: Observer geometry, optical and radar predictions, weights and residuals.
---

The observation model connects an orbit to the measurements used to fit it.
For each observation, Trajecta prepares the reception time and observer geometry,
predicts the measurement from the current trajectory, and forms an observed
minus computed residual. The optimizer uses these residuals and their
uncertainties to adjust the orbit and any fitted non-gravitational parameters.

## Input and preparation

Optical observations provide a direction, normally right ascension and
declination; radar observations can provide round-trip delay, Doppler shift or
both. The data can also include station and catalogue codes, observing
technique, time precision, magnitude and measurement uncertainties.

The input reader converts observations into a common representation. The fit
preparation resolves station data, converts the reception epoch to TDB, obtains
Earth's state from the planetary ephemeris and prepares weights for each
measurement. Observations are ordered by epoch. The original submitted optical
coordinates are retained when catalogue corrections are applied.

Radar measurements participate in the fit, but the initial-orbit observation
pool is built from optical observations. For controls and supported files, see
the [Observations guide](../../guide/observations/).

## Observer geometry

For a ground station, Trajecta uses the longitude and parallax constants from
the observatory catalogue. A supported roving observation can instead supply
its own terrestrial coordinates. Spacecraft observations use the supplied
geocentric observer offset.

For optical and radar stations, Trajecta uses the Earth PCK to transform the
terrestrial position and velocity from ITRF93 to J2000, then rotates them into
ecliptic axes. Earth's DE440 state is combined with the station offset. The
radar calculation evaluates the receiver at reception and the transmitter at
transmission.

When the PCK does not cover an epoch, the station calculation can use observed
C04 or historical UT1 data with the IAU 1976/1980 rotation model. The historical
path uses a nominal pole. Without suitable time data, the calculation reports
an error. Historical and predictive PCKs have different coverage and
accuracy; a prediction is not a measurement. See
[Dynamics and time](../dynamics/#earth-orientation) for the data sources.

## Optical prediction

Optical coordinates refer to the direction from the observer at reception to
the object at emission. Trajecta solves the one-way light-time equation
iteratively. The calculation accounts for the Sun's displacement between
emission and reception when forming the line of sight from heliocentric states.

The fitted optical direction also includes solar gravitational light deflection.
Its treatment depends on the observation's astrometric convention:

- **Catalogue-referenced:** the model accounts for the finite-distance
  deflection of the object relative to the deflection of reference stars.
- **Absolute:** the model retains the object's solar deflection without
  subtracting the reference-star term.
- **Gaia FPR:** an absolute convention is used, with the photocentre treatment
  described below.

These are astrometric directions. The optical fit does not apply a separate
stellar-aberration correction to them. They should not be confused with an
apparent direction in a reference system that includes stellar aberration.
The ephemeris feature has its own output options and conventions; see
[Ephemerides and close approaches](../../ephemerides-and-close-approaches/).

## Stellar occultations

An occultation record supplies an astrometric position derived from an
asteroid passing in front of a star. Trajecta uses that position in the optical
fit; it does not reconstruct the position from disappearance and reappearance
times or fit the asteroid's shape to individual chords.

The MPC 80-column reader recognizes occultation observations by technique
code `E`. ADES XML and PSV records can provide either an RA/Dec position or
the occulted star's `raStar` and `decStar` together with `deltaRA` and
`deltaDec`. For the latter, the reader adds the declination offset and converts
`deltaRA`, expressed as ΔRA cos(Dec), to an RA offset. The offset components
are in arcseconds; the star coordinates are in degrees. ADES `dist`/`pa`
offsets are not supported. An ADES record must use either direct RA/Dec
coordinates or the star-and-offset representation; mixing the two is not
supported. See the [MPC field descriptions](https://docs.minorplanetcenter.net/mpc-ops-docs/data-and-services/schema/obs-sbn/).

The resulting direction is predicted using the observer geometry in the
record, one-way light time and the optical deflection convention described
above. When both supplied angular sigmas are valid, occultation observations
use them even if **Prefer observer-supplied sigmas** is off; an available valid
RA/Dec correlation is also used. Otherwise the ordinary weighting profile
applies. The separate timing-uncertainty term is not added for occultations.

## Radar prediction

Radar prediction models the transmitted and received signal separately.
It iterates the down-leg time from the object to the receiver and the up-leg
time from the transmitter to the object. The transmitter and receiver can be
different stations. The object is propagated from its reception-epoch state
to the bounce epoch with a short two-body propagation.

The computed delay includes the light-travel time and the implemented solar
Shapiro and tropospheric corrections, as well as a clock-scale correction.
The Doppler model uses the transmitted frequency, line-of-sight motion and
its implemented relativistic and delay-rate corrections. A surface return
uses the supplied object radius; if it is absent, preparation uses a **1 km**
fallback and logs a warning. A centre return does not use this radius.

Radar delay is represented in seconds and Doppler in hertz. Their residuals
are reported with the units shown in the observation table. Radar prediction
requires valid station geometry, ephemerides and Earth-orientation data.

## Residuals

Residuals use **O−C**: the observed measurement, after applicable corrections,
minus the computed prediction. For optical observations the two components are:

```text
Δα cos δ = wrap(α observed − α computed) · cos(δ submitted)
Δδ       = δ observed − δ computed
```

The wrap keeps the right-ascension difference on the appropriate side of the
0°/360° boundary. When catalogue debiasing is enabled, the residual uses the
corrected coordinates; the observation table displays the submitted RA and
Dec. Radar delay and Doppler residuals also use O−C.

The fit's optical RMS is computed from the two angular residual components of
included optical observations. It is separate from the normalized residuals
used with measurement uncertainties.

## Catalogue debias

Catalogue debiasing adjusts optical coordinates for estimated catalogue-dependent
systematic offsets. Trajecta looks up a correction using the submitted position,
catalogue code and observation epoch, then subtracts it from the submitted
coordinates. Radar measurements are not debiased this way.

The current fit preparation selects the EFCC18 table. If the table is absent
or cannot be loaded, no catalogue correction is applied; observations without
an applicable correction also remain unchanged. The applied correction can be
inspected in observation details. See [Data sources](../../data-sources/) for
the table's availability.

## Measurement weights

Optical weighting uses the selected profile:

| Profile | Adopted uncertainty |
| --- | --- |
| `vfcc17.rules` | Looks up RA and Dec sigmas using station, technique, catalogue, program code and date. |
| Uniform | Uses 0.5 arcsec for each optical coordinate. |

When VFCC17 rules are unavailable or no rule resolves an observation, the
fallback is **1 arcsec per coordinate**. Observer-supplied optical RA and Dec
sigmas can take precedence when enabled and valid. Gaia FPR and occultation
observations have special handling for supplied sigmas; invalid or unsupported
values are not used as ordinary observer sigmas. A valid supplied correlation
can describe the orientation of the optical error ellipse.

With VFCC17, eligible closely spaced observations from the same station can be
downweighted as a group. For a batch of more than four observations, the
residual multiplier is `sqrt(4 / N)` on eligible coordinates. Observer-supplied
sigmas and explicitly forced coordinate weights are not adjusted by this
step. This is a weighting adjustment, not removal of observations.

Radar delay and Doppler use the uncertainties supplied for their respective
measurement channels.

## Timing uncertainty

When timing uncertainty is enabled, an optical observation's time uncertainty
contributes along its computed sky-motion direction to the error ellipse.
Trajecta uses a supplied timing RMS when available, otherwise stated time
accuracy, otherwise a configured default. Video observations use the
implemented 10 ms timing value. The timing term is not applied to spacecraft,
occultation or Gaia station 258 observations in this preparation path.

This modifies the residual's uncertainty model; it does not change the reported
observation epoch. The result depends on both the timing uncertainty and the
object's apparent angular motion.

## Normalized residual statistics

The observation table shows **χ (meas.)**; the details panel can also show
**χ (post-fit)**. Both are square roots of their respective χ² statistics.
They describe individual observations, not the reduced chi-square of the
complete fit.

For a weighted optical observation, let `r` be its two-component angular
residual and let `W` be the transform that normalizes it using the adopted
measurement uncertainty. The transform can include correlation and timing
uncertainty. With `z = W r`:

```text
χ (meas.) = sqrt(zᵀ z)
```

The post-fit statistic also accounts for the fitted parameter covariance.
Let `J` be the weighted residual Jacobian for the observation and `C` the
formal fit-parameter covariance. Trajecta uses:

```text
Q = J C Jᵀ
S = I − Q  for an observation included in the fit
S = I + Q  for an observation excluded from the fit
χ (post-fit) = sqrt(zᵀ S⁻¹ z)
```

The post-fit statistic uses the linearized observation model and formal
parameter covariance consistent with the adopted weights. It does not account
for unmodelled systematic errors. The value is shown only when the calculation
produces a valid statistic. For radar, the measurement statistic combines the
available normalized delay and Doppler residuals; the post-fit statistic
similarly accounts for the fitted covariance.

## Observation selection

When enabled, sigma-clip selection uses post-fit residual statistics and can
refit after changing the observation mask. It uses separate rejection and
recovery thresholds and safeguards for the retained observations and arc.
Manually locked selections are respected. The automatic rejection step does
not reject radar observations. The accepted mask is part of the solution and
is reused by Observation MC.

## Gaia FPR

For suitable numbered asteroids, the application can request Gaia FPR
observations and re-solve. These are additional spacecraft observations with
their supplied observer offsets, sigmas and correlation. The current model
uses a per-observation approximation for transit systematics.

When a physical radius is available, the Gaia path includes a phase-dependent
photocentre correction and an additional uncertainty term. The correction is
not applied without a positive radius. The loader checks for potentially
overlapping existing Gaia data before adding FPR observations. Availability
and the data returned by the archive determine whether a re-solve can proceed.
