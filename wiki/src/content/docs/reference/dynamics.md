---
title: Dynamics and time
description: Coordinates, force models, numerical propagation and time scales.
---

Trajecta propagates an orbit by integrating its equations of motion with the
selected force model. The same dynamics is used by orbit fitting, ephemerides
and risk calculations, subject to the options of the particular calculation.

## State and reference frame

The orbital state consists of heliocentric position **r** and velocity **v**
in the J2000 ecliptic frame. Core propagation uses AU for position, AU/day for
velocity and AU/day² for acceleration. Results shown or exported by a feature
may use other units or an equatorial frame; check the accompanying labels.

With an available DE440 source, `BarycentricIntegrator` translates the initial
heliocentric state to the Solar System barycentre, integrates it there and
translates the resulting state back. The Sun's barycentric position, velocity
and acceleration come from the planetary ephemeris. This change of origin does
not change the reference-plane orientation.

## Equations of motion

The propagated state follows `dr/dt = v` and `dv/dt = a`. In the heliocentric
force evaluation, the acceleration combines:

| Contribution | Description |
| --- | --- |
| Sun | Central Newtonian solar gravity. |
| Planets, Moon and Pluto | Direct point-mass accelerations for the enabled bodies. Earth and Moon are separate perturbers. |
| Solar System origin | Correction for the barycentric acceleration of the Sun supplied by DE440. |
| Massive asteroids | Optional direct perturbations from the SB441-N16 kernel. |
| Relativity | Einstein–Infeld–Hoffmann (EIH) correction for the Sun and enabled planetary bodies. |
| Solar oblateness | Solar J2, with a fixed tilted pole in the J2000 frame. |
| Earth oblateness | Earth J2, J3 and J4 within 0.1 AU of Earth. |
| Non-gravitational | An optional acceleration defined by the selected model and parameters. |

DE440 supplies body states and gravitational parameters. If the selected
calculation enables massive asteroids, SB441-N16 must also be available. Its
asteroid point-mass terms are separate from the EIH correction. A body's
selection in **Perturbations** controls its direct contribution; the Sun's
DE440 origin acceleration is still used for the heliocentric transformation.

The Earth zonal terms are applied only while Earth is enabled and the object
is within 0.1 AU of its centre. Their axis follows the modelled Earth pole at
the calculation epoch. Solar J2 uses the DE440 J2 and reference-radius pair.

For an object identified as one of the perturbers, the corresponding direct
self-perturbation can be excluded through its target identifier. Compare force
settings and kernels when comparing trajectories produced by different tools.

## Non-gravitational acceleration

The model resolves acceleration into radial, transverse and orbit-normal
directions. With unit vectors **R**, **T** and **N**, it has the form:

```text
a_ng = g(r) · (A1 R + A2 T + A3 N)
```

`A1`, `A2` and `A3` are acceleration coefficients in AU/day². The radial
direction points away from the Sun, the normal follows `r × v`, and the
transverse direction completes the orbital frame.

| Model | Distance law and parameters |
| --- | --- |
| Inverse square | `g(r) = 1 / r²`, with `r` in AU; supports A1, A2, A3 and an optional time shift DT. |
| Yarkovsky | The same inverse-square law with a transverse A2 term. |
| Marsden H2O | Water-sublimation law with A1, A2, A3 and an optional time shift DT. |
| Marsden CO | Carbon-monoxide law with A1, A2, A3 and an optional time shift DT. |

The Marsden laws use `g(r) = α (r/r0)^(-m) [1 + (r/r0)^n]^(-k)`.
For a nonzero DT in a supported model, the implementation evaluates the
distance for `g` at a shifted point on the instantaneous two-body orbit. The
acceleration directions are evaluated from the current state. Which coefficients
are estimated is part of the fit configuration; a fitted coefficient is not
automatically a direct measurement of a single physical mechanism.

## Numerical integration

OrbitalCore provides the fifteenth-order Gauss–Radau integrators **RA15** and
**IAS15**. The application fit settings expose RA15 and an initial step; the
integrator adapts its subsequent steps. Encounter monitoring can further
restrict steps near a target body. An initial step is not a fixed output interval.

The current core constants are `RA15Integrator::kAccuracy = 1e-12` and
`IAS15Integrator::kEpsilon = 1e-9`. They control different step-selection
schemes and are not interchangeable error tolerances. IAS15's adaptive minimum
is 0.001 day; an encounter cap or a shortened final step can be smaller.
State and elapsed-time accumulation use compensated arithmetic.

When a calculation requests derivatives, the integrator advances the state,
state transition matrix and non-gravitational sensitivities together. Step
control uses the physical state rather than the size of the derivative blocks.
Dense output can evaluate a state between accepted integration steps without
making the output spacing the integration step size.

## State transition and covariance

The state transition matrix **Φ** describes how a change in the initial state
affects the propagated state. With `x = (r, v)` and `f = (v, a)`, the
variational equations are:

```text
dΦ/dt = (∂f/∂x) Φ,    Φ(t0) = I
```

The force pipeline supplies position and velocity derivatives of its
acceleration. When non-gravitational coefficients are fitted, their
sensitivities are propagated as additional columns. A covariance can then be
transported through the corresponding linearized state and parameter mapping.
The matrix is a local description around the nominal trajectory; a strongly
nonlinear uncertainty region may require sampling instead.

## Time scales

The core uses **TDB** for dynamical epochs, propagation and ephemeris queries.
MJD is a numerical representation of an epoch, not a time scale: an MJD must
still be identified as TDB, TT, UTC or UT1. Internal propagation paths use MJD
TDB; API or export formats may also express epochs as Julian dates.

| Scale | Role |
| --- | --- |
| UTC | Civil observation and interface times, including leap seconds. |
| TT | Terrestrial uniform time used in time and Earth-orientation transformations. |
| TDB | Dynamical time used for orbit integration and DE/SPK ephemerides. |
| UT1 | Earth-rotation time used where terrestrial orientation is needed. |

`TimeScale` converts at input and output boundaries. UTC-to-TT uses the leap
second history and `TT = TAI + 32.184 s`; TT-to-TDB uses the core's periodic
approximation. UT1 conversion uses available Earth-orientation or historical
time data. The use of a future UTC date does not provide knowledge of future
leap seconds or future measured Earth orientation.

## Earth orientation

The default Earth-orientation profile uses C04/IAU 1980 data for pole and
nutation-related corrections. The Earth zonal model evaluates the pole at the
current epoch. Outside the available C04 records, the implementation uses its
extrapolation and correction model; those values are not measured EOP data.

For TT–UT1 conversion, Trajecta first uses `finals2000A.all` where applicable,
then the historical `ET-UT.dat` table. Other code paths can use a documented
fallback when orientation data are unavailable. The actual data coverage and
frame model matter when comparing topocentric results or close Earth encounters.
See [Data sources](../../data-sources/) for the files and their status.

## Related checks

The [integrator checks](../../verification/integrator-round-trip/) include
two-body and forward/backward calculations. The
[Holman](../../verification/dynamics-holman/),
[Parijskij and Ceres](../../verification/dynamics-ceres/) and
[2023 BU](../../verification/dynamics-2023bu/) pages compare selected dynamics
and variational calculations with independent implementations. Their settings
and results are specific to the cases described on those pages.
