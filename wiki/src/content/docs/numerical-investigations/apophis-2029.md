---
title: 'Apophis 2029: numerical sensitivity'
description: A 1000-day propagation comparison between Trajecta and ASSIST.
---

This test compares Trajecta and ASSIST propagation of (99942) Apophis through
its April 2029 close approach to Earth. Both calculations start from the same
state and use the same ephemeris files and selected force components. The
propagation spans **1000 days**, from 2029-01-01 to 2031-09-28, in TDB.

## Summary

Trajecta and ASSIST agree closely through the Apophis encounter under the
selected force models: their position differences are 2.36 mm at day 103 and
5.25 cm at day 110. Interprogram and internal comparisons are sensitive to
numerical settings after the encounter. Convergence at 1000 days and the
absolute propagation error of either code have not been established.

The comparison uses the test force configuration described below.

## Initial State And Model

The initial heliocentric state is taken from JPL Horizons solution **JPL #220**.
Horizons supplies the initial conditions, not the reference model for the
comparison. Both programs use DE440, the same SB441-N16 kernel, the Sun, eight
planets, the Moon and Pluto. The enabled components are Newtonian gravity,
solar J2, and EIH relativity with the outer sum restricted to the Sun. The same
inverse-square non-gravitational parameters are used: A1 = 4.999999873689e-13
and A2 = -2.901766720242e-14 AU/day²; A3 = 0.

The first comparison disables Earth's non-central gravity in **both** programs.
The following comparison enables its zonal harmonics with the pole fixed to
J2000. Trajecta's solar-only outer EIH and the Earth-harmonic configurations
are test-only changes; ASSIST source was not changed.

Trajecta used RA15 at its default settings. ASSIST 1.2.3 used REBOUND 4.6.0
and IAS15 with the `GLOBAL` step controller, epsilon `1e-9` and a minimum step
of `0.001` days, explicitly set as in the ASSIST Apophis example. A new
REBOUND 4.6.0 simulation reports a default minimum step of `0` days.
Daily ASSIST outputs use
`integrate_or_interpolate`; Trajecta also uses dense output between accepted
steps. The reference vectors and ASSIST script are provided below.

## Direct Comparison

The table gives Euclidean differences between Trajecta and ASSIST, in the same
heliocentric J2000 ecliptic frame and at the same TDB epochs. These are direct
vector differences, not differences between two distances from Horizons.

| Epoch, TDB | Elapsed | Position difference | Velocity difference |
| --- | ---: | ---: | ---: |
| 2029-04-13 00:00 | 102 days | 0.065 mm | 9.85e-10 m/s |
| 2029-04-14 00:00 | 103 days | 2.36 mm | 8.47e-8 m/s |
| 2029-04-21 00:00 | 110 days | 5.25 cm | 8.13e-8 m/s |
| 2030-01-01 00:00 | 365 days | 3.96 m | 4.74e-7 m/s |
| 2031-09-28 00:00 | 1000 days | 10.65 m | 2.39e-6 m/s |

The same test was repeated with Trajecta IAS15 through the production
barycentric adapter. All values below are position-vector differences from
the same ASSIST run and use the same selected forces and output epochs.

| Elapsed | Trajecta RA15 vs ASSIST | Trajecta IAS15 vs ASSIST | RA15 vs IAS15 |
| ---: | ---: | ---: | ---: |
| 102 days | 0.065 mm | 0.077 mm | 0.017 mm |
| 103 days | 2.36 mm | 2.37 mm | 0.017 mm |
| 110 days | 5.25 cm | 5.27 cm | 0.43 mm |
| 365 days | 3.96 m | 4.00 m | 4.66 cm |
| 1000 days | 10.65 m | 10.82 m | 0.17 m |

At 1000 days, Trajecta's RA15 and IAS15 results differ by 0.17 m. The
remaining Trajecta–ASSIST difference has not been attributed to a specific
force or numerical effect.

## Numerical Sensitivity

The following experiments explore numerical settings and output scheduling.

### ASSIST Numerical Sensitivity

With interpolation and no Earth zonal harmonics, the 1000-day position
difference from Trajecta was 48.06 m at epsilon `1e-9` with no minimum step,
and 88.66 m at epsilon `1e-10` with no minimum step. Compare these with the
10.65 m above.

At epsilon `1e-9`, the 110-day ASSIST positions obtained with minimum steps
of 0.0005 and 0.002 days differed from the 0.001-day run by 2.1 mm and
84.4 mm respectively. Without a minimum step, the same 110-day run used
523,354 steps and differed by 93.5 mm from the 0.001-day run.
[jorbit's documentation](https://jorbit.readthedocs.io/en/latest/core/data/constants.html)
also reports very small IAS15 steps for the Apophis encounter without a
minimum step.

### RA15 Step-Control Sensitivity

For the no-Earth-harmonics model above, a Trajecta-only experiment varied
`RA15Integrator::kAccuracy`, keeping the initial state, forces, 10-day initial
step and output epochs unchanged. This is an internal compile-time constant,
not a user-facing setting; its normal value is `1e-12`.

RA15 uses this constant to propose an adaptive step from the highest-order
coefficient of its acceleration polynomial. Specifically, the proposed
step magnitude is `(kAccuracy / H)^(1/9)`, where `H` is derived from the
largest scaled absolute component of `B[7]` for the physical state, the current
step and the Radau weight. The constant is not a tolerance on position error or a
guaranteed bound on accumulated error. The values were varied temporarily for
this experiment; the normal default remains unchanged.

Position differences below are relative to the default `1e-12` run at 1000
days, not absolute errors.

| `RA15Integrator::kAccuracy` | Difference at 1000 days |
| ---: | ---: |
| `1e-10` | 0.049 m |
| `1e-11` | 0.416 m |
| `1e-13` | 0.170 m |
| `1e-14` | 0.040 m |
| `1e-15` | 0.398 m |
| `1e-16` | 0.224 m |

The default run used 48 accepted steps whose endpoints fell between days 100
and 105; the shortest was approximately 1,139 seconds. In that run, temporary
instrumentation recorded no corrector exits on stagnation across the full
integration. The parameter sweep is not monotonic, and the cause of the
variations remains open.

### Close-Encounter Sensitivity

Further tests indicate that the close encounter strongly amplifies tiny
differences in the propagated state. In the full test model, the `1e-11` and
`1e-12` RA15 runs differ in position by 0.017 mm at day 90 and 0.416 m at day
1000. Restarting both calculations from the same day-90 state still leaves a
0.105 m difference at day 1000, so the part of the discrepancy generated
after the restart cannot be attributed solely to the earlier integration.

As a separate finite-precision sensitivity experiment, increasing one
heliocentric position component at day 90 by one representable `double` step
(0.0166 mm) changed the day-1000 position by 0.191 m.
With planets and asteroids disabled, the `1e-11`, `1e-13` and `1e-14` runs
were within 0.42 mm of the `1e-12` run at day 1000; an Earth-only planetary
run showed differences up to 0.435 m. These runs use different physical models.

The experiments show why very small numerical differences can become visible
after this encounter, and that the effect is not restricted to the massive
asteroids or the non-gravitational force. They do not identify a specific
rounding operation or rule out truncation, interpolation or force-evaluation
effects.

### Sensitivity To Output Epochs

Before switching to interpolation, an exploratory comparison used exact
integration to each requested output epoch and no minimum step. A second run
kept the same initial state, forces and integrator settings, but requested only days 0–12,
90, 101–104, 110, 365, 730 and 1000. The Trajecta–ASSIST distance at day
1000 changed from 13.48 m to 68.52 m. At the common epochs, the two Trajecta
outputs were identical at the recorded precision; the two ASSIST outputs
differed by 0.41 mm on day 102, 8.34 mm on day 103 and 82.00 m on day 1000.

In that experiment the script called REBOUND `integrate` at each requested
epoch with its default exact-finish behavior. Changing the output schedule therefore changes
the integration step sequence. These measurements show sensitivity to that
schedule, amplified after the encounter.

## Earth Zonal Harmonics

This comparison uses the same interpolated outputs and ASSIST minimum step of
`0.001` days as the main table. Both programs use the same fixed J2000 pole
and matching Earth J2, J3, J4 and
reference-radius values. The J2-only comparison sets J3 and J4 to zero in both
programs; the J2/J3/J4 comparison enables all three. Other forces and initial
conditions remain as above.

ASSIST applies these terms without a distance cutoff. The first column for each
model keeps Trajecta's usual **0.1 AU cutoff**; the second disables that cutoff
in Trajecta for the comparison. The fixed pole and no-cutoff variants are test
configurations, not changes to the production defaults.

| Epoch, TDB | J2, 0.1 AU | J2, no cutoff | J2/J3/J4, 0.1 AU | J2/J3/J4, no cutoff |
| --- | ---: | ---: | ---: | ---: |
| 2029-04-13 00:00 | 1.055 mm | 0.130 mm | 1.061 mm | 0.088 mm |
| 2029-04-14 00:00 | 3.16 mm | 2.76 mm | 2.91 mm | 2.50 mm |
| 2029-04-21 00:00 | 7.54 cm | 6.14 cm | 7.14 cm | 5.71 cm |
| 2030-01-01 00:00 | 7.44 m | 5.00 m | 6.89 m | 4.42 m |
| 2031-09-28 00:00 | 21.79 m | 13.99 m | 19.93 m | 12.04 m |

Each value is the direct position-vector distance between Trajecta and ASSIST
at the same TDB epoch. At 1000 days the corresponding RA15/IAS15 differences
within Trajecta are 0.37 m and 0.05 m for J2, and 0.41 m and 0.07 m for
J2/J3/J4, respectively. ASSIST took 1044 accepted steps with J2 and 1019
with J2/J3/J4.

## Reproduction

[Trajecta C++ driver](../../verification-data/dynamics-apophis/dynamics_apophis_benchmark_model.cpp),
[ASSIST comparison script](../../verification-data/dynamics-apophis/assist_apophis_solar_eih.py),
and [Horizons vectors](../../verification-data/dynamics-apophis/horizons-apophis-1000-days.json).
The Trajecta driver uses the test-only core configuration described above;
the driver alone does not restrict EIH to the Sun or disable Earth harmonics.

For ASSIST, set `DE440_FILE` and `SB441_FILE` to the local kernel paths, then run:

```bash
python assist_apophis_solar_eih.py horizons-apophis-1000-days.json
```

The script defaults to the main comparison settings. Use `--help` for sampling,
minimum step, epsilon and Earth-harmonic options. For the earlier exact-output
experiment, use `--sampling exact --min-dt 0`. The script reports ASSIST's
vector differences from the saved Horizons values. The interprogram distances
above require comparing these components with the Trajecta driver's output at
the same epochs. The local calculations use
DE440; the saved Horizons response identifies DE441 and supplies initial
conditions and output epochs, not an independently matched force model.
