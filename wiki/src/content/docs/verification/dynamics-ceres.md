---
title: 'Dynamics: Parijskij and Ceres'
description: Trajecta and ASSIST comparison through an asteroid encounter with Ceres.
---

This comparison propagates asteroid (5303) Parijskij through its 1996
encounter with Ceres. The initial state and sampling follow the
[ASSIST Ceres example](https://assist.readthedocs.io/en/latest/jupyter_examples/5303_Ceres/).
Trajecta and ASSIST use the same initial state and selected force components;
their position difference at the sampled encounter is **0.12 mm**, increasing
to **2.66 cm** after 3,653 days. These are interprogram differences, not
established absolute propagation errors.

## Setup

The initial state is the example's barycentric equatorial ICRF state at
JD `2449718.5` TDB (1995-01-01), in AU and AU/day. ASSIST receives that
state directly. Trajecta's public propagation interface uses heliocentric
ecliptic J2000 inputs and outputs; its `BarycentricIntegrator` converts the
state for barycentric integration. The ASSIST output is transformed to the
same heliocentric basis for comparison.

Both programs use the same **DE440** and **SB441-N16** files. The example's
prose mentions DE441, but its code loads DE440; this comparison follows the
code and does not use DE441. The selected components are the Sun, eight
planets, Moon, Pluto, 16 massive asteroids including Ceres, solar J2 and EIH
relativity. ASSIST's EIH outer sum is set to 11 sources for comparison with
Trajecta's full EIH model; this is not the example's default solar-only outer
sum. No non-gravitational acceleration is used.

Earth's J2/J3/J4 harmonics are disabled in both programs; Earth's point-mass
gravity remains enabled.

The run used Trajecta commit `0b20b5c6de852e4674964e2c4aa9cc349d3f0cf4`,
default RA15 settings, ASSIST 1.2.3 and REBOUND 4.6.0 with IAS15 `GLOBAL`,
epsilon `1e-9` and minimum step `0`. ASSIST used
`integrate_or_interpolate`. Trajecta's existing VM core was temporarily rebuilt
without Earth harmonics for this test; its normal source and library were
restored afterward. No ASSIST source was modified. The calculation ran on an
Ubuntu VM.

## Sampled Encounter

The example uses 3,653 equally spaced output points over 3,653 days. The
minimum below is selected from those points, not obtained by refining a
continuous closest-approach time. Ceres positions for both distance
calculations come from the same ASSIST-loaded SB441-N16 ephemeris.

| Quantity | Trajecta | ASSIST |
| --- | ---: | ---: |
| Sampled minimum epoch, JD TDB | 2450337.669496 | 2450337.669496 |
| Sampled Ceres distance | 0.005560704128626 AU | 0.005560704128627 AU |

Both select the same grid point, approximately **831,869.5 km** from Ceres,
consistent with the rounded values in the ASSIST example.

## Trajectory Comparison

| Quantity | Difference |
| --- | ---: |
| Position at the sampled encounter | 0.115 mm |
| Ceres distance at that same sample | 0.095 mm |
| Maximum position difference across all 3,653 samples | 2.69 cm |
| Position after 3,653 days | 2.66 cm |
| Velocity after 3,653 days | 1.01e-9 m/s |

As a separate Trajecta-only diagnostic, excluding Ceres from the asteroid
perturbers changed the position by **45.4 km** at the selected encounter
sample and **45,118 km** at the final sample, showing the effect of Ceres in
this scenario.

## Numerical Sensitivity

The following checks use the same initial state and force configurations to
explore sensitivity to the numerical settings.

### Trajecta RA15

The production RA15 accuracy constant was not changed. Instead, a driver used
the existing integrator's fixed-step callback, with the indicated initial step
and a fixed step after the first accepted step. The final step can be shortened
to reach the requested epoch. These runs requested the final epoch only; the
adaptive baseline used the 3,653-point output grid. Fixed-step corrector
exhaustion is not reported by this diagnostic and was not independently assessed.

| Fixed step after first accepted step | Final position difference from adaptive RA15 | Accepted steps |
| ---: | ---: | ---: |
| 8 days | 212.155 m | 457 |
| 4 days | 95.56 mm | 914 |
| 2 days | 26.21 mm | 1,827 |
| 1 day | 26.20 mm | 3,653 |

The 2-day and 1-day fixed-step positions differ by **0.017 mm**. The 1-day
result differs from the default ASSIST result by **0.51 mm**. The former is
comparable to the spacing of binary64 coordinates at AU-scale magnitudes;
it should not be read as micrometre-level integration accuracy.

The sensitivity of the adaptive calculation was also checked by varying
`RA15Integrator::kAccuracy` temporarily on the VM. This compile-time step
control parameter is not an absolute position tolerance; the initial step was
10 days and the ASSIST reference was unchanged. The production value remains
`1e-12`.

| RA15 `kAccuracy` | Final position difference from default ASSIST | Accepted steps |
| ---: | ---: | ---: |
| `1e-12` | 26.64 mm | 179 |
| `1e-13` | 0.62 mm | 233 |
| `1e-14` | 0.51 mm | 305 |
| `1e-15` | 0.49 mm | 397 |

The following are direct vector differences between the corresponding RA15
final states, not differences between their distances from ASSIST:

| RA15 `kAccuracy` pair | Final position difference | Final velocity difference |
| --- | ---: | ---: |
| `1e-13` vs `1e-14` | 0.116 mm | 4.51e-12 m/s |
| `1e-14` vs `1e-15` | 0.069 mm | 7.57e-13 m/s |
| `1e-13` vs `1e-15` | 0.149 mm | 5.26e-12 m/s |

At the final epoch, these runs differ by only a few representable coordinate
increments for this initial step and force model.

At `1e-12`, changing the initial step between 0.001 and 20 days gave final
differences from ASSIST between 15.15 and 42.06 mm. At `1e-13`, the same
initial-step values gave 0.41 to 1.19 mm. On the original sample grid the
default RA15 difference from ASSIST was 0.12 mm at the sampled Ceres encounter,
1.68 mm around day 1,000 and 14.62 mm around day 2,000. The observed growth
is not a sudden jump at the sampled encounter. The reduction with a stricter
parameter suggests sensitivity to the adaptive step sequence in this case.

As a separate internal control, excluding Ceres from Trajecta's perturbers
changed the final difference between the `1e-12` and `1e-13` runs from
26.06 mm to 3.54 mm. With `1e-12` versus `1e-14`, the corresponding values
were 26.17 mm and 4.06 mm. This control uses a different force model and
shows sensitivity to the encounter.

### Output-Schedule Control

To check whether the difference above was caused by the requested output
epochs, the adaptive and 1-day fixed-step calculations were each run twice:
once with the original 3,653-point grid and once requesting only the final
epoch. The initial state and test force configuration were unchanged.

| Comparison at day 3,653 | Position difference | Velocity difference |
| --- | ---: | ---: |
| Adaptive RA15: grid versus final-only | 0 m | 0 m/s |
| Fixed-step RA15: grid versus final-only | 0 m | 0 m/s |
| Adaptive versus fixed, both on grid | 26.197 mm | 9.49e-10 m/s |
| Adaptive versus fixed, both final-only | 26.197 mm | 9.49e-10 m/s |

Zero means the exported double-precision state components were identical in
these runs. The output schedule therefore does not explain the observed
26.197 mm difference in this experiment. The adaptive and fixed runs also
use different initial steps and step sequences; the cause of their difference
has not been established.

An earlier version of the fixed-step diagnostic bypassed the barycentric
adapter and therefore integrated heliocentric states. Its numerical results
are superseded by the barycentric comparisons above.

### Trajecta IAS15

The same case was run through Trajecta's barycentric adapter using IAS15
`GLOBAL`, initial step 10 days and the production minimum adaptive step of
0.001 day. ASSIST used `GLOBAL`, epsilon `1e-9` and minimum step 0. The
Trajecta epsilon values below were changed only for this test; the production
value is `1e-9` and was restored afterward.

| Trajecta IAS15 epsilon | Final position difference from default ASSIST | Accepted steps |
| ---: | ---: | ---: |
| `1e-9` | 1.50 mm | 235 |
| `1e-10` | 0.50 mm | 335 |
| `1e-11` | 1.04 mm | 473 |

With the default epsilon, the difference at the sampled Ceres encounter is
0.082 mm. Changing the initial step from 0.001 to 10 days gives final
differences from default ASSIST between 1.42 and 3.23 mm. Grid and final-only
outputs give identical exported double-precision states for the default run.
Changing epsilon does not produce a monotonic difference from the fixed
ASSIST reference.

### ASSIST IAS15

With the same interpolated output grid and minimum step `0`, epsilon was
changed using ASSIST's existing configuration. Differences are relative to
its `1e-9` run at the final epoch.

| IAS15 epsilon | Final position difference from `1e-9` | Accepted steps |
| ---: | ---: | ---: |
| `1e-9` | reference | 244 |
| `1e-10` | 1.268 mm | 349 |
| `1e-11` | 1.265 mm | 497 |

The `1e-10` and `1e-11` results differ from each other by **2.481 mm**. The
differences do not decrease monotonically against the default run.

## Reproduction And Limits

The [Trajecta driver](../../verification-data/dynamics-ceres/trajecta_parijskij.cpp)
and [ASSIST comparison script](../../verification-data/dynamics-ceres/compare_parijskij.py)
are provided. Compile the driver against the same OrbitalCore build with
compatible compiler and Eigen alignment flags. Set `DE440_FILE` and
`SB441_FILE`, then run:

```sh
python compare_parijskij.py /path/to/trajecta_parijskij --convergence --output-schedule-check --ias15-check --ra15-check
```

The Trajecta driver requires a core built without the Earth-harmonic force,
for example using this [test implementation](../../verification-data/dynamics-apophis/EarthJ2_disabled_probe.cpp).
It does not disable the force in a normal production build, which retains
Earth harmonics.

The DE440 SHA-256 was
`29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062`;
the SB441-N16 SHA-256 was
`919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90`.
The ephemeris files are not included.

The encounter distance here is sampled on a grid. A continuous
closest-approach search and a fit to observations would be separate checks.
