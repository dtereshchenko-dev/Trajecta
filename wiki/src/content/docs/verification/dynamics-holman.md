---
title: 'Dynamics and STM: 3666 Holman'
description: A matched-initial-state comparison of Trajecta and ASSIST.
---

This comparison uses the initial state from the [ASSIST Holman example](https://assist.readthedocs.io/en/latest/jupyter_examples/VariationalEquations/)
to compare a propagated state and its six-dimensional state transition matrix
(STM). It covers **10,000 days** from 2020-01-01 TDB.

## Setup

The example's initial state is barycentric, equatorial ICRF, in AU and AU/day,
at JD `2458849.5` TDB. ASSIST uses it directly; for Trajecta, the Sun's
barycentric position and velocity are subtracted and the result is rotated to
heliocentric ecliptic J2000. Both use the same DE440 and SB441-N16 files.
The comparison uses Trajecta commit `0b20b5c6de852e4674964e2c4aa9cc349d3f0cf4`,
ASSIST 1.2.3 and REBOUND 4.6.0 on an Ubuntu VM.

The selected components are the Sun, eight planets, Moon, Pluto, 16 massive
asteroids, solar J2 and EIH relativity. ASSIST's EIH outer sum uses 11 sources;
no non-gravitational acceleration is applied. Earth's J2/J3/J4 harmonics are
disabled in both programs; Earth's point-mass gravity remains enabled.

Trajecta uses its default RA15 settings; ASSIST uses IAS15 with `GLOBAL`,
epsilon `1e-9` and minimum step `0`. ASSIST outputs use
`integrate_or_interpolate`. Its six initial variational vectors are rotated
into the equatorial frame, and its final STM is rotated back to the same
heliocentric ecliptic basis used by Trajecta.

## State Comparison

The table gives Euclidean differences between the two programs. There is no
known exact trajectory for this comparison.

| Elapsed days | Position difference | Velocity difference |
| ---: | ---: | ---: |
| 1 | 0.067 mm | 7.57e-13 m/s |
| 100 | 0.068 mm | 1.67e-11 m/s |
| 1,000 | 0.470 mm | 2.86e-11 m/s |
| 5,000 | 5.19 mm | 1.91e-10 m/s |
| 10,000 | 8.17 mm | 2.59e-10 m/s |

The same comparison was repeated with Trajecta IAS15 using the production
barycentric adapter and its default settings. The propagation includes the
STM; the table reports the state difference from ASSIST.

| Elapsed days | Position difference | Velocity difference |
| ---: | ---: | ---: |
| 1 | 0.067 mm | 7.57e-13 m/s |
| 100 | 0.070 mm | 1.29e-11 m/s |
| 1,000 | 0.537 mm | 2.16e-11 m/s |
| 5,000 | 3.16 mm | 1.40e-10 m/s |
| 10,000 | 10.02 mm | 3.27e-10 m/s |

## STM Comparison

All six columns and all four 3×3 blocks are compared in the same state basis.
Since position and velocity derivatives have different units, the tables
report their differences as responses to specified initial perturbations. Each
value is the largest Euclidean response difference among the three individual
axes of that input type; the two columns may attain their maxima on different
axes.

**Position response** (`∂r/∂r₀` and `∂r/∂v₀`):

| Elapsed days | Initial position perturbation: `1e-8 AU` | Initial velocity perturbation: `1e-10 AU/day` |
| ---: | ---: | ---: |
| 100 | 1.71e-13 m | 2.19e-13 m |
| 1,000 | 1.99e-11 m | 6.46e-11 m |
| 5,000 | 5.72e-9 m | 1.68e-8 m |
| 10,000 | 2.62e-8 m | 8.11e-8 m |

**Velocity response** (`∂v/∂r₀` and `∂v/∂v₀`):

| Elapsed days | Initial position perturbation: `1e-8 AU` | Initial velocity perturbation: `1e-10 AU/day` |
| ---: | ---: | ---: |
| 100 | 2.95e-19 m/s | 3.09e-19 m/s |
| 1,000 | 1.45e-18 m/s | 4.31e-18 m/s |
| 5,000 | 2.56e-16 m/s | 7.52e-16 m/s |
| 10,000 | 7.86e-16 m/s | 2.41e-15 m/s |

These values compare linearized responses to the stated perturbations; they
are not uncertainties or absolute errors.

As a separate consistency check at day 10,000, each program's STM was compared
with central differences of its own independently propagated perturbed states.
For initial ecliptic `x` perturbations of `1e-7` and `1e-8 AU`, Trajecta's
position-response residuals were 2.13 mm and 0.784 mm; ASSIST's were 10.64 mm
and 2.13 mm. For initial ecliptic `vx` perturbations of `1e-9` and
`1e-10 AU/day`, the corresponding residuals were 0.471 mm and 0.551 mm for
Trajecta, and 0.060 mm and 12.89 mm for ASSIST. The checks cover two
directions and position responses, and include finite-difference, integration
and rounding effects.

## Reproduction

The [Trajecta driver](../../verification-data/dynamics-holman/trajecta_holman.cpp)
and [ASSIST comparison script](../../verification-data/dynamics-holman/compare_holman.py)
are provided. Compile the driver against the same OrbitalCore build and use
compatible compiler and Eigen alignment flags. Set `DE440_FILE` and `SB441_FILE`
to the local files, then run:

```sh
python compare_holman.py /path/to/trajecta_holman --check-fd
python compare_holman.py /path/to/trajecta_holman --ias15
```

For this comparison, the existing VM core was temporarily rebuilt with the
Earth-harmonic force disabled using the
[test implementation](../../verification-data/dynamics-apophis/EarthJ2_disabled_probe.cpp).
The driver itself does not disable that force: a normal production build
retains it. The original VM source and library were restored after the run.
No ASSIST source code was modified.

The DE440 SHA-256 was
`29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062`;
the SB441-N16 SHA-256 was
`919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90`.
Neither kernel is included.

This main-belt case uses the force components described above. Other objects,
force configurations and close encounters require their own checks.
