---
title: 'Dynamics: 2023 BU'
description: A 1000-day propagation comparison between Trajecta and ASSIST.
---

This test compares Trajecta and ASSIST propagation of **2023 BU** through
its close Earth approach in January 2023. Both calculations start from the
same state and use the same ephemeris files and selected force components.
The propagation spans **1000 days**, from 2023-01-20 to 2025-10-16, in TDB.

This is a consistency check between programs under the settings described
below. It does not establish an absolute error bound for close encounters.

## Initial State And Model

The initial heliocentric state is from JPL Horizons solution **JPL #23**.
Horizons supplies the initial conditions, not the reference force model.
Both calculations use DE440, the same SB441-N16 kernel, the Sun, eight
planets, the Moon and Pluto. Newtonian gravity, solar J2 and EIH relativity
are enabled, with the outer EIH sum restricted to the Sun. Non-gravitational
forces are disabled.

Trajecta uses its default RA15 and IAS15 settings and the barycentric
integrator adapter. ASSIST 1.2.3 uses REBOUND 4.6.0 IAS15 with the GLOBAL
step controller, epsilon `1e-9` and an explicitly set minimum adaptive step
of `0.001` days. This is the setting used in the ASSIST Apophis example,
not the REBOUND 4.6.0 default: a new simulation reports `min_dt = 0`.
Both programs sample the trajectories daily, using interpolation between
integration steps where applicable. A final shortened step may be smaller
than the minimum adaptive step.

Trajecta's solar-only outer EIH and the Earth-harmonic configurations below
are test-only changes; ASSIST source was not changed.

## Earth Harmonics Disabled

The table gives Euclidean differences in the heliocentric J2000 ecliptic
frame at the same TDB epochs. Earth's non-central gravity is disabled in both
programs.

| Epoch, TDB | Elapsed | RA15 vs ASSIST | IAS15 vs ASSIST | RA15 vs IAS15 |
| --- | ---: | ---: | ---: | ---: |
| 2023-01-27 00:00 | 7 days | 0.037 mm | 0.158 mm | 0.148 mm |
| 2023-01-28 00:00 | 8 days | 1.76 mm | 1.71 mm | 0.605 mm |
| 2023-03-01 00:00 | 40 days | 0.0558 m | 0.0488 m | 0.0175 m |
| 2024-01-20 00:00 | 365 days | 0.576 m | 0.915 m | 0.345 m |
| 2025-10-16 00:00 | 1000 days | 1.118 m | 1.728 m | 0.613 m |

ASSIST took 867 accepted steps. The RA15/IAS15 column compares the two
Trajecta solutions.

## Earth Zonal Harmonics

Both programs use a fixed J2000 pole and matching Earth J2, J3, J4 and
reference-radius values. The J2-only comparison sets J3 and J4 to zero in
both programs; the J2/J3/J4 comparison enables all three. Other selected
forces and initial conditions remain as above.

ASSIST applies these terms without a distance cutoff. The first column for
each model keeps Trajecta's usual **0.1 AU cutoff**; the second disables that
cutoff in Trajecta for the comparison. Configurations retaining Trajecta's
0.1 AU cutoff therefore do not have identical Earth-harmonic force models.
The fixed pole and no-cutoff variants are test configurations, not changes
to production defaults.

The values are direct Trajecta RA15 minus ASSIST position-vector distances:

| Epoch, TDB | J2, 0.1 AU | J2, no cutoff | J2/J3/J4, 0.1 AU | J2/J3/J4, no cutoff |
| --- | ---: | ---: | ---: | ---: |
| 2023-01-27 00:00 | 0.037 mm | 0.037 mm | 0.083 mm | 0.083 mm |
| 2023-01-28 00:00 | 1.88 mm | 1.88 mm | 2.40 mm | 2.40 mm |
| 2023-03-01 00:00 | 0.0526 m | 0.0526 m | 0.0732 m | 0.0732 m |
| 2024-01-20 00:00 | 1.176 m | 1.152 m | 0.697 m | 0.673 m |
| 2025-10-16 00:00 | 2.230 m | 2.176 m | 1.359 m | 1.304 m |

At 1000 days the corresponding RA15/IAS15 differences within Trajecta are
0.836 m and 0.837 m for J2, and 0.417 m and 0.417 m for J2/J3/J4,
respectively. ASSIST took 891 and 930 accepted steps for J2 and J2/J3/J4.
Its Earth-harmonic implementation used here includes J2, J3 and J4, not J5.

## Numerical Sensitivity

In a separate ASSIST run with no minimum step and interpolated outputs, the
1000-day IAS15 position differences from Trajecta IAS15 were 17.69 m with
Earth harmonics off and 17.81 m with J2/J3/J4. With the explicit
`0.001`-day minimum, the corresponding differences are 1.73 m and 1.77 m.
The accepted step sequences differ between these runs.

Earlier exploratory runs forced an exact integration step at every output
epoch. In the no-harmonics case, changing the requested output schedule
changed the Trajecta–ASSIST distance at 1000 days from 28.34 m to 4.85 m.
Those runs used different output and minimum-step settings from the tables
above, so the difference cannot be attributed to the force model alone.

## Reproduction

[Trajecta driver](../../verification-data/dynamics-2023bu/dynamics_2023bu.cpp),
[ASSIST comparison script](../../verification-data/dynamics-2023bu/assist_2023bu.py)
and [Horizons vectors](../../verification-data/dynamics-2023bu/horizons-2023bu-1000-days.json).
The local calculations use DE440; the saved Horizons response identifies DE441.
The script expects two reference paths; passing the same 1000-day file twice
uses its records once. Its third argument selects `none`, `j2` or `j234`;
the optional fourth selects `interpolate` (default) or `exact`, and
the optional fifth specifies the minimum step in days (default `0.001`).
Set `DE440_FILE` and `SB441_FILE` to the local kernel paths. For the
no-harmonics ASSIST run, for example:

```sh
python assist_2023bu.py horizons-2023bu-1000-days.json \
  horizons-2023bu-1000-days.json none
```

The script prints ASSIST's position and velocity differences from the saved
Horizons vectors. To reproduce the interprogram distances in the tables,
compare its vector components with the Trajecta driver's output at matching
epochs; the script does not calculate those distances itself.

Trajecta's solar-only outer EIH and the disabled or fixed-pole Earth-harmonic
models require the test-only core configuration described above. The J2/J3/J4
variant uses J3 `-2.532410e-6` and J4 `-1.619898e-6`; the no-cutoff
variant additionally disables the Earth-harmonic distance cutoff. The
standard VM source and library were restored after the runs.
