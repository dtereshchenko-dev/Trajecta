---
title: 'Integrators: forward and backward'
description: Analytic and round-trip checks of RA15 and IAS15.
---

How closely do RA15 and IAS15 return to their starting point? To find out, we
propagated each orbit forward and then back to its initial epoch. The tables
show the differences in position and velocity between the returned and
initial states. We also compared a circular orbit with its analytic solution.

A round trip checks internal consistency, not the absolute accuracy of the
forward trajectory: errors made on the way out can cancel on the way back.

## Two-Body Checks

We start at `(1, 0, 0) AU` at MJD 60000 and use the two-body acceleration
`-k²r/|r|³`. For the circular orbit, the initial velocity is
`(0, k, 0) AU/day`. The analytic solution is evaluated at the actual
represented endpoint, after one period. Each return leg starts with a fresh
integrator.

| Check | RA15 position | IAS15 position | RA15 velocity | IAS15 velocity |
| --- | ---: | ---: | ---: | ---: |
| Circular, forward vs analytic | 0.111 mm | 0.185 mm | 2.06e-11 m/s | 3.56e-11 m/s |
| Circular, forward and back | 0.179 mm | 0.067 mm | 4.28e-11 m/s | 6.72e-12 m/s |
| Eccentric, forward and back | 18.38 mm | 12.60 mm | 2.79e-9 m/s | 1.93e-9 m/s |

For the eccentric orbit, the initial velocity is
`(0, 1.303840 k, 0) AU/day`. We propagated it for `100 × 2π/k` days,
or about 100 circular-orbit years, before returning to the start. This is not
100 periods of the eccentric orbit.

Both return distances are below the `1e-11 AU` (about 1.50 m) threshold
used by the existing IAS15 test. This is a test threshold, not an accuracy
guarantee. The submillimetre values are also close to the scale of
floating-point rounding.

The two-body checks call the integrators directly; they do not use
`BarycentricIntegrator`.

## Close-Approach Scenarios

For [2023 BU](../dynamics-2023bu/) and
[Apophis](../../numerical-investigations/apophis-2029/), we used the same
initial states as in the linked studies, but ran Trajecta's production force
model. It includes DE440, the Sun, planets, Moon, Pluto, the 16 SB441-N16
asteroids, relativity, solar J2, and Earth J2/J3/J4 with a rotating pole and
a 0.1 AU cutoff. The Apophis run uses the inverse-square A1/A2 parameters
from its study; non-gravitational forces are off for 2023 BU.

Both legs run through `NBodyPropagator` and `BarycentricIntegrator`.
The adapter converts the heliocentric input to a barycentric state at the
start of each leg and converts the result back at its target epoch. The
return starts a fresh propagator from the forward result, so it does not reuse
the forward step sequence.

Each row is a separate round trip with only the endpoint requested. Closure
does not necessarily worsen as the interval grows.

| Object | Forward interval | RA15 position | IAS15 position | RA15 velocity | IAS15 velocity |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2023 BU | 7 days | 3.88 mm | 0.49 mm | 6.90e-9 m/s | 8.78e-10 m/s |
| 2023 BU | 8 days | 3.52 mm | 1.07 mm | 6.42e-9 m/s | 1.90e-9 m/s |
| 2023 BU | 40 days | 1.25 mm | 3.03 mm | 2.08e-9 m/s | 5.30e-9 m/s |
| 2023 BU | 365 days | 6.05 mm | 4.98 mm | 9.06e-9 m/s | 7.81e-9 m/s |
| 2023 BU | 1000 days | 12.71 mm | 44.16 mm | 1.98e-8 m/s | 6.57e-8 m/s |
| Apophis | 102 days | 0.028 mm | 0.035 mm | 3.58e-12 m/s | 4.64e-12 m/s |
| Apophis | 103 days | 2.14 mm | 0.80 mm | 4.50e-10 m/s | 1.35e-10 m/s |
| Apophis | 110 days | 0.53 mm | 29.87 mm | 4.75e-11 m/s | 7.19e-9 m/s |
| Apophis | 365 days | 39.52 mm | 133.26 mm | 8.78e-9 m/s | 3.16e-8 m/s |
| Apophis | 1000 days | 1.272 m | 44.52 mm | 3.01e-7 m/s | 1.05e-8 m/s |

## Configuration And Reproduction

We used the existing Linux VM build of OrbitalCore and the default integrator
settings: RA15 `kAccuracy = 1e-12`, IAS15 `epsilon = 1e-9`, and an IAS15
minimum adaptive step of `0.001` days. The initial step is 10 days; encounter
caps and the final step can be smaller than the adaptive minimum.

The reference Trajecta checkout was at
[`0b20b5c6de852e4674964e2c4aa9cc349d3f0cf4`](https://github.com/dtereshchenko-dev/Trajecta/commit/0b20b5c6de852e4674964e2c4aa9cc349d3f0cf4).
The VM used GCC 13.3.0 and a Release build of OrbitalCore. The library's
SHA-256 was
`a396e0443dd72e1b767ff15781bda90b450c819c43db4d7031ae899ef0fb2742`.

- [C++ driver](../../verification-data/integrator-round-trip/round_trip.cpp)
  and [full numerical output](../../verification-data/integrator-round-trip/results.txt).
- Initial states: [2023 BU Horizons vectors](../../verification-data/dynamics-2023bu/horizons-2023bu-1000-days.json)
  and [Apophis Horizons vectors](../../verification-data/dynamics-apophis/horizons-apophis-1000-days.json).
  Horizons supplies the initial states.

To compile the driver on Linux, set `TRAJECTA` to the source root,
`CORE_BUILD` to the OrbitalCore build directory, and `EIGEN_INCLUDE` and
`JSON_INCLUDE` to the dependency include directories:

```sh
g++ -O3 -DNDEBUG -std=c++20 -mavx2 -mfma \
  -I"$TRAJECTA/OrbitalCore/src" -I"$TRAJECTA/OrbitalCore/api" \
  -I"$CORE_BUILD/generated" -I"$EIGEN_INCLUDE" -I"$JSON_INCLUDE" \
  "$TRAJECTA/wiki/public/verification-data/integrator-round-trip/round_trip.cpp" \
  -L"$CORE_BUILD" -Wl,-rpath,"$CORE_BUILD" -lOrbitalCore -o round_trip
```

Run the resulting executable from a working directory containing the
`earth_orientation` directory:

```sh
/path/to/round_trip /path/to/linux_p1550p2650.440 /path/to/sb441-n16.bsp \
  "$TRAJECTA/wiki/public/verification-data/dynamics-2023bu/horizons-2023bu-1000-days.json" \
  "$TRAJECTA/wiki/public/verification-data/dynamics-apophis/horizons-apophis-1000-days.json"
```

The SHA-256 hashes of the DE440 and SB441-N16 files used for this run are
`29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062` and
`919d612ce3c72a78fc7158f9120156542d0f21e6b8b052e4c1339c759747fd90`,
respectively. For a comparable rerun, the library and Earth-orientation
inputs also need to match.

The VM loaded 19,756 Earth-orientation records using the default IAU1980
profile. Its C04 file ends on 2026-01-05, before the Apophis interval, so
that run depends on extrapolation. The files present in `earth_orientation`
had these SHA-256 hashes:

| File | SHA-256 |
| --- | --- |
| `eopc04.62-now` | `40347f18a10eb72b15e4f84e647b5f1b7da4709534d935e2e1289697bdfa79df` |
| `eopc04.cnc` | `f93279874565fd3025418b153c1b60bcc6b6ed84e9eec4bf2c72045821d3865c` |
| `finals2000A.all` | `715ff5ce1d3ce6733e4461d8c813105208cb66e4fe1fef3882fb0ff848be6542` |
| `ET-UT.dat` | `d631a0b1ec8e270395f8b8aa7d7e6685b6711f01bed5ef54ae851ca6b799d12f` |

Not every file necessarily contributes to every force evaluation.
