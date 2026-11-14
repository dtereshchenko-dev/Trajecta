# Radar Prediction Fixtures

Generated on 2026-10-15 using an independent two-way radar predictor, not our
prediction implementation. `RadarObservation.SameStatePredictionsMatchIndependentReference`
consumes these files; regenerate only through an independently built reference.

## Cases

`radar_cases.dat` begins with the case count (10). Each case contains:

1. reception MJD TT, state epoch MJD TDB, transmit frequency Hz;
2. heliocentric ECLM J2000 position AU and velocity AU/day at reception;
3. receiver and transmitter body-fixed Cartesian positions AU.

The five reception epochs are the actual radar epochs in
`data/examples/NEO/apophis.txt`. Each is tested with TX=RX=251 and with
TX=253 / RX=251, at 2380 MHz. The reception states are propagated from the
gravity-only Apophis fit, then supplied identically to both predictors.

`radar_expected.dat` has one row per case: one-based index, round-trip
delay seconds, Doppler Hz, followed by the native receiver ground position AU
and velocity AU/day in equatorial J2000. The current golden test asserts the
observables; derivative correctness is covered separately by finite-difference
tests.

## Reference Generation

The reference used isolated inputs with exact epochs and a shared DE440/SB441
ephemeris bridge. The harness verified SHA-256 hashes of the linked objects,
archives and ephemeris kernels before execution.

Local evidence and harnesses: `build/radar-20261015/`, including
`native-manifest.json`, `native-comparison.json`, `cases.dat`, `native.dat`
and `ours.dat`. Full Apophis reference fit is separate (`a2-fit2/` remotely,
`native-fit.log`, `fit.dat`, `radar-check.dat` locally). It restores the original
MPC fractional reception epochs after reading the legacy integer-second `.rad`
format. The earlier `a2-fit/` experiment had a driver dimension error and is not
reference evidence.

Observed maximum differences: 3.33387e-11 s and 0.000487477 Hz. Test tolerances
are 5e-10 s and 0.001 Hz. The implementations are not bit-identical; station
velocity differentiation differs. The tolerances remain far below the supplied
Apophis measurement sigmas (4-4.5 microseconds, 0.1-0.25 Hz).

## PCK comparison with GRSS

`RadarObservation.PckComparisonWithGrss` uses a geometric Apophis-like state at
MJD 53398.0306038786 TDB, DE440, and the body-fixed station coordinates above.
The same state was propagated with GRSS 4.5.9 on Ubuntu using the identical
Earth PCK files. It checks the receiver state and both monostatic and bistatic
delay and Doppler at 2380 MHz. The 1e-7 s / 0.005 Hz comparison limits allow
for differences in the two radar models, including atmospheric and relativistic
corrections. The test is skipped if PCK does not cover the epoch.
