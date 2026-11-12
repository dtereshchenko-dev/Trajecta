---
title: Verification map
description: Existing reference comparisons, internal tests and open verification work.
---

This page maps existing evidence for selected parts of Trajecta. A reference
comparison applies only to its documented inputs, conventions and tolerances;
it is not a claim that the whole application has a particular accuracy.

This is an inventory of tests and recorded results. Some tests require data
files and are skipped when those files are missing; check the output of the
build you are using.

## Reference Comparisons

| Component | Evidence and criterion | Scope and remaining work |
| --- | --- | --- |
| Body-fixed orientations | [CSPICE reference matrices](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/data/body_orientation.md): 56 matrices for seven planets and the Moon at seven epochs, with an elementwise test tolerance of `2e-10`. | Covers these transformations; Earth orientation and trajectory calculations require separate checks. |
| Equinoctial coordinates | [Independent Fortran reference](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/data/equinoctial_reference.md): five reference conversions with analytic Jacobians. | The nearly circular diagnostic still needs an independent analytic-Jacobian reference. |
| Radar observables | [Radar reference fixtures](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/data/radar_reference.md): ten matched-state two-way delay and Doppler cases. Recorded maximum differences were `3.34e-11 s` and `0.000488 Hz`; test tolerances are `5e-10 s` and `0.001 Hz`. | Covers predictions at the selected states; a complete independently fitted orbit and radar derivatives need separate checks. |
| Ephemerides | [`EphemerisTest.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/EphemerisTest.cpp) contains selected Horizons comparisons, including Chesley observer quantities and 2023 BU near its encounter. The Chesley sky-direction criterion is below `0.02 arcsec`. | The force models differ. Further comparisons should cover more stations, epochs, uncertainty quantities and apparent coordinates. |
| Dynamics and STM | [3666 Holman](../dynamics-holman/): matched initial state and selected force components compared with ASSIST over 10,000 days, with additional finite-difference checks. | One main-belt case, with Earth harmonics disabled in both programs. |
| Massive-asteroid encounter | [Parijskij and Ceres](../dynamics-ceres/): matched-initial-state comparison with ASSIST over 3,653 days, including the sampled encounter. | The minimum is grid-sampled; both distance calculations use the same Ceres ephemeris. |
| MOV | [2019 MO](../2019-mo/): comparison of sampled Trajecta and OrbFit distributions, including quantiles and marginal CDFs. | OrbFit was modified for diagnostic output; impact-probability calibration needs separate study. |
| Orbit and atmosphere | [2023 CX1](../2023-cx1/): comparisons with ESA and JPL orbit information and published meteorite recoveries. | Orbit comparisons and atmospheric recovery comparisons answer different questions. The atmospheric comparison is retrospective. |
| Atmospheric sensitivity | [2008 TC3](../2008-tc3/): fixed-entry scenarios and published recovery coordinates. | Shows dependence on the chosen material assumptions. |

## Internal And Analytic Checks

These checks cover selected calculations within the implementation.

| Component | Existing tests | What they establish |
| --- | --- | --- |
| RA15 | [`NBodyPropagatorTest.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/NBodyPropagatorTest.cpp), including circular-orbit and state/STM checks. | Analytic circular-orbit cases and internal consistency under the test tolerances. |
| IAS15 | [`IAS15IntegratorTest.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/IAS15IntegratorTest.cpp): circular orbits, invariants, round trips and RA15 comparisons. | Covers selected analytic cases and consistency between implementations. |
| RA15 and IAS15 round trips | [Forward/backward checks](../integrator-round-trip/) for two-body cases, 2023 BU and Apophis. | Measures position and velocity closure, not the absolute error of the forward trajectory. |
| Earth J2/J3/J4 | [`EarthJ2Test.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/EarthJ2Test.cpp): compares the implemented acceleration with a finite-difference gradient of the zonal potential and checks partials. | Covers selected states; the complete model and rotating-pole chain need further independent checks. |
| Non-gravitational forces | [`NonGravModelTest.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/NonGravModelTest.cpp) compares analytic partials with numerical differences. | Checks derivatives for selected models and states, not the physical fidelity of a fitted non-gravitational model. |
| STM and covariance | [`NBodyPropagatorTest.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/NBodyPropagatorTest.cpp) includes STM versus finite differences and propagation consistency checks. | Provides local numerical checks, not independent long-term uncertainty calibration. |
| Time and Earth orientation | [`EarthOrientationTest.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/EarthOrientationTest.cpp) and [`EphemerisTest.cpp`](https://github.com/dtereshchenko-dev/Trajecta/blob/master/OrbitalCore/tests/EphemerisTest.cpp) cover interpolation, time conversions and related edge cases. | More topocentric comparisons are needed across epochs and EOP models. |

## Numerical Investigations

[Apophis 2029](../../numerical-investigations/apophis-2029/) documents
Trajecta/ASSIST agreement around a close encounter and sensitivity after it.
Convergence at 1000 days and absolute errors have not been established.
[2023 BU](../dynamics-2023bu/) compares selected configurations and records
their sensitivity to output sampling and minimum step.

## Priorities

1. Compare individual force accelerations at identical states against independent implementations.
2. Expand analytic and high-precision integrator benchmarks, with explicit acceptance criteria.
3. Establish numerical convergence on matched force models, including close encounters.
4. Expand geocentric and topocentric ephemeris comparisons with matched definitions, time scales and Earth-orientation inputs.
5. Validate propagated uncertainty and rare-impact estimates against independent references and convergence studies.
