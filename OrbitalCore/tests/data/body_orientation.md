# Body-Fixed Orientation Reference

`BodyFixedFrameTest.cpp` checks planetary and lunar body-fixed transforms
against independent CSPICE matrices. The test fixture does not require Python,
network access or the downloaded reference kernels. Current Earth orientation
uses CSPICE and separate Earth PCK files.

## Models

- Earth uses binary PCK ITRF93 when covered. The IAU 1976/1980 UT1/GAST
  transformation without polar motion remains a fallback when time data cover
  the epoch. Earth is not part of the 56-matrix fixture described below.
- Mercury, Venus, Mars, Jupiter, Saturn, Uranus and Neptune use the complete
  pole and prime-meridian model from NAIF `pck00011.tpc`, including periodic
  terms. Coefficients are generated into `src/hazard/IauRotationData.h`.
- The Moon uses numerical PA librations from the application's existing DE440
  binary, followed by the fixed PA-to-ME transformation in
  `moon_de440_250416.tf`. The output frame is `MOON_ME_DE440_ME421`, not PA.
  Missing librations or an unsupported ephemeris version produce an invalid
  transform rather than silently selecting an approximate lunar orientation.

The public MJD entry points take TDB. Matrices map equatorial J2000 vectors to
body-fixed vectors. The ecliptic-vector wrapper remains available. Latitudes/longitudes use the
selected reference surface, with east-positive longitude.

## Reference Data

Download these official reference kernels into a temporary development folder:

- <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc>
- <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc>
- <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf>

Run `generate_body_orientation.py <kernel-directory>` in an environment with
SpiceyPy. It regenerates the production coefficient header and
`body_orientation_expected.json`. The JSON records source hashes and toolkit
version. No kernel downloads occur during tests.

The fixture has 56 matrices: seven planets plus the Moon at seven TDB epochs
from 1900 through 2100. The elementwise tolerance is 2e-10; the tests also check
orthogonality, determinant, retrograde rotation, unavailable data, concurrent
libration queries and preservation of the ephemeris position cache.

## Verification on 2026-10-20

These results describe the configuration used on that date, before the Earth
PCK transition. Evidence: `build/body-orientation-20261020/`.

- 100 selected Core tests passed, including all five `BodyFixedFrame` tests.
- Lunar MC control: changing libration angles changes longitude but leaves
  clone crossing times and planet-relative inertial positions exactly equal.
- 48 workspace UI checks passed; six opt-in tests were skipped.
- Three existing report JavaScript tests passed.
- Earth coordinate probe: 20 epoch/vector combinations are byte-identical
  before and after the change.
- Fresh 2020 VV / RA15 / 30-year preservation gate:
  `build/radius-covariance-20261020/body-orientation-before/` and
  `body-orientation-after/`. Nominal, all 17,647 audit records, 981 completed VA
  records and all 16 risk rows agree exactly.
- Core, CLI and UI were rebuilt with matching DLLs. RA15/NBody and orientation
  header dependencies are nonzero and tracked by Ninja.

This is not a new independent comparison, a force-model change, or a performance
benchmark. Full Core tests and Linux were not run. Planetary gravity harmonics
and non-Earth map imagery were not added. The MC interface still exposes Earth,
Moon and Mars; the orientation API supports all eight planets and the Moon.
Old MC archives retain their original coordinates and are labelled as legacy
when orientation provenance is absent. Regenerating HTML alone cannot correct
old lunar or planetary coordinates.
