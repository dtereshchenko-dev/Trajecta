# Earth Orientation Data

All Earth-orientation files live here and deploy under `earth_orientation/`
beside each executable. The common `EarthOrientation` provider currently uses
the **IAU1980 (C04) pole profile** by default. `IAU2000` selects
the finals2000A pole explicitly, not a second service or an automatic fallback.
There is no environment-variable profile override.

The frozen C04 datasets are identified by these checksums:

| File | SHA-256 |
|---|---|
| `eopc04.62-now` | `91c60678740044b6f3b32395b63bc886d94ad3b5592bd50929f031e9a9ebfb9b` |
| `eopc04.cnc` | `cad01734f7d75a35fd20cfa98ceee5450c42504a45f657ec59f2af4884f71051` |

C04 supplies daily UTC epochs, polar motion in arcseconds, UT1-UTC in seconds,
and IAU-1980 dpsi/deps in arcseconds. It is not interchangeable with the dX/dY
columns of `finals2000A.all`. Both readers convert UTC epochs to TT MJD and interpolate
continuous TT-UT1, including across leap seconds.

`finals2000A.all` uses Bulletin A xp/yp (arcseconds), UT1-UTC (seconds), and
IAU2000A dX/dY (milliarcseconds). Its reader uses linear interpolation, accepts
observed and predicted records, and keeps missing celestial offsets separate
from missing rotation data. Outside coverage it returns an invalid sample with
zero observed pole corrections and UT1=UTC, never C04 data. Field definitions:
[USNO format description](https://maia.usno.navy.mil/ser7/readme.finals2000A).

The compatibility policy uses the native 15th-degree Chebyshev windows, three
unused intervals, three overlap intervals, C6 blending, and C2 end extensions.
After the extended range the observed pole/nutation terms vanish and TT-UT1
continues with the series' secular slope. `eopc04.cnc` is applied afterward.
Its historical IERS-1993 consistency trend is **not a physical prediction of
future polar motion**. Replacing a dataset changes the force model and requires
a new documented mirror baseline.

`EarthOrientationProfile::IAU1980` explicitly reproduces the native 80-column input
buffer; `FullC04` preserves the complete dEps column. Neither convention reads
future EOP from orbit-fit results. The selected convention, source path and
extrapolation policy are logged as `[EarthOrientation]`.

Missing or invalid C04 files produce an explicit warning and a nominal-pole
fallback; it does not reproduce measured pole corrections. No network download occurs at runtime.

The optical clock is independent of the dynamical pole profile.
`EarthOrientation::ttToUt1Mjd` uses the finals2000A TT-UT1 series with the
historical cnvtim-style iteration on daily epochs. Optical observers apply
this UT1 correction to GAST but do not apply polar motion. Observer diagnostics
report `ut1FromEop` separately from the available (not applied) xp/yp values.
Missing time data leaves the caller's UTC fallback unchanged; it never silently
substitutes the dynamical C04 clock and its consistency corrections.
