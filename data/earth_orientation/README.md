# Earth orientation and time data

The Earth binary PCK files `earth_historic.bpc`, `earth_latest.bpc` and
`earth_predict.bpc` supply ITRF93 orientation for the Earth gravity pole,
terrestrial stations and Earth-fixed crossing coordinates. Their observed and predictive coverage and
accuracy differ. Source: [NAIF generic PCK kernels](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/).

The current CI and installer inputs use the following NAIF files, renamed for
use by Trajecta. The binary data are not modified. These large files are ignored
by Git; the CI workflows download them and check the hashes below.

| Local file | NAIF source file | SHA-256 |
| --- | --- | --- |
| `earth_historic.bpc` | `earth_620120_260806.bpc` | `652b02c0b3a5a283d393497f7f96e1fe4aa47ee6dcc035a089eb9e6b5df86643` |
| `earth_latest.bpc` | `earth_000101_270208_261112.bpc` | `eb72aef57bcedb247e282510a881a608fe74e4d9e7869aaa130474f5efe8dc54` |
| `earth_predict.bpc` | `earth_2026_260806_2126_predict.bpc` | `9f565ab5aaccd33afe54f3a1b20cb8b136636c31a0bd4c0986fac722b931c280` |

The desktop application checks for newer versions in Data sources and lets the
user update each file separately. Updates are staged and applied on restart in
the user data directory. On Windows and Linux, the installer places all three
PCKs there and preserves existing files. Consult the actual file coverage rather
than assuming a future epoch is covered.

On Windows and Linux, `eopc04.62-now` and `eopc04.cnc` are installed in the user Earth
orientation directory with the PCK files. These C04 data supply the gravity-pole
fallback when PCK does not cover the epoch. C04 is also used for TT-UT1
conversion within its observed coverage and as a station-orientation fallback
within that coverage. The gravity fallback applies the consistency corrections
in `eopc04.cnc` and extrapolates outside the observed C04 span. Explicit time conversion
instead interpolates the uncorrected observed TT-UT1 values. The supplied C04
series is used from 1972-01-19 through 2023-01-15; earlier records in the file
are not loaded.

The frozen C04 inputs have these SHA-256 hashes:

| File | SHA-256 |
| --- | --- |
| `eopc04.62-now` | `91c60678740044b6f3b32395b63bc886d94ad3b5592bd50929f031e9a9ebfb9b` |
| `eopc04.cnc` | `cad01734f7d75a35fd20cfa98ceee5450c42504a45f657ec59f2af4884f71051` |

The gravity-pole interpolation uses 15th-degree Chebyshev windows with
overlap blending and endpoint extensions. The IAU1980 profile preserves the
native 80-column input convention; FullC04 reads the full nutation column.
The consistency trend is not a physical prediction of future polar motion.
Changing these files changes the force model and requires renewed comparison.

`ET-UT.dat` is installed in the same user directory. It supplies historical TT-UT1 values and supports the historical
orientation fallback with a nominal pole. Explicit UT1 conversions fail
outside the coverage of both sources. The Earth PCK supplies an orientation,
not an independent conversion between UT1 and TT.
