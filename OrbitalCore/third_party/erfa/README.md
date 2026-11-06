# ERFA Earth Orientation Subset

Pinned upstream: [liberfa/erfa v2.0.1](https://github.com/liberfa/erfa/releases/tag/v2.0.1),
commit `9915ba38c9365f8b0738269b8c2ac1fdd5f8dee3`.

This is an **unmodified 35-file transitive C source subset**, not the full ERFA
distribution. Entry points: `eraXys06a`, `eraC2ixys`, `eraEra00`, `eraSp00`,
`eraPom00`, `eraC2tcio`. `eraS06` is already included as a dependency.
Both upstream headers are retained unchanged; declarations outside this subset
are not implemented here. No external runtime or system installation is needed.

## Provenance and License

- `LICENSE`: upstream BSD-3-Clause terms and SOFA-derived attribution. Preserve
  this file and individual source notices; redistribute notices with binaries.
- `manifest.json`: full commit, dependency graph, SHA-256 of each vendored file.
- Downloaded commit archive SHA-256:
  `f6f3613d7d78a7cc9f2bbe5a859c8fe94e73feb5db2ca9688363f597ee985a40`.
- Local additions only: this README, manifest, `sources.txt`, `sources.cmake`.
- Upstream source and regression tests are available in the pinned commit archive.

## CMake Integration

`sources.cmake` only exposes lists; including it creates no targets. Enable C in
the enclosing project, then adapt the following at the intended ownership level:

```cmake
enable_language(C)
include("${CMAKE_CURRENT_SOURCE_DIR}/third_party/erfa/sources.cmake")
add_library(orbital_erfa STATIC ${ORBITAL_ERFA_SOURCES})
target_include_directories(orbital_erfa PUBLIC "${ORBITAL_ERFA_INCLUDE_DIR}")
set_target_properties(orbital_erfa PROPERTIES
    C_STANDARD 99 C_STANDARD_REQUIRED YES POSITION_INDEPENDENT_CODE YES)
if(MSVC)
    target_compile_options(orbital_erfa PRIVATE /fp:precise)
else()
    target_compile_options(orbital_erfa PRIVATE -fno-fast-math)
    target_link_libraries(orbital_erfa PUBLIC m)
endif()
target_link_libraries(OrbitalCore PRIVATE orbital_erfa)
```

The enclosing repository CMake files are deliberately not modified by this
sidecar. Compile `.c` as C, not by renaming it to C++. `erfa.h` handles C++ linkage.
The independent reference uses stricter FP flags, recorded in its local CMake.

## API Contract

- `eraXys06a(2400000.5, mjdTt, ...)` gives IAU 2006/2000A CIP X/Y and CIO `s` in
  GCRS, already including frame bias. Do not apply a second bias correction.
- Add modern celestial-pole offsets dX/dY **in radians**, not dpsi/deps.
- Pass TT to XYS and `eraSp00`; pass **UT1**, preferably split day/fraction, to
  `eraEra00`. Do not use TT or UTC as UT1 inside table coverage.
- `eraPom00` uses xp/yp in radians. `eraC2tcio` returns **GCRS -> ITRS**.
  The ITRS body-fixed z-axis expressed in GCRS is row 2, equivalently
  `transpose(rc2t) * (0,0,1)`. It differs from the CIP when PM is nonzero.
- For the SOFA cookbook 5.5 convention, keep the original model `s` after adding
  dX/dY. `eraS06` evaluates a series for `s + XY/2`; recomputing it with corrected
  X/Y is a distinguishable convention, not a mandatory missing correction.

References: [ERFA XYS](https://github.com/liberfa/erfa/blob/v2.0.1/src/xys06a.c),
[matrix convention](https://github.com/liberfa/erfa/blob/v2.0.1/src/c2tcio.c),
[SOFA Earth Attitude cookbook, sections 5.5-5.6](https://www.iausofa.org/s/sofa_pn_c.pdf),
[ERFA s06](https://github.com/liberfa/erfa/blob/v2.0.1/src/s06.c).
