# Third-Party Components

Trajecta uses the following third-party libraries and data:

## Qt Framework
- **Modules:** Quick, Quick Controls, Quick 3D, Charts, Widgets, Network,
  State Machine, and their dependencies.
- **Licenses:** Qt modules have different licensing terms. Qt Quick 3D and
  Qt Charts are available under GPL v3 or a commercial Qt license; other
  modules may be available under LGPL v3, GPL, or commercial terms. Consult
  the license and attribution information for the exact Qt version and
  modules distributed with the application.
- **Source:** https://www.qt.io and https://doc.qt.io/qt-6/licensing.html

## ERFA
- **License:** BSD-3-Clause
- **Source:** https://github.com/liberfa/erfa
- **Usage:** Selected routines for astronomical reference systems and Earth
  orientation. ERFA is derived, with permission, from the IAU SOFA library;
  SOFA is not responsible for the ERFA implementation.
- **Notice:** The bundled license is in `licenses/ERFA.txt` in Windows packages
  and `OrbitalCore/third_party/erfa/LICENSE` in the source tree.

## Eigen
- **License:** MPL 2.0
- **Website:** https://eigen.tuxfamily.org
- **Usage:** Linear algebra, matrix operations

## Fonts
- **Inter:** SIL Open Font License 1.1; copyright The Inter Project Authors.
  Source: https://github.com/rsms/inter
- **JetBrains Mono:** SIL Open Font License 1.1; copyright The JetBrains Mono
  Project Authors. Source: https://github.com/JetBrains/JetBrainsMono
- **Source Sans 3:** SIL Open Font License 1.1; copyright Adobe.
  Source: https://github.com/adobe-fonts/source-sans
- **Cascadia Mono:** SIL Open Font License 1.1; copyright Microsoft.
  Source: https://github.com/microsoft/cascadia-code
- Their full notices are included in `licenses/` in Windows packages and in
  `TrajectaUI/resources/fonts/` in the source tree.

## Other Libraries
- **nlohmann/json:** JSON parsing and serialization.
- **libcurl and OpenSSL:** HTTP and HTTPS communication, including Gaia and
  SBDB requests and weather and terrain data acquisition.
- **GDAL and PROJ:** Raster formats and coordinate reference support.
- **OpenJPEG, libpng, libtiff, libgeotiff, json-c, SQLite and compression
  libraries:** Deployed when required by the native dependency graph.
- **Plotly.js, Leaflet and Lucide:** Embedded HTML report resources.

The Windows package includes notices for these libraries in `licenses/`.
GDAL and PROJ data resources are included under `share/`. Qt itself may contain
additional third-party components; consult the Qt distribution's notices for
the exact modules and binaries shipped.

## JPL DE440 Planetary Ephemeris
- **License:** Public Domain
- **Source:** NASA Jet Propulsion Laboratory
- **Website:** https://ssd.jpl.nasa.gov/planets/eph_export.html
- **Usage:** High-precision planetary positions

## OpenMP
- **Usage:** Parallel computation in orbit optimization and other calculations.
- The OpenMP runtime and its redistribution terms depend on the compiler and
  toolchain used to build the application.

## Project Pluto jpl_eph (jpleph.cpp / jpl_int.h / jpleph.h)
- **License:** GPL v2 or later
- **Author:** Bill J. Gray, Project Pluto (translated from Fortran by P.A. Dybczynski)
- **Source:** https://www.projectpluto.com/jpl_eph.htm — https://github.com/Bill-Gray/jpl_eph
- **Usage:** Reading binary JPL DE440 / DE441 planetary ephemeris kernels.
- **Modifications:** File I/O replaced with caller-supplied memory-mapped
  view (`jpl_init_ephemeris_from_view` instead of `jpl_init_ephemeris`),
  removing per-record `fread`/`fseek` syscalls on the LM-Jacobian hot
  path.  The original GPL v2 header is preserved verbatim; a dated
  modification notice was added at the top of each touched file.
  Distributed under GPL v3 with the rest of Trajecta ("any
  later version" clause of the original).

## Astronomical Data and Services
- **SB441-N16:** Optional JPL small-body ephemeris used for massive-asteroid
  perturbations. Source: https://ssd.jpl.nasa.gov/planets/eph_export.html
- **Earth orientation and time:** `finals2000A.all` from the IERS Rapid Service
  / Prediction Center, C04 and its consistency data, and historical ET-UT data.
  See `data/earth_orientation/README.md` in the source tree for file sources,
  conventions, coverage and fallbacks. Optical timing uses UT1 data; the radar
  observation model also uses Earth-orientation data. The Earth zonal force
  uses the C04-based pole orientation when Earth-orientation data are available.
- **Astrometric catalog debias tables:** JPL-hosted catalogue corrections used
  when debiasing is enabled. Source: https://ssd.jpl.nasa.gov/ftp/ssd/debias/
- **Observatory data:** Optical station information from the MPC and a local
  radar catalogue derived from JPL radar coordinates with MPC station
  identifiers. Sources: https://minorplanetcenter.net/iau/lists/ObsCodes.html
  and https://ssd-api.jpl.nasa.gov/doc/sb_radar.html
- **Gaia FPR and JPL SBDB:** Optional online sources of astrometry and physical
  parameters. Sources: https://gea.esac.esa.int/archive/ and
  https://ssd-api.jpl.nasa.gov/doc/sbdb.html

These datasets and services have their own provenance and conditions; consult
the respective providers when redistributing or updating them.
