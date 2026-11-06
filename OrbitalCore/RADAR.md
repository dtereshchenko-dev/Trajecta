# Radar Observations

## Contract

`ObsNode::radar` and `ParsedObservation::radar` preserve one radar record with
one or two independent scalar measurements. Receiver is the observation station
code; transmitter is explicit. A missing delay or Doppler channel is absent,
not a zero measurement. RA/Dec are not defined for radar records.

- Epoch: reception UTC, converted to the existing internal TT convention.
- Delay and its sigma: round-trip seconds in the public Core contract.
- Doppler and its sigma: Hz, received minus transmitted frequency.
- Transmit frequency: Hz.
- ADES input: delay in seconds, rmsDelay in microseconds, frq in MHz.
- MPC R/r (Q/q) pairs: implied decimal fields converted on ingestion; paired
  identity, epoch and station fields must agree. Malformed radar fails loading.

`RadarReturn::Surface` is retained on import. Preparation requires an explicit
positive `surfaceRadiusKm`; it never silently treats a surface echo as a
center-of-mass measurement. Import does not infer a size from photometry.

## Observation Model

Radar stations use deployed `RADCODE.dat`, with cylindrical geocentric
coordinates in 1e-10 AU and longitude in degrees. These coordinates are distinct
from the rounded optical observatory catalogue. Unknown radar stations fail
preparation instead of falling back to optical coordinates.

The local file starts with `! RADCODE format=2`; every station is a single row:
`code jpl_code mpc_code longitude_deg rho_xy_1e-10AU rho_z_1e-10AU name`.
The last field can contain spaces. A missing verified JPL/MPC identifier is `-`.
`ObservatoryCatalog` builds lookup indices from these columns; there is no C++
correspondence table or separate alias record. Duplicate identifiers, conflicting
keys, unsupported versions, malformed coordinates and unknown stations are errors.
RWO v3 JPL identifiers use this shared lookup. This is our versioned format,
not a change to the legacy external catalogue format.

The canonical keys retain the legacy RADCODE convention, not a universal
MPC namespace: EISCAT is legacy radar key `273`, current MPC `259`, JPL `-73`.
The optical MPC `273` still denotes Euclid. Radar loading and preparation never
use the optical catalogue for station geometry. The file covers all 18 stations
in the JPL API response of 2026-10-18 and the documented DSS-25 (`-25`), with
current MPC identifiers verified by name and coordinates. All 18 positions now
use the API precision, replacing the older rounded values where present.
DSS25, SRT, Irbene and MeerKAT are absent from that response: their legacy
coordinates are retained and explicitly marked as not API-verified. DSS25's JPL
identifier is documented by the API; identifiers for the other three are not
invented. Coverage is of this astrometry catalogue, not every radar facility.
An external orbit solver still needs its own station catalogue to know
any additional stations; extending our file does not modify that installation.

Maintenance uses the public, unfiltered radar-observation endpoint
`https://ssd-api.jpl.nasa.gov/sb_radar.api?coords=1&c-coords=1`, not an orbit solution.
Check its signature/version and the units before comparing `coords` by JPL code.
Verify new code assignments against the MPC catalogue by name and location;
never infer an MPC ID from a JPL integer. Retain stations absent from the response,
record the date/provenance, and run the parser/Apophis interop tests after updating
the file. Updating RADCODE requires deployment/restarting the app, not recompiling
the lookup logic. Runtime loading, fitting and import/export never call this API.

`obs::predictRadar` uses a two-way model: moving transmitter and
receiver, iterated bounce/transmit times, solar Shapiro delay, tropospheric delay,
TT/TDB clock correction and relativistic Doppler terms. Earth orientation uses
C04 and the full IAU1980 nutation series on this path. The short light-time arc
uses two-body propagation; the fitted trajectory uses the requested force model.
The source must provide Earth and barycentric Sun states. This is not a shape,
spin, ionosphere or plasma model.

Local state derivatives use numerical differentiation of the same prediction.
They are chained through the propagated state/STM and non-gravitational
sensitivities, including velocity dependence. Radar sigmas whiten the native
channels; an unweighted mixed-unit objective is rejected. Delay/Doppler errors
are currently treated as independent; no channel covariance is imported.

## Fit and Output

IOD and ATT-EXT use optical observations. The final fit includes valid radar
records alongside the optical selection. Radar-only IOD is not implemented.
An existing solution can be refined with radar and changed non-gravitational
settings without repeating IOD. Stable observation IDs and manual overrides
apply to whole radar records, including both channels when present.

The objective and covariance degrees of freedom count scalar measurements:
two per optical record, one or two per radar record. Full orbital/non-grav
cross-covariance is retained. `FitResult::cost` remains optical unweighted RMS
in arcseconds; it must not be interpreted as a combined radar/optical RMS.
Radar residuals expose observed, computed, O-C, sigma and chi-square separately.

RWO v2 output includes optical records and radar range/range-rate records, with
separate column headers. RWO v2 and ESA v3 are accepted on import. Observer
positions use standard S/S and S/s pairs; versioned comment metadata preserves
original optical techniques and precision, radar frequency and native values.
Other readers can ignore those comments. A native RWO without transmitter
frequency uses a marked reference frequency for unit conversion, and the UI
shows range rate rather than an invented physical Doppler frequency. Radar CSV
output is removed. Imported fit selection is archived, not a manual override
for a subsequent calculation.

Observational MC perturbs each present radar channel by its supplied sigma and
the requested noise multiplier; the optical-only RNG draw order is unchanged.
The user-facing covariance MC path does not perturb observations.

## RADCODE v2 Verification, 2026-10-18

Evidence: `../build/apophis-export-20261018/radcode-v2/`.
`stations-audit.json` records the live API signature/coordinates, exact agreement
for all 18 returned stations, MPC identity/position checks and the four retained
stations absent from the response. The file contains 22 canonical stations.

- 68 focused Core tests pass, including the versioned parser, malformed data,
  namespace collisions, station geometry and OEF/RWO round trips; three focused
  workspace UI tests also pass.
- Full Apophis Yarkovsky fit converges: 9381/9473 selected, identical mask to
  before the coordinate update. A2 changes from -2.9127422518199802e-14 to
  -2.9124488574064904e-14 AU/day^2 with the updated station coordinates.
- An independent unmodified RWO reader reads all 9473 exported RWO rows, all 50 radar channels
  and 1364 observer positions. Selection and radar station pairs are preserved;
  native observer-coordinate read error is below 7e-8 m (roundoff).
- Native CAR/EQU/KEP OEF readers preserve A2 and all 7x7 covariance entries,
  including orbital/non-grav cross terms. Maximum error normalized by the two
  parameter sigmas is 1.44e-16.
- Core imports the full reference RWO v2 (9469 rows) and ESA v3 RWO (9577 rows).
  Each contains 50 radar records and 1364 satellite/occultation positions;
  values, sigmas, station identities, geometry and archived selection were
  checked against each respective source, not inferred from matching counts.
- The VM's original RADCODE and both native reader/fit binary hashes are
  unchanged. Only test inputs/outputs in its work directory were written.

These are import/export and station-catalogue checks, not a new comparison of
the two fit engines with identical updated station coordinates. The reference
solver installation deliberately retains its original rounded catalogue.

## Verification, 2026-10-15

See `../build/radar-20261015/` for frozen inputs, binaries, commands and outputs.

- Apophis: 9533 optical records plus 5 radar records / 7 radar channels;
  fitted A2 = -2.903000e-14 AU/day^2, sigma = 2.17214e-15.
- Independent live reference fit: A2 = -2.92594874114067e-14, sigma = 2.16741e-15.
  Difference is about 0.11 of our formal sigma. Rejection is independent, so
  this is not a byte-exact full-fit mirror.
- Ten identical-state monostatic/bistatic cases agree within 3.34e-11 s in
  delay and 4.88e-4 Hz in Doppler. See the fixture provenance below.
- 214 selected Core tests pass, covering parsers, native units, compact
  residuals, derivatives, non-grav refinement and manual radar inclusion.
- 40 workspace UI tests pass; opt-in native 3D and live refinement tests are
  skipped. Radar table/inspector/dialog captures cover dark and light themes.
- Optical 2020 VV, RA15 / 30 years: nominal, COV, RWO, residuals and all
  12661 hazard audit records plus four risk rows match the pre-change binary.
- Covariance MC, 2023 CX1: 32 clones, seed 42, 100 km boundary; result archive
  and all crossing IDs/states/times match the pre-change binary byte for byte.

This is not validation of all radar stations/epochs or a complete Core/Linux
test run. Observational-MC end-to-end radar distributions have not been validated.
Full numerical fixtures: `tests/data/radar_reference.md`.
