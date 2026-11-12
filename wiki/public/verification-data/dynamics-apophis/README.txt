Apophis dynamics comparison
===========================

The comparison uses the saved Horizons JPL #220 initial state and reference
vectors, Trajecta commit 0b20b5c6de852e4674964e2c4aa9cc349d3f0cf4,
ASSIST 1.2.3 and REBOUND 4.6.0. Calculations were run on Ubuntu 24.04 x86-64.

External kernels (not included):

* DE440: linux_p1550p2650.440
  SHA-256 29915576d0a6555766b99485ac3056ee415e86df4fce282611c31afb329ad062
* SB441-N16: sb441-n16s.bsp
  SHA-256 815981128c8ecf06547c4562a9b2ae76e4ae17b9a8831743afdd74ea3f823ecf

No DE441 kernel was downloaded or used locally. The saved Horizons response
identifies DE441 as its reference model.

Both drivers use A1 = 4.999999873689e-13, A2 = -2.901766720242e-14 AU/day^2
and A3 = 0. Trajecta uses RA15 and ASSIST uses IAS15, with default integrator
controls. Their Earth models differ: Trajecta uses its production rotating
pole, J2/J3/J4 and 0.1 AU cutoff; ASSIST uses its own Earth harmonics with J3
and J4 set to zero. The Trajecta run loaded C04 data with the last observed
record dated 2026-01-05; the 2029 sample is extrapolated.

Place the C04 files in the earth_orientation directory under the driver's
working directory. This run used:

  eopc04.62-now  SHA-256 40347f18a10eb72b15e4f84e647b5f1b7da4709534d935e2e1289697bdfa79df
  eopc04.cnc     SHA-256 f93279874565fd3025418b153c1b60bcc6b6ed84e9eec4bf2c72045821d3865c

For the test, Trajecta's EihRelativity.h was replaced with the provided
EihRelativity_sun_outer.h and the existing OrbitalCore build was rebuilt.
This test-only header restricts the outer EIH sum to the Sun while retaining
planetary terms in the inner sums. It is not a production configuration.
The driver dynamics_apophis_benchmark_model.cpp was linked to that core
without replacing the production Earth force for the first run. Two further
runs replaced EarthJ2.cpp with the provided J2-only rotating-pole and
fixed-pole implementations, respectively, and rebuilt the same core. Both
retain the 0.1 AU cutoff and use the same EIH and integrator settings. The
rotating-pole output includes all reference epochs; the other output excerpts
contain the three reported epochs.

ASSIST's driver uses its public gr_eih_sources = 1 setting without changing
ASSIST code. It defaults to no Earth harmonics, interpolated outputs, epsilon
1e-9 and a 0.001-day minimum step. Set DE440_FILE and SB441_FILE to the
existing kernels, and run the 1000-day comparison with:

  python assist_apophis_solar_eih.py horizons-apophis-1000-days.json

Use --earth-harmonics j2 or j234 to select the fixed-pole Earth terms.
The earlier exact-output runs used --sampling exact --min-dt 0.

The Horizons reference was generated on 2026-11-10. External kernels and
executable builds are not included.
