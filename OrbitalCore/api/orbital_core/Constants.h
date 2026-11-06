#pragma once

namespace Constants {
	static constexpr const double SECONDS_PER_DAY = 86400.0;
	static constexpr const double AU_PER_DAY = 173.1446326742408;  /// speed of light in AU/day
	static constexpr const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923;
	static constexpr const double TWOPI = (2.0 * PI);

	static constexpr const double DEG_TO_RAD = PI / 180.0;
	static constexpr const double RAD_TO_DEG = 180.0 / PI;

	static constexpr const double J2000_JD = 2451545.0;
	static constexpr const double JULIAN_CENTURY_DAYS = 36525.0;

	static constexpr const double GAUSS_K = 0.01720209895;  /// Gaussian gravitational constant (rad/day)
	static constexpr const double EARTH_RADIUS = 6371.0;    /// km
	// The unit the MPC parallax constants (rho cos phi', rho sin phi') are
	// expressed in, mirroring the reference's reau = eradkm / aukm.  A third
	// Earth radius with a third job: it scales observatory positions, so it
	// belongs to astrometry, not to geometry or to the gravity field.  See
	// EARTH_EQ_RADIUS_KM and EARTH_J2_REF_RADIUS_KM further down.
	//
	// Collapsing it onto the WGS84 radius would move every observatory by
	// 0.7 m - about 1e-6 arcsec at 1 AU, five orders under the residual
	// floor - and would remove a real split, since MPECParser builds ADES
	// observatory ECEF straight from WGS84 instead.  Left alone here because
	// that is a decision about parity with the reference, not a cleanup.
	static constexpr const double EARTH_PARALLAX_RADIUS_KM = 6378.1363;
	static constexpr const double EARTH_RADIUS_IN_METERS = 6371000.0;  /// m
	static constexpr const double EARTH_LOCAL_DYNAMICS_AU = 0.1;
	static constexpr const double EARTH_TARGET_PLANE_AU = 0.1;
	static constexpr const double EARTH_CLOSE_APPROACH_REPORT_AU = 0.2;
	static constexpr int ENCOUNTER_STEP_POINTS = 100;

	static constexpr const double MU_IN_KM = 1.3271244004127942e11;      /// km^3/s^2
	static constexpr const double MU_IN_AU = 3.964015993175041e-14;      /// AU^3/s^2
	static constexpr const double MU_IN_METERS = 1.3271244004127942e20;  /// m^3/s^2
	static constexpr const double AU_IN_METERS = 1.495978707e11;         /// m

	static constexpr const double MU_IN_UNIT = 1.0;
	static constexpr const double EARTH_INCL = ((23.0 * 3600.0 + 26.0 * 60.0 + 21.448) * PI /
	                                            (180.0 * 3600.0));  /// reference EQUM00 -> ECLM00 obliquity (rad)
	static constexpr const double GM_SUN = 2.9591220828411956e-4;   /// DE440 GMS, AU^3/day^2
	static constexpr const double ARCSEC_RAD = PI / (180.0 * 3600.0);

	static constexpr const double EARTH_J2 = 1.0826267e-3;
	static constexpr const double EARTH_MOON_MASS_RATIO = 81.300568221497215;  /// DE440 EMRAT

	static constexpr const double EARTH_EQ_RADIUS_KM = 6378.137;
	static constexpr const double EARTH_EQ_RADIUS_AU = EARTH_EQ_RADIUS_KM / (AU_IN_METERS / 1000.0);

	// WGS84 flattening, the companion of the radius above.
	static constexpr const double EARTH_INV_FLATTENING = 298.257223563;
	static constexpr const double EARTH_FLATTENING = 1.0 / EARTH_INV_FLATTENING;

	// The dynamical one: the reference radius the EGM96 / JGM-3 gravity field
	// is expanded about.  It is a matched pair with EARTH_J2 below (EGM96
	// gives 1.08262668e-3 at this radius), so the two must move together or
	// not at all - substituting the WGS84 radius here silently changes the
	// force model.  Use only for gravity-field terms.
	static constexpr const double EARTH_J2_REF_RADIUS_KM = 6378.1363;
	static constexpr const double EARTH_J2_REF_RADIUS_AU = EARTH_J2_REF_RADIUS_KM / (AU_IN_METERS / 1000.0);

	static constexpr const double EARTH_RADIUS_IN_AU = EARTH_PARALLAX_RADIUS_KM / (AU_IN_METERS / 1000.0);

	static constexpr const double LUNAR_DISTANCE_KM = 384399.0;
	static constexpr const double LUNAR_DISTANCE_AU = 2.569551e-3;  // 384399 / 149597870.7

	static constexpr const double GM_EARTH_PLUS_MOON = 8.9970113929473466e-10;  /// DE440 GMB, AU^3/day^2
	static constexpr const double GM_EARTH_ALONE = 8.8876924467071022e-10;      ///< GMB split by DE440 EMRAT.
	static constexpr const double SUN_EARTH_MASS_RATIO = GM_SUN / GM_EARTH_ALONE;
	static constexpr const double SUN_EARTH_MOON_MASS_RATIO = GM_SUN / GM_EARTH_PLUS_MOON;
}  // namespace Constants
