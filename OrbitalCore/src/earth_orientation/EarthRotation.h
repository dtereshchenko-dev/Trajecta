#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/OrbitalCoreExport.h"

namespace orbital_core {
	struct EarthOrientationSample;
}

namespace EarthRotation {

	// Independent IAU 2000B nutation and IAU 2006 precession implementation.
	// The current observation and Earth J2 paths do not use this model.
	namespace Iau2000B {
		// Epochs are JD TT; angles are radians.
		ORBITAL_CORE_API double meanObliquity(double jd_tt);
		ORBITAL_CORE_API void nutationAngles(double jd_tt, double& dpsi, double& deps);
		// Leading term dpsi * cos(mean obliquity), without complementary terms.
		ORBITAL_CORE_API double equationOfEquinoxes(double jd_tt);
		// direction +1: mean to true of date; -1: inverse.
		ORBITAL_CORE_API void applyNutation(double jd_tt, Vec3D& v, int direction);
		// direction +1: J2000 to true of date; -1: inverse.
		ORBITAL_CORE_API void applyPrecessionNutation(double jd_tt, Vec3D& v, int direction);
	}  // namespace Iau2000B

	namespace Iau1976_1980 {

		// Legacy observer and body-fixed chain, without measured pole offsets.
		ORBITAL_CORE_API double gmstRadians(double jd_ut1);
		ORBITAL_CORE_API double gmstRadiansMjd(double mjdUt1);
		ORBITAL_CORE_API double equationOfEquinoxes(double jd_tt);
		ORBITAL_CORE_API void applyPrecessionNutation(double jd_tt, Vec3D& v, int direction);

		// Nominal IAU1976/1980 pole retained for reference comparisons.
		// Production Earth J2 uses the EOP-aware chain below.
		ORBITAL_CORE_API Vec3D dynamicalPoleEquatorial(double jdTt) noexcept;

	}  // namespace Iau1976_1980

	namespace WithEop {
		// Earth-fixed to equatorial, selecting the model from the typed EOP sample.
		ORBITAL_CORE_API Vec3D terrestrialToEquatorialMjd(double mjdTt, const Vec3D& bodyFixed,
		                                                  const orbital_core::EarthOrientationSample& eop) noexcept;
		// Full BF pole for Earth J2. The typed EOP offsets select IAU1976/1980 or
		// IAU2006/2000A; the latter returns GCRS/ICRF axes. Epoch is always MJD TT.
		ORBITAL_CORE_API Vec3D poleEquatorialMjd(double mjdTt,
		                                         const orbital_core::EarthOrientationSample& eop) noexcept;
	}  // namespace WithEop

}  // namespace EarthRotation
