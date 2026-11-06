#include "orbital_core/hazard/PhaClassifier.h"
#include <limits>
#include <cmath>
#include "orbital_core/hazard/Moid.h"
#include "orbital_core/Orbit.h"

namespace orbital_core::hazard {

	PhaCriteria classifyPha(const Orbit& orbit, double absoluteMag, const PlanetMoidTable& moids) noexcept
	{
		PhaCriteria c;
		c.absoluteMag = absoluteMag;
		c.semiMajorAxisAu = orbit.getSemiMajorAxis();
		c.perihelionAu = orbit.getType() == OrbitType::Parabolic ? orbit.getPerigee()
		                                                         : c.semiMajorAxisAu * (1.0 - orbit.getExcentricity());
		c.moidEarthAu = std::numeric_limits<double>::quiet_NaN();

		c.perihelionOk =
		    std::isfinite(c.perihelionAu) && c.perihelionAu > 0.0 && c.perihelionAu < kPhaPerihelionThresholdAu;

		c.absoluteMagOk = std::isfinite(absoluteMag) && absoluteMag < kPhaAbsoluteMagThreshold;

		if (!c.perihelionOk) {
			c.isPha = false;
			return c;
		}

		c.moidEarthAu = moids[Planet::Earth];

		c.moidOk = std::isfinite(c.moidEarthAu) && c.moidEarthAu < kPhaMoidThresholdAu;

		c.isPha = c.perihelionOk && c.moidOk && c.absoluteMagOk;
		return c;
	}
}  // namespace orbital_core::hazard
