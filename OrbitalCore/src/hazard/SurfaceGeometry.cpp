#include "orbital_core/hazard/SurfaceGeometry.h"

#include <algorithm>
#include <cmath>

namespace orbital_core::hazard {

	SurfaceGeodetic geodeticOnReferenceEllipsoidKm(const Vec3D& equatorialKm, double equatorialRadiusKm,
	                                               double flattening) noexcept
	{
		SurfaceGeodetic out;
		if (equatorialRadiusKm <= 0.0)
			return out;

		const double r = equatorialKm.magnitude();
		if (r < 1.0e-12)
			return out;

		out.longitudeRad = std::atan2(equatorialKm.y, equatorialKm.x);

		if (flattening <= 0.0) {
			out.latitudeRad = std::asin(std::clamp(equatorialKm.z / r, -1.0, 1.0));
			out.altitudeKm = r - equatorialRadiusKm;
			out.valid = true;
			return out;
		}

		const double e2 = flattening * (2.0 - flattening);
		const double xy = std::hypot(equatorialKm.x, equatorialKm.y);

		double t = 0.0;
		double previousAltitude = 0.0;
		double altitude = 0.0;
		double sinLat = 0.0;

		for (int iter = 0; iter < 40; ++iter) {
			const double zShifted = equatorialKm.z + t;
			const double d = std::hypot(xy, zShifted);
			if (d < 1.0e-12)
				return out;

			sinLat = std::clamp(zShifted / d, -1.0, 1.0);
			const double n = equatorialRadiusKm / std::sqrt(1.0 - e2 * sinLat * sinLat);
			t = n * e2 * sinLat;
			altitude = d - n;

			if (iter > 0 && std::abs(altitude - previousAltitude) < 1.0e-5)
				break;
			previousAltitude = altitude;
		}

		out.latitudeRad = std::asin(sinLat);
		out.altitudeKm = altitude;
		out.valid = true;
		return out;
	}

	double referenceSurfaceSignedAltitudeKm(const Vec3D& bodyFixedKm, double equatorialRadiusKm,
	                                        double flattening) noexcept
	{
		const SurfaceGeodetic g = geodeticOnReferenceEllipsoidKm(bodyFixedKm, equatorialRadiusKm, flattening);
		return g.valid ? g.altitudeKm : -equatorialRadiusKm;
	}
}  // namespace orbital_core::hazard
