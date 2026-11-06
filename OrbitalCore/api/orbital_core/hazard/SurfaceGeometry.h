#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"

namespace orbital_core::hazard {

	struct ORBITAL_CORE_API SurfaceGeodetic {
		double latitudeRad{};
		double longitudeRad{};
		double altitudeKm{};
		bool valid{false};
	};

	[[nodiscard]] ORBITAL_CORE_API SurfaceGeodetic geodeticOnReferenceEllipsoidKm(const Vec3D& equatorialKm,
	                                                                              double equatorialRadiusKm,
	                                                                              double flattening) noexcept;

	[[nodiscard]] ORBITAL_CORE_API double referenceSurfaceSignedAltitudeKm(const Vec3D& bodyFixedKm,
	                                                                       double equatorialRadiusKm,
	                                                                       double flattening) noexcept;
}  // namespace orbital_core::hazard
