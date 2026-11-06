#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/BPlane.h"
#include "orbital_core/hazard/Planet.h"

#include <Eigen/Core>

namespace orbital_core::hazard {

	struct OpikTargetPlanePoint {
		Eigen::Vector2d coordsKm{Eigen::Vector2d::Zero()};
		BPlaneBasis basis{};
		// Focused impact-disk radius divided by the physical body radius.
		// Earth's LOV boundary is a 6420 km sphere; focusing uses that same
		// radius, independently of the surface and J2 reference radii.
		double gamma{1.0};
		double vInfKmS{0.0};
		bool valid{false};
	};

	struct OpikTargetPlaneDifferential {
		OpikTargetPlanePoint point;
		Eigen::Matrix<double, 2, 6> dCoordsDGeoState{Eigen::Matrix<double, 2, 6>::Zero()};
		bool valid{false};
	};

	[[nodiscard]] ORBITAL_CORE_API OpikTargetPlanePoint computeOpikTargetPlanePoint(const Vec3D& rGeoKm,
	                                                                                const Vec3D& vGeoKmS,
	                                                                                Planet planet) noexcept;

	// Analytic first derivative of the reference tpcar_mtpcar + opik_tpcar
	// transformation. Columns are with respect to (r[km], v[km/s]); output
	// coordinates remain in km.
	[[nodiscard]] ORBITAL_CORE_API OpikTargetPlaneDifferential
	computeOpikTargetPlaneDifferential(const Vec3D& rGeoKm, const Vec3D& vGeoKmS, Planet planet) noexcept;

}  // namespace orbital_core::hazard
