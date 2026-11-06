#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/BPlane.h"
#include "orbital_core/hazard/Planet.h"

#include <Eigen/Core>

namespace orbital_core::hazard {

	enum class EncounterPlaneKind {
		OpikTargetPlane,
		ScaledBPlane,
		ModifiedTargetPlane,
	};

	struct EncounterPlanePoint {
		Eigen::Vector2d coordsKm{Eigen::Vector2d::Zero()};
		BPlaneBasis basis{};
		double gamma{1.0};
		double vInfKmS{0.0};
		EncounterPlaneKind kind{EncounterPlaneKind::OpikTargetPlane};
		bool valid{false};
	};

	struct EncounterPlaneJacobian {
		Eigen::Vector2d coordsKm{Eigen::Vector2d::Zero()};
		Eigen::Matrix<double, 2, 6> dCoordsDTargetState{Eigen::Matrix<double, 2, 6>::Zero()};
		BPlaneBasis basis{};
		double gamma{1.0};
		EncounterPlaneKind kind{EncounterPlaneKind::OpikTargetPlane};
		bool valid{false};
	};

	[[nodiscard]] EncounterPlanePoint encounterPlaneAtHeliocentricState(const Vec3D& rTargetAu,
	                                                                    const Vec3D& vTargetAuPerDay,
	                                                                    const Vec3D& rPlanetAu,
	                                                                    const Vec3D& vPlanetAuPerDay, Planet planet,
	                                                                    EncounterPlaneKind kind) noexcept;

	[[nodiscard]] EncounterPlaneJacobian encounterPlaneJacobianAtHeliocentricState(
	    const Vec3D& rTargetAu, const Vec3D& vTargetAuPerDay, const Vec3D& rPlanetAu, const Vec3D& vPlanetAuPerDay,
	    Planet planet, EncounterPlaneKind kind) noexcept;

}  // namespace orbital_core::hazard
