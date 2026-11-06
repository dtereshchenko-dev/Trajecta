#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/CloseApproach.h"
#include "orbital_core/hazard/Planet.h"

#include <Eigen/Core>
#include <optional>

namespace orbital_core::hazard {

	[[nodiscard]] double circularEarthVInfinityKmS(const Orbit& orbit) noexcept;

	[[nodiscard]] double closeApproachVInfinityKmS(const Vec3D& vRelAuPerDay, double rAu, Planet planet) noexcept;

	void refreshCloseApproachKinematics(CloseApproachEvent& event) noexcept;

	struct InteriorPericenterMapping {
		Vec3D relativePosition;
		Vec3D relativeVelocity;
		Eigen::Matrix<double, 6, 6> stateTransition{Eigen::Matrix<double, 6, 6>::Identity()};
		double deltaDays{};
		bool valid{false};

		// Derivative of each input's own pericenter state, not a fixed-time STM.
		// Singular for a degenerate radial minimum; never use for general transport.
		[[nodiscard]] ORBITAL_CORE_API std::optional<Eigen::Matrix<double, 6, 6>> eventJacobian(
		    Planet planet) const noexcept;
	};

	// Reference close_app/tp_fser regularization for an encounter that crosses
	// the body's reference sphere. The osculating planetocentric two-body
	// orbit is continued from the entry state to pericenter together with
	// its fixed-time state transition matrix.
	[[nodiscard]] ORBITAL_CORE_API InteriorPericenterMapping mapInteriorEncounterToPericenter(
	    const Vec3D& relativePosition, const Vec3D& relativeVelocity, Planet planet) noexcept;

}  // namespace orbital_core::hazard
