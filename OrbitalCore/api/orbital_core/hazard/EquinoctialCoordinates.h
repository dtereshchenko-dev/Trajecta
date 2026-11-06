#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <Eigen/Core>
#include <optional>

namespace orbital_core::hazard::equinoctial {
	using Vector6 = Eigen::Matrix<double, 6, 1>;
	using Matrix6 = Eigen::Matrix<double, 6, 6>;

	// Elliptic prograde chart [a,h,k,p,q,lambda], AU/radians; CAR is AU, AU/day.
	[[nodiscard]] ORBITAL_CORE_API std::optional<Vector6> fromCartesian(const Vector6& state);
	[[nodiscard]] ORBITAL_CORE_API std::optional<Vector6> toCartesian(const Vector6& elements,
	                                                                  Matrix6* jacobian = nullptr);
	[[nodiscard]] ORBITAL_CORE_API std::optional<Matrix6> cartesianToElementsJacobian(const Vector6& elements);
}  // namespace orbital_core::hazard::equinoctial
