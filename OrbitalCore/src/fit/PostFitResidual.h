#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <Eigen/Core>

namespace orbital_core::fit {

	struct PostFitResidualStatistic {
		double chiSquared{0.0};
		bool valid{false};
	};

	// Both the residual and projected formal parameter covariance are whitened
	// with the same measurement metric as the fit. No a-posteriori scaling.
	[[nodiscard]] ORBITAL_CORE_API PostFitResidualStatistic postFitResidualStatistic(
	    const Eigen::Vector2d& residual, const Eigen::Matrix2d& projectedCovariance, bool included) noexcept;

}  // namespace orbital_core::fit
