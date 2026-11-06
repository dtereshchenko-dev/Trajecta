#include "fit/PostFitResidual.h"

#include <Eigen/Cholesky>
#include <cmath>

namespace orbital_core::fit {

	PostFitResidualStatistic postFitResidualStatistic(const Eigen::Vector2d& residual,
	                                                  const Eigen::Matrix2d& projectedCovariance,
	                                                  bool included) noexcept
	{
		if (!residual.allFinite() || !projectedCovariance.allFinite())
			return {};
		const Eigen::Matrix2d symmetric = 0.5 * (projectedCovariance + projectedCovariance.transpose()).eval();
		const Eigen::Matrix2d covariance = Eigen::Matrix2d::Identity() + (included ? -1.0 : 1.0) * symmetric;
		const Eigen::LLT<Eigen::Matrix2d> factor(covariance);
		if (factor.info() != Eigen::Success)
			return {};
		const double chiSquared = residual.dot(factor.solve(residual));
		if (!std::isfinite(chiSquared) || chiSquared < 0.0)
			return {};
		return {chiSquared, true};
	}

}  // namespace orbital_core::fit
