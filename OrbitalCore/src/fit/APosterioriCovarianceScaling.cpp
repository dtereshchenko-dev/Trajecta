#include "fit/APosterioriCovarianceScaling.h"

#include <algorithm>
#include <cmath>

namespace orbital_core::fit {

	double aPosterioriVarianceFactor(double chiSquare, int numUsedResiduals, int numSolvedParameters) noexcept
	{
		if (!std::isfinite(chiSquare) || chiSquare < 0.0 || numUsedResiduals <= numSolvedParameters) {
			return 1.0;
		}

		const double used = static_cast<double>(numUsedResiduals);
		const double dof = static_cast<double>(numUsedResiduals - numSolvedParameters);
		return std::max(used, chiSquare) / dof;
	}

}  // namespace orbital_core::fit
