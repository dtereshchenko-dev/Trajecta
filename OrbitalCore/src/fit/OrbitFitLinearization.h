#pragma once

#include "obs/AstrometricResidual.h"
#include "orbital_core/FitResult.h"
#include <Eigen/Core>
#include <vector>

namespace orbital_core::fit {

	struct PostFitDiagnostics {
		Residuals residuals;
		FitParameterCovariance formalCovariance;
		bool valid{false};
	};

	// Canonical velocity parameters; raw angular residuals in radians.
	struct LinearizationSnapshot {
		std::vector<double> params;
		std::vector<double> rawResiduals;
		Eigen::MatrixXd rawJacobian;
		std::vector<obs::AstrometricResidualTransform> transforms;
		bool valid{false};

		[[nodiscard]] bool matches(const std::vector<double>& values) const { return valid && params == values; }
	};

	// Whitened observation-only system, with velocity columns in AU/day.
	struct WeightedObservationLinearization {
		std::vector<double> residuals;
		Eigen::MatrixXd jacobian;
		double sumSquares{0.0};
		bool valid{false};
	};

}  // namespace orbital_core::fit
