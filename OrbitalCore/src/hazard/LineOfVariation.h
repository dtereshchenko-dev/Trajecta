#pragma once

#include <Eigen/Core>

#include <vector>

namespace orbital_core::hazard {

	struct LovWeakDirection {
		Eigen::VectorXd parameterStepPerSigma;
		double covarianceEigenvalue{0.0};
		bool valid{false};
	};

	struct LovSamplingConfig {
		double sigmaMax{};
		double sigmaStepMax{};
		double probabilityStep{};
		double diskRadiusBodyRadii{};
		int maxPositiveSamples{};
		bool probabilitySampling{false};
	};

	[[nodiscard]] LovWeakDirection weakDirectionFromCovariance(const Eigen::MatrixXd& covariance);

	[[nodiscard]] std::vector<double> sampleLovSigmas(const LovSamplingConfig& config);

}  // namespace orbital_core::hazard
