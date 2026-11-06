#pragma once

#include <Eigen/Core>

namespace orbital_core::hazard {

	struct ObservationNormalMetric {
		Eigen::VectorXd eigenvalues;
		Eigen::MatrixXd eigenvectors;
		Eigen::MatrixXd covariance;
		double eigenvalueFloor{};
		double covarianceScale{1.0};
		bool valid{};
	};

	[[nodiscard]] ObservationNormalMetric makeObservationNormalMetric(const Eigen::MatrixXd& jacobian,
	                                                                  double chiSquare);

}  // namespace orbital_core::hazard
