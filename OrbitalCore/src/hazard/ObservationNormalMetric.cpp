#include "hazard/ObservationNormalMetric.h"
#include "fit/APosterioriCovarianceScaling.h"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbital_core::hazard {
	ObservationNormalMetric makeObservationNormalMetric(const Eigen::MatrixXd& jacobian, double chiSquare)
	{
		ObservationNormalMetric out;
		const Eigen::Index numResiduals = jacobian.rows();
		const Eigen::Index numParameters = jacobian.cols();
		if (numResiduals <= numParameters || numParameters <= 0 || !jacobian.array().isFinite().all() ||
		    !std::isfinite(chiSquare) || chiSquare < 0.0) {
			return out;
		}

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, Eigen::ComputeThinV);
		if (svd.info() != Eigen::Success)
			return out;

		const Eigen::VectorXd singular = svd.singularValues();
		if (singular.size() != numParameters || !(singular(0) > 0.0)) {
			return out;
		}

		const double singularFloor = singular(0) * std::numeric_limits<double>::epsilon() *
		                             static_cast<double>(std::max(numResiduals, numParameters));
		out.eigenvalueFloor = singularFloor * singularFloor;
		if (!(singular(numParameters - 1) > singularFloor))
			return out;

		out.eigenvalues = singular.reverse().array().square().matrix();
		out.eigenvectors = svd.matrixV().rowwise().reverse();
		out.covarianceScale =
		    fit::aPosterioriVarianceFactor(chiSquare, static_cast<int>(numResiduals), static_cast<int>(numParameters));

		out.covariance = Eigen::MatrixXd::Zero(numParameters, numParameters);
		for (Eigen::Index i = 0; i < out.eigenvalues.size(); ++i) {
			const double value = out.eigenvalues(i);
			if (!std::isfinite(value) || value <= out.eigenvalueFloor)
				continue;
			const Eigen::VectorXd axis = out.eigenvectors.col(i);
			out.covariance.noalias() += (out.covarianceScale / value) * (axis * axis.transpose());
		}

		out.valid = out.covariance.array().isFinite().all();
		return out;
	}

}  // namespace orbital_core::hazard
