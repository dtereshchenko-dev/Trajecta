#pragma once

#include "NBodyPropagator.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <optional>

namespace orbital_core::uncertainty {
	enum class CovariancePolicy { PositiveDefinite, PositiveSemidefinite };

	inline std::optional<Eigen::MatrixXd> covarianceFactor(const Eigen::MatrixXd& covariance, int dimension,
	                                                       CovariancePolicy policy = CovariancePolicy::PositiveDefinite)
	{
		if (covariance.rows() != dimension || covariance.cols() != dimension || !covariance.allFinite() ||
		    (covariance.diagonal().array() <= 0.0).any())
			return std::nullopt;
		const Eigen::VectorXd scale = covariance.diagonal().cwiseSqrt();
		const Eigen::MatrixXd correlation =
		    scale.cwiseInverse().asDiagonal() * covariance * scale.cwiseInverse().asDiagonal();
		if (!correlation.allFinite() || !correlation.isApprox(correlation.transpose(), 1e-10))
			return std::nullopt;
		if (policy == CovariancePolicy::PositiveSemidefinite) {
			const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(correlation);
			if (eigen.info() != Eigen::Success || eigen.eigenvalues().minCoeff() < -1e-10)
				return std::nullopt;
			return scale.asDiagonal() * eigen.eigenvectors() *
			       eigen.eigenvalues().cwiseMax(0.0).cwiseSqrt().asDiagonal();
		}
		const Eigen::LLT<Eigen::MatrixXd> llt(correlation);
		if (llt.info() != Eigen::Success)
			return std::nullopt;
		return scale.asDiagonal() * Eigen::MatrixXd(llt.matrixL());
	}

	inline Eigen::MatrixXd stateJacobian(const NBodyPropagator::PropagationResultFull& state, const NonGravModel& ng)
	{
		Eigen::MatrixXd jacobian(6, 6 + ng.fitParameterCount());
		jacobian.leftCols<6>() = state.Phi;
		for (int i = 0; i < ng.fitParameterCount(); ++i)
			jacobian.col(6 + i) = state.SensA.col(ng.sensitivityColumnForFitParameter(i));
		return jacobian;
	}
}  // namespace orbital_core::uncertainty
