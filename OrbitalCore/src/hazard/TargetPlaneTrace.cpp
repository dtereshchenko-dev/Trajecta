#include "hazard/TargetPlaneTrace.h"

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <cmath>

namespace orbital_core::hazard {

	TargetPlaneTracePoint makeTargetPlaneTracePoint(const TargetPlaneDifferential& differential)
	{
		TargetPlaneTracePoint out;
		out.sigma = differential.sigma;
		out.jd = differential.jd;
		out.bodyRadiusKm = differential.bodyRadiusKm;
		out.impactCrossSectionBodyRadii = differential.impactCrossSectionBodyRadii;

		if (differential.bodyRadiusKm <= 0.0 || differential.impactCrossSectionBodyRadii <= 0.0 ||
		    differential.dCoordsDParams.rows() != 2 ||
		    differential.dCoordsDParams.cols() != differential.parameterStepPerSigma.size() ||
		    differential.parameterCovariance.rows() != differential.parameterStepPerSigma.size() ||
		    differential.parameterCovariance.cols() != differential.parameterStepPerSigma.size() ||
		    !differential.coordsKm.array().isFinite().all() || !differential.dCoordsDParams.array().isFinite().all() ||
		    !differential.parameterCovariance.array().isFinite().all() ||
		    !differential.parameterStepPerSigma.array().isFinite().all())
			return out;

		out.coordsKm = differential.coordsKm;
		out.dCoordsDSigmaKm = differential.dCoordsDParams * differential.parameterStepPerSigma;

		out.distanceKm = out.coordsKm.norm();
		out.distanceBodyRadii = out.distanceKm / differential.bodyRadiusKm;
		out.stretchLovKmPerSigma = out.dCoordsDSigmaKm.norm();
		out.stretchLovBodyRadiiPerSigma = out.stretchLovKmPerSigma / differential.bodyRadiusKm;

		const Eigen::MatrixXd symParameterCovariance =
		    0.5 * (differential.parameterCovariance + differential.parameterCovariance.transpose());
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> covarianceSolver(symParameterCovariance);
		if (covarianceSolver.info() != Eigen::Success)
			return out;

		const Eigen::VectorXd covarianceRootValues = covarianceSolver.eigenvalues().cwiseMax(0.0).cwiseSqrt();
		const Eigen::MatrixXd targetPlaneRoot =
		    differential.dCoordsDParams * covarianceSolver.eigenvectors() * covarianceRootValues.asDiagonal();
		Eigen::JacobiSVD<Eigen::MatrixXd> targetPlaneSvd(targetPlaneRoot, Eigen::ComputeFullU);
		if (targetPlaneSvd.info() != Eigen::Success || targetPlaneSvd.singularValues().size() != 2)
			return out;

		out.stretchKm = targetPlaneSvd.singularValues()(0);
		out.widthKm = targetPlaneSvd.singularValues()(1);
		out.widthBodyRadii = out.widthKm / differential.bodyRadiusKm;
		out.stretchBodyRadii = out.stretchKm / differential.bodyRadiusKm;

		Eigen::Vector2d majorAxis = targetPlaneSvd.matrixU().col(0);
		if (majorAxis.dot(out.dCoordsDSigmaKm) < 0.0)
			majorAxis = -majorAxis;
		out.alphaRad = std::atan2(-majorAxis.x(), majorAxis.y());

		if (out.stretchLovKmPerSigma > 0.0) {
			out.alphaLovRad = std::atan2(-out.dCoordsDSigmaKm.x(), out.dCoordsDSigmaKm.y());
		}

		out.dDistance2DSigmaKm2 = 2.0 * out.coordsKm.dot(out.dCoordsDSigmaKm);
		out.dDistance2DSigmaBodyRadii2 =
		    out.dDistance2DSigmaKm2 / (differential.bodyRadiusKm * differential.bodyRadiusKm);
		out.impact = out.distanceKm <= differential.impactCrossSectionBodyRadii * differential.bodyRadiusKm;
		out.valid = std::isfinite(out.distanceKm) && std::isfinite(out.stretchLovKmPerSigma) &&
		            std::isfinite(out.widthKm) && out.stretchKm > 0.0 && out.widthKm > 0.0 &&
		            std::isfinite(out.dDistance2DSigmaKm2);
		return out;
	}

}  // namespace orbital_core::hazard
