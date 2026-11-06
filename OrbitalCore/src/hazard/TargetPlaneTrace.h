#pragma once

#include <Eigen/Core>

namespace orbital_core::hazard {

	struct TargetPlaneDifferential {
		Eigen::Vector2d coordsKm{Eigen::Vector2d::Zero()};
		Eigen::MatrixXd dCoordsDParams;
		Eigen::MatrixXd parameterCovariance;
		Eigen::VectorXd parameterStepPerSigma;
		double bodyRadiusKm{};
		double impactCrossSectionBodyRadii{1.0};
		double sigma{};
		double jd{};
	};

	struct TargetPlaneTracePoint {
		Eigen::Vector2d coordsKm{Eigen::Vector2d::Zero()};
		Eigen::Vector2d dCoordsDSigmaKm{Eigen::Vector2d::Zero()};
		double sigma{};
		double jd{};
		double bodyRadiusKm{};
		double impactCrossSectionBodyRadii{1.0};
		double distanceKm{};
		double distanceBodyRadii{};
		double stretchLovKmPerSigma{};
		double stretchLovBodyRadiiPerSigma{};
		double stretchKm{};
		double widthKm{};
		double stretchBodyRadii{};
		double widthBodyRadii{};
		double alphaRad{};
		double alphaLovRad{};
		double dDistance2DSigmaKm2{};
		double dDistance2DSigmaBodyRadii2{};
		bool impact{};
		bool valid{};
	};

	[[nodiscard]] TargetPlaneTracePoint makeTargetPlaneTracePoint(const TargetPlaneDifferential& differential);

}  // namespace orbital_core::hazard
