#pragma once

#include "orbital_core/Types.h"
#include "orbital_core/hazard/ImpactRisk.h"

#include <Eigen/Core>

#include <limits>

namespace orbital_core::hazard {

	struct VirtualImpactorRefinementInput {
		const ImpactRiskConfig* cfg{};
		Eigen::MatrixXd dCoordsDParams;
		FitParameterCovariance parameterCovariance;
		Eigen::Vector2d currentCoordsKm{Eigen::Vector2d::Zero()};
		Eigen::Vector2d targetCoordsKm{Eigen::Vector2d::Zero()};
		Eigen::VectorXd initialParameters;
	};

	struct VirtualImpactorRefinementResult {
		Eigen::VectorXd parameters;
		FitParameterCovariance parameterCovariance;
		double observationChiSquare{std::numeric_limits<double>::quiet_NaN()};
		bool attempted{};
		bool valid{};
		bool badConditioning{};
		const char* failureStage{"none"};
	};

	[[nodiscard]] Eigen::Vector2d impactBoundaryTarget(const Eigen::Vector2d& coordsKm, double impactRadiusKm) noexcept;

	[[nodiscard]] VirtualImpactorRefinementResult refineVirtualImpactorWithObservations(
	    const VirtualImpactorRefinementInput& input);

}  // namespace orbital_core::hazard
