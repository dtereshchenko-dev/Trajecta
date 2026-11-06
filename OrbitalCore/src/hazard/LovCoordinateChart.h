#pragma once

#include "orbital_core/FitResult.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/hazard/LovCoordinateSystem.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <Eigen/Core>
#include <optional>

namespace orbital_core::hazard {

	struct LovCoordinateChart {
		Eigen::VectorXd nominalCoordinates;
		FitParameterCovariance coordinateCovariance;
		NonGravModel nonGrav;
		double epochTdb{};
		LovCoordinateSystem coordinateSystem{LovCoordinateSystem::Equinoctial};
		bool valid{false};
	};

	[[nodiscard]] ORBITAL_CORE_API Eigen::VectorXd lovCoordinateUnits(LovCoordinateSystem coordinateSystem,
	                                                                  const Eigen::VectorXd& coordinates, bool scaling);

	[[nodiscard]] ORBITAL_CORE_API LovCoordinateChart makeEquinoctialLovChart(
	    const Orbit& orbit, const FitParameterCovariance& physicalCovariance, const NonGravModel& nonGrav);

	[[nodiscard]] ORBITAL_CORE_API LovCoordinateChart makeKeplerianLovChart(
	    const Orbit& orbit, const FitParameterCovariance& physicalCovariance, const NonGravModel& nonGrav);

	[[nodiscard]] ORBITAL_CORE_API LovCoordinateChart makeCometaryLovChart(
	    const Orbit& orbit, const FitParameterCovariance& physicalCovariance, const NonGravModel& nonGrav);

	// Select COT for near-parabolic and open conics, otherwise preserve the
	// requested chart. The chosen chart remains fixed throughout a LOV trace.
	[[nodiscard]] ORBITAL_CORE_API LovCoordinateChart makeAutomaticLovChart(
	    const Orbit& orbit, const FitParameterCovariance& physicalCovariance, const NonGravModel& nonGrav,
	    LovCoordinateSystem preferred = LovCoordinateSystem::Equinoctial);

	[[nodiscard]] ORBITAL_CORE_API std::optional<Eigen::VectorXd> cometaryCoordinates(const Orbit& orbit,
	                                                                                  const NonGravModel& nonGrav);

	[[nodiscard]] ORBITAL_CORE_API std::optional<Eigen::VectorXd> keplerianCoordinates(const Orbit& orbit,
	                                                                                   const NonGravModel& nonGrav);

	[[nodiscard]] ORBITAL_CORE_API std::optional<Eigen::VectorXd> equinoctialCoordinates(const Orbit& orbit,
	                                                                                     const NonGravModel& nonGrav);

	[[nodiscard]] ORBITAL_CORE_API std::optional<Eigen::VectorXd> physicalParametersFromChart(
	    const LovCoordinateChart& chart, const Eigen::VectorXd& coordinates);

	[[nodiscard]] ORBITAL_CORE_API bool physicalParametersFromChart(const LovCoordinateChart& chart,
	                                                                const Eigen::VectorXd& coordinates,
	                                                                Eigen::Ref<Eigen::VectorXd> output);

	[[nodiscard]] ORBITAL_CORE_API std::optional<Eigen::MatrixXd> physicalJacobianFromChart(
	    const LovCoordinateChart& chart, const Eigen::VectorXd& coordinates);

	[[nodiscard]] ORBITAL_CORE_API bool physicalJacobianFromChart(const LovCoordinateChart& chart,
	                                                              const Eigen::VectorXd& coordinates,
	                                                              Eigen::Ref<Eigen::MatrixXd> output);

	[[nodiscard]] ORBITAL_CORE_API std::optional<Eigen::VectorXd> physicalStepPerSigmaFromChart(
	    const LovCoordinateChart& chart, const Eigen::VectorXd& coordinates,
	    const Eigen::VectorXd& coordinateStepPerSigma);

}  // namespace orbital_core::hazard
