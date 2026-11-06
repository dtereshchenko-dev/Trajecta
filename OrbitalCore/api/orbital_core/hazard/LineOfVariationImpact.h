#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/hazard/ImpactRisk.h"
#include "orbital_core/hazard/LovCoordinateSystem.h"

#include <functional>
#include <vector>

class IEphemerisSource;
class NBodyPropagator;
class NonGravModel;
class Orbit;

namespace orbital_core::hazard {

	struct ORBITAL_CORE_API LineOfVariationImpactConfig {
		CloseApproachConfig closeApproach;
		ImpactRiskConfig impactRefinement;
		double sigmaMax{5.0};
		double sigmaStepMax{0.01};  // Fixed sigma step, or maximum for probability sampling.
		double probabilityStep{1.0e-7};
		bool probabilitySampling{false};  // False selects uniform sigma nodes.
		int maxPositiveSamples{10000};
		int maxRows{-1};
		int encounterStepPoints{
		    Constants::ENCOUNTER_STEP_POINTS};  // reference npoint; zero disables encounter step control.
		LovCoordinateSystem coordinateSystem{LovCoordinateSystem::Equinoctial};
		bool coordinateScaling{true};  // Scale coordinates independently of posterior covariance scaling.
	};

	[[nodiscard]] ORBITAL_CORE_API std::vector<ImpactRiskRow> traceLineOfVariationImpacts(
	    const Orbit& nominalOrbit, const StateCov6x6& cov0, const NBodyPropagator& propagator,
	    const LineOfVariationImpactConfig& cfg = {}, IEphemerisSource* eph = nullptr,
	    std::function<bool(int, int)> progressCb = {});

	[[nodiscard]] ORBITAL_CORE_API std::vector<ImpactRiskRow> traceLineOfVariationImpacts(
	    const Orbit& nominalOrbit, const FitParameterCovariance& parameterCovariance,
	    const NonGravModel& nominalNonGrav, const NBodyPropagator& propagator,
	    const LineOfVariationImpactConfig& cfg = {}, IEphemerisSource* eph = nullptr,
	    std::function<bool(int, int)> progressCb = {});

}  // namespace orbital_core::hazard
