#pragma once

#include "orbital_core/OrbitSolution.h"
#include "orbital_core/IntegratorKind.h"
#include "orbital_core/hazard/ObsMonteCarlo.h"
#include <stop_token>

namespace orbital_core::hazard {

	struct SolutionPropagationConfig {
		NBodySettings forces;
		IntegratorKind integrator{IntegratorKind::RA15};
		CloseApproachConfig encounter;
		ObsMcConfig sampling;
		unsigned workers{15};
	};

	struct SolutionPropagationResult {
		std::vector<CloseApproachEvent> events;
		ObsMcResult mc;
	};

	// Orchestrates the existing encounter scan and covariance MC, without refitting.
	// Application owns initialized default ephemeris sources for the entire call.
	// Cancellation is checked around the scan and between MC clones.
	[[nodiscard]] ORBITAL_CORE_API SolutionPropagationResult
	propagateSolution(const OrbitSolution& solution, const SolutionPropagationConfig& config,
	                  std::stop_token cancel = {}, std::function<bool(int, int)> progress = {});

}  // namespace orbital_core::hazard
