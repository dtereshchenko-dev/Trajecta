#pragma once
#include "orbital_core/FitPipeline.h"

namespace orbital_core {
	struct SolutionRefinementRequest {
		const FitPipelineResult& previous;
		FitConfig config;
		ObservationPreparationOptions preparation;
		std::span<const ObservationOverride> observationOverrides{};
	};

	[[nodiscard]] ORBITAL_CORE_API FitPipelineResult refineSolution(const SolutionRefinementRequest& request,
	                                                                const FitRunControl& control = {});
}  // namespace orbital_core
