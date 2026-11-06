#pragma once

#include "hazard/TargetPlaneTrace.h"
#include "orbital_core/OrbitalCoreExport.h"

namespace orbital_core::hazard {

	struct LovImpactProbabilityInput {
		TargetPlaneTracePoint tracePoint;
		double impactCrossSectionBodyRadii{1.0};
	};

	struct LovImpactProbabilityResult {
		double probability{};
		double sigmaImpact{};
		bool candidate{};
		bool valid{};
	};

	[[nodiscard]] ORBITAL_CORE_API LovImpactProbabilityResult
	lovImpactProbability(const LovImpactProbabilityInput& input);

}  // namespace orbital_core::hazard
