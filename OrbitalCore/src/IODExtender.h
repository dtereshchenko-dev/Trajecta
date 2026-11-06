#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "OrbitOptimizer.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Calculations {

	struct ORBITAL_CORE_API ExtensionResult {
		Orbit_p orbit;
		std::vector<ObsNode> keptObs;
		std::vector<std::uint8_t> includedMask;
		std::vector<std::size_t> keptObsIndices;
		std::size_t firstIncluded{0};
		std::size_t lastIncluded{0};
		// Input range the extension linked.  Exclusions inside it are outliers of
		// the auto-reject; first/lastIncluded can lie inside it when an edge
		// observation is one of them.
		std::size_t arcFirst{0};
		std::size_t arcLast{0};
		std::uint16_t configuredPerturberMask{0};
		std::uint16_t effectivePerturberMask{0};
		std::uint16_t automaticPerturberMask{0};
		std::uint16_t blockedPerturberMask{0};
		int iterations{0};
		double finalRMS{0.0};
		bool converged{false};
	};

	// True when the extension linked the whole input, so every exclusion is an
	// outlier decision rather than an arc restriction.
	[[nodiscard]] inline bool spansInputRange(const ExtensionResult& result, std::size_t observationCount) noexcept
	{
		return observationCount > 0 && result.arcFirst == 0 && result.arcLast + 1 >= observationCount;
	}

	// Scoring, auto-reject and subset fits consume rows selected from the pipeline's
	// full prepared pool. Do not prepare again or recount batches after filtering
	// radar/manual exclusions from the optical seed pool.
	// Each subset fit runs the backend to completion. A converged fit is reused
	// while its included rows and effective force mask remain unchanged; detected
	// perturbers require another fit only when they change that effective mask.
	// Failed fits cannot be promoted by a score plateau or an auto-reject round.
	[[nodiscard]] ORBITAL_CORE_API ExtensionResult attemptExtensions(
	    const Orbit_p& seed, const std::vector<ObsNode>& allObs, std::size_t subarcStart, std::size_t subarcEnd,
	    const OrbitOptimizerConfig& baseCfg, std::span<const orbital_core::obs::PreparedObservation> prepared);

}  // namespace Calculations
