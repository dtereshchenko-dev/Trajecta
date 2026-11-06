#pragma once

#include "obs/PreparedObservation.h"
#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <cstddef>
#include <vector>

namespace orbital_core::obs {

	struct SameStationReweightFactor {
		double ra{1.0};
		double dec{1.0};
	};

	struct SameStationReweightStats {
		std::size_t observationCount{0};
		std::size_t downweightedCount{0};
		double maxEffectiveNearby{1.0};
		double minFactor{1.0};
		double meanFactor{1.0};
		double maxFactor{1.0};
	};

	[[nodiscard]] ORBITAL_CORE_API std::vector<SameStationReweightFactor> computeSameStationOverObservationFactors(
	    const std::vector<PreparedObservation>& observations, AstrometricWeightProfile profile);

	[[nodiscard]] ORBITAL_CORE_API SameStationReweightStats
	summarizeSameStationOverObservationFactors(const std::vector<SameStationReweightFactor>& factors);

}  // namespace orbital_core::obs
