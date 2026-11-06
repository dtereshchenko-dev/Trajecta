#include "obs/OverObservationReweight.h"
#include <string>
#include <vector>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <string_view>

namespace orbital_core::obs {
	namespace {

		constexpr double kSameStationMaxGapDays = 8.0 / 24.0;
		constexpr std::size_t kVfcc17BatchCeiling = 4;

		[[nodiscard]] bool isEligible(const PreparedObservation& observation)
		{
			constexpr std::string_view techniques = "CcPANMTEeVSOHnB";
			return !observation.radar && !observation.observer.obsCode.empty() &&
			       techniques.find(observation.note2) != std::string_view::npos;
		}


	}  // namespace

	std::vector<SameStationReweightFactor> computeSameStationOverObservationFactors(
	    const std::vector<PreparedObservation>& observations, AstrometricWeightProfile profile)
	{
		const std::size_t count = observations.size();
		if (count == 0 || profile != AstrometricWeightProfile::Vfcc17)
			return {};

		std::vector<SameStationReweightFactor> factors(count);
		std::map<std::string, std::vector<std::size_t>> byStation;
		for (std::size_t i = 0; i < count; ++i) {
			if (isEligible(observations[i]))
				byStation[observations[i].observer.obsCode].push_back(i);
		}

		bool hasDownweightedObservation = false;

		for (auto& [station, indices] : byStation) {
			(void)station;
			std::stable_sort(indices.begin(), indices.end(), [&observations](std::size_t lhs, std::size_t rhs) {
				return observations[lhs].jdTdb < observations[rhs].jdTdb;
			});

			std::size_t batchBegin = 0;
			while (batchBegin < indices.size()) {
				if (observations[indices[batchBegin]].isSatellite) {
					++batchBegin;
					continue;
				}
				std::size_t batchEnd = batchBegin + 1;
				while (batchEnd < indices.size()) {
					const double gap =
					    observations[indices[batchEnd]].jdTdb - observations[indices[batchEnd - 1]].jdTdb;
					if (gap > kSameStationMaxGapDays)
						break;
					++batchEnd;
				}

				const std::size_t batchSize = batchEnd - batchBegin;
				if (batchSize >= kVfcc17BatchCeiling + 1) {
					const double factor =
					    std::sqrt(static_cast<double>(kVfcc17BatchCeiling) / static_cast<double>(batchSize));
					for (std::size_t k = batchBegin; k < batchEnd; ++k) {
						const std::size_t index = indices[k];
						const bool observerSigma = observations[index].observerSigma;
						if (!observations[index].forceWeightRa && !observerSigma)
							factors[index].ra = factor;
						if (!observations[index].forceWeightDec && !observerSigma)
							factors[index].dec = factor;
						hasDownweightedObservation =
						    hasDownweightedObservation || factors[index].ra < 1.0 || factors[index].dec < 1.0;
					}
				}

				batchBegin = batchEnd;
			}
		}

		return hasDownweightedObservation ? factors : std::vector<SameStationReweightFactor>{};
	}

	SameStationReweightStats summarizeSameStationOverObservationFactors(
	    const std::vector<SameStationReweightFactor>& factors)
	{
		SameStationReweightStats stats;
		stats.observationCount = factors.size();
		if (factors.empty())
			return stats;

		double sum = 0.0;
		stats.minFactor = std::numeric_limits<double>::infinity();
		stats.maxFactor = 0.0;

		for (const auto& factor : factors) {
			const double minCoordinate = std::min(factor.ra, factor.dec);
			const double maxCoordinate = std::max(factor.ra, factor.dec);
			stats.minFactor = std::min(stats.minFactor, minCoordinate);
			stats.maxFactor = std::max(stats.maxFactor, maxCoordinate);
			sum += 0.5 * (factor.ra + factor.dec);
			if (minCoordinate < 0.999999) {
				++stats.downweightedCount;
				const double effectiveNearby =
				    static_cast<double>(kVfcc17BatchCeiling) / (minCoordinate * minCoordinate);
				stats.maxEffectiveNearby = std::max(stats.maxEffectiveNearby, effectiveNearby);
			}
		}

		stats.meanFactor = sum / static_cast<double>(factors.size());
		return stats;
	}

}  // namespace orbital_core::obs
