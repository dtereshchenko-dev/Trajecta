#include "fit/ObservationSelection.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace orbital_core::fit {
	ObservationSelection resolveObservationSelection(std::span<const ObsNode> observations,
	                                                 std::span<const std::uint8_t> initial,
	                                                 std::span<const ObservationOverride> overrides)
	{
		if (!initial.empty() && initial.size() != observations.size())
			throw std::invalid_argument("Observation selection mask size mismatch");
		ObservationSelection result;
		result.included.assign(observations.size(), 1);
		for (std::size_t i = 0; i < initial.size(); ++i)
			result.included[i] = initial[i] ? 1 : 0;
		for (std::size_t i = 0; i < observations.size(); ++i) {
			if (observations[i].note2 != 'X' && observations[i].note2 != 'x')
				continue;
			if (result.locked.empty())
				result.locked.assign(observations.size(), 0);
			result.included[i] = 0;
			result.locked[i] = 1;
		}
		if (!overrides.empty()) {
			std::unordered_map<ObservationId, std::size_t> indices;
			for (std::size_t i = 0; i < observations.size(); ++i) {
				if (observations[i].id.empty() || !indices.emplace(observations[i].id, i).second)
					throw std::invalid_argument("Observations require unique non-empty identities");
			}
			if (result.locked.empty())
				result.locked.assign(observations.size(), 0);
			std::unordered_set<ObservationId> seen;
			for (const auto& override : overrides) {
				const auto found = indices.find(override.id);
				if (found == indices.end() || !seen.insert(override.id).second)
					throw std::invalid_argument("Unknown or duplicate observation override");
				const auto i = found->second;
				switch (override.disposition) {
					case ObservationDisposition::Auto:
						break;
					case ObservationDisposition::Include: {
						const auto& o = observations[i];
						if (!std::isfinite(o.JD) ||
						    (o.radar ? !o.radar->valid() : (!std::isfinite(o.RA) || !std::isfinite(o.DEC))))
							throw std::invalid_argument("Cannot include invalid astrometry");
						result.included[i] = 1;
						result.locked[i] = 1;
						break;
					}
					case ObservationDisposition::Exclude:
						result.included[i] = 0;
						result.locked[i] = 1;
						break;
					default:
						throw std::invalid_argument("Invalid observation disposition");
				}
			}
		}
		if (std::ranges::count(result.included, std::uint8_t{1}) < 3)
			throw std::invalid_argument("Selection keeps fewer than 3 observations");
		return result;
	}
}  // namespace orbital_core::fit
