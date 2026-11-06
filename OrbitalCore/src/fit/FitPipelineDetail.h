#pragma once

#include "orbital_core/ObsNode.h"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>

namespace orbital_core::fit {

	// Inclusive indices; timestamp ties at either boundary remain in the arc.
	struct ObservationSubarc {
		std::size_t first;
		std::size_t last;
	};

	[[nodiscard]] inline std::optional<ObservationSubarc> mapIodSubarc(std::span<const ObsNode> iod,
	                                                                   std::span<const ObsNode> all,
	                                                                   ObservationSubarc subarc)
	{
		if (subarc.first > subarc.last || subarc.last >= iod.size())
			return std::nullopt;

		const auto first = std::ranges::lower_bound(all, iod[subarc.first].JD, {}, &ObsNode::JD);
		const auto end = std::ranges::upper_bound(all, iod[subarc.last].JD, {}, &ObsNode::JD);
		if (first == all.end() || end == all.begin() || first >= end)
			return std::nullopt;

		return ObservationSubarc{static_cast<std::size_t>(first - all.begin()),
		                         static_cast<std::size_t>(end - all.begin() - 1)};
	}

}  // namespace orbital_core::fit
