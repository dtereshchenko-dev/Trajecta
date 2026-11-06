#pragma once

#include "orbital_core/ObsNode.h"

#include <cstdint>
#include <iterator>
#include <set>
#include <span>
#include <utility>

namespace orbital_core::fit {

	// An included observation is the last in its
	// opposition when no other included observation lies within 180 days.
	// This is a local time-window test, not a transitive clustering of the arc.
	// Update after every mask change; rejected neighbours cannot protect a row.
	// Manual exclusions are supplied as absent and are never recovered here.
	class OppositionCoverage {
	   public:
		OppositionCoverage(std::span<const ObsNode> observations, std::span<const std::uint8_t> included)
		    : observations_(observations)
		{
			for (std::size_t i = 0; i < observations.size(); ++i)
				if (included[i])
					include(i);
		}

		[[nodiscard]] bool isLast(std::size_t index) const
		{
			const auto row = active_.find({observations_[index].JD, index});
			if (row == active_.end())
				return false;
			if (row != active_.begin() && row->first - std::prev(row)->first <= 180.0)
				return false;
			const auto next = std::next(row);
			return next == active_.end() || next->first - row->first > 180.0;
		}

		void include(std::size_t index) { active_.emplace(observations_[index].JD, index); }

		void exclude(std::size_t index) { active_.erase({observations_[index].JD, index}); }

	   private:
		std::span<const ObsNode> observations_;
		std::set<std::pair<double, std::size_t>> active_;
	};

}  // namespace orbital_core::fit
