#pragma once

#include "orbital_core/hazard/BoundaryCrossing.h"
#include <optional>
#include <string>
#include <vector>

namespace orbital_core::hazard {
	struct SemilinearContour {
		int sigma{};
		// Missing points split the curve; they must not be joined across a miss.
		std::vector<std::optional<BoundaryCrossing>> points;
		int missed{}, failed{};
	};

	struct SemilinearContours {
		std::string status{"not_computed"};
		std::string reason;
		std::vector<SemilinearContour> contours;
	};
}  // namespace orbital_core::hazard
