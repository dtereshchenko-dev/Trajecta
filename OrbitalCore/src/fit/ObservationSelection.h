#pragma once
#include "orbital_core/FitSession.h"
#include <vector>

namespace orbital_core::fit {
	struct ObservationSelection {
		std::vector<std::uint8_t> included;
		std::vector<std::uint8_t> locked;
	};

	// Replaced discovery records (MPC note 2 'X'/'x') start excluded and locked.
	// They never enter IOD, ATT-EXT, a fit or automatic recovery unless a manual
	// Include override selects them.
	ORBITAL_CORE_API ObservationSelection resolveObservationSelection(std::span<const ObsNode> observations,
	                                                                  std::span<const std::uint8_t> initial,
	                                                                  std::span<const ObservationOverride> overrides);
}  // namespace orbital_core::fit
