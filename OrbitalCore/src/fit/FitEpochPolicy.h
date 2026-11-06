#pragma once

#include "orbital_core/FitConfig.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <span>

namespace orbital_core::fit {

	// Existing reporting convention: choose a JD ending in .5 from the arc's
	// TDB endpoints, biased toward its last observation. Empty input returns 0.
	[[nodiscard]] double defaultOutputEpochJdTdb(std::span<const ObsNode> observations);

	// Midpoint of the observed arc in TDB. The Mjd form is the exact one; the
	// JD form is its lossy view, kept for callers still on that convention.
	[[nodiscard]] double arcMidpointEpochTdbMjd(std::span<const ObsNode> observations);

	[[nodiscard]] double arcMidpointEpochTdb(std::span<const ObsNode> observations);

	[[nodiscard]] Orbit_p rebaseInitialOrbitToFitEpoch(const Orbit_p& initial, std::span<const ObsNode> observations,
	                                                   const FitConfig& cfg, const NonGravModel* warmStartNonGrav,
	                                                   bool disableNonGrav);

}  // namespace orbital_core::fit
