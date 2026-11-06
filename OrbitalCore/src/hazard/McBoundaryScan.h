#pragma once

#include "NBodyPropagator.h"
#include "orbital_core/hazard/IntervalMonteCarlo.h"

namespace orbital_core::hazard::mc {

	// Presentation-only separation, not a propagation/search parameter.
	inline constexpr double eventMaxGapDays = 1.0;

	// Each worker owns its propagator. The surface-stop implementation alone
	// refines physical boundary crossings; no close-approach catalogue is built.
	void scan(McClone& clone, const IntervalMcConfig& config, NBodyPropagator& propagator, IEphemerisSource& ephemeris,
	          std::stop_token cancel);

	void summarize(IntervalMcResult& result);

}  // namespace orbital_core::hazard::mc
