#pragma once

#include "orbital_core/FitSession.h"
#include "obs/PreparedObservation.h"

namespace orbital_core::fit {

	// Internal entry for a pipeline-owned, fully prepared observation pool.
	// Row order and weighting configuration must match request/config exactly.
	// The public FitSession entry prepares its own pool and uses the same runner.
	[[nodiscard]] FitResult runPreparedFitSession(const FitConfig& config, const FitRequest& request,
	                                              std::span<const obs::PreparedObservation> prepared,
	                                              const FitRunControl& control);

}  // namespace orbital_core::fit
