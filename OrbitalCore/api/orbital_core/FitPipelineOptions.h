#pragma once

namespace orbital_core {

	enum class ObservationSelectionPolicy : int {
		// Explicit opt-in to seeding enabled clipping with ATT-EXT's working
		// mask. Excluded rows remain eligible for recovery. Without clipping,
		// include all rows. Manual observation choices always take precedence.
		ArcExtension = 0,
		// Start final selection from all eligible rows, regardless of force
		// model. ATT-EXT supplies an orbit, not the final rejection history.
		AllOriginal = 1,
	};

	enum class PerturberSelectionPolicy : int {
		AutomaticFromArc = 0,
		Configured = 1,
	};

	// Arc-extension policies, independent of numerical fit settings.
	struct FitPipelineOptions {
		ObservationSelectionPolicy observationSelectionPolicy{ObservationSelectionPolicy::AllOriginal};
		// Keep the configured force model by default; automatic arc-based
		// reduction can materially change the nominal orbit and covariance.
		PerturberSelectionPolicy perturberSelectionPolicy{PerturberSelectionPolicy::Configured};
	};

}  // namespace orbital_core
