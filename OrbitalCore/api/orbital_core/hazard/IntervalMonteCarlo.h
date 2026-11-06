#pragma once

#include "orbital_core/FitPipeline.h"
#include "orbital_core/hazard/Planet.h"
#include <array>
#include <filesystem>
#include <functional>
#include <stop_token>

namespace orbital_core::hazard {

	enum class McMethod { Covariance, Observations };
	enum class McCloneStatus { Unstarted, Complete, FitFailed, PropagationFailed, Cancelled };

	struct IntervalMcConfig {
		McMethod method{McMethod::Covariance};
		Planet target{Planet::Earth};
		double startMjdTdb{};
		double endMjdTdb{};
		double altitudeKm{100.0};
		int samples{10000};
		std::uint64_t seed{42};
		unsigned workers{15};
	};

	struct McBoundaryEntry {
		double mjdTdb{};
		Vec3D relativePositionKm;
		Vec3D relativeVelocityKmS;
		double latitudeDeg{};
		double longitudeDeg{};
	};

	struct McClone {
		int index{-1};  // -1 is the nominal, never a random trial.
		McCloneStatus status{McCloneStatus::Unstarted};
		double epochMjdTdb{};
		std::array<double, 6> state{};  // heliocentric ecliptic AU, AU/day
		NonGravModel nonGrav;
		bool entryBeforeInterval{};
		std::optional<McBoundaryEntry> entry;
		std::string error;
	};

	// Temporal groups of physical boundary entries: split only at a gap greater
	// than 24 hours between consecutive crossings, never at UTC midnight.
	// This presentation convention changes no trajectories, selection or search.
	// Nominal entries do not create trial events or contribute to probabilities.
	struct McEvent {
		double startMjdTdb{};
		double endMjdTdb{};
		int entryClones{};
		double probability{};
		double probabilityLow{};
		double probabilityHigh{};
	};

	struct IntervalMcResult {
		IntervalMcConfig config;
		// Captured during computation, never inferred from the runtime at export.
		std::string bodyFixedFrame;
		std::string bodyFixedModel;
		McClone nominal;
		std::vector<McClone> clones;
		std::vector<McEvent> events;
		bool cancelled{};
		int completed{};
		int failedFit{};
		int failedPropagation{};
		int cancelledClones{};
		int unstarted{};
		int entryClones{};
		int entriesBeforeInterval{};
	};

	// Borrows the immutable accepted solution for the synchronous call. Samples at
	// fitEpochOrbit with full objective-scaled COV or measurement-noise refits with
	// the accepted fixed mask. No IOD, clipping, prior or covariance rescaling.
	// All times are MJD TDB, and both bounds are inclusive. Per-clone failures remain
	// explicit; only Complete trials enter empirical probabilities (Wilson 95%).
	// Cancellation invalidates partial trajectories, but retains completed trials.
	// The target is one body; this is not an all-body first-collision competition.
	[[nodiscard]] ORBITAL_CORE_API IntervalMcResult runIntervalMonteCarlo(const FitPipelineResult& solution,
	                                                                      const IntervalMcConfig& config,
	                                                                      std::stop_token cancel = {},
	                                                                      std::function<void(int, int)> progress = {});

	// Writes a new immutable numerical archive and an interactive HTML report.
	// No legacy corridor format is impersonated. Caller supplies a new directory.
	ORBITAL_CORE_API void saveIntervalMonteCarlo(const std::filesystem::path& directory,
	                                             const FitPipelineResult& source, const IntervalMcResult& result);

	// Rebuilds only report.html from an existing result.json, without propagation.
	// Includes complete-trial boundary coordinates, UTC times and a body-fixed map.
	// Numerical files are never modified by report rendering.
	ORBITAL_CORE_API void renderIntervalMonteCarloReport(const std::filesystem::path& directory);

}  // namespace orbital_core::hazard
