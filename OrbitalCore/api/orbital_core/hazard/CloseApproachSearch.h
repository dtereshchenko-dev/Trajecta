#pragma once

#include "orbital_core/FitPipeline.h"
#include "orbital_core/FitRunControl.h"
#include "orbital_core/hazard/CloseApproach.h"
#include <optional>
#include <string>

namespace orbital_core::hazard {

	struct CloseApproachSearchRequest {
		double startMjdTdb{};
		double endMjdTdb{};
		// Empty selects per-body reporting limits, not a force-model cutoff.
		std::optional<double> maximumDistanceAu;
	};

	struct NominalCloseApproach {
		CloseApproachEvent encounter;
		double mjdTdb{};
		double relativeSpeedKmS{};
		// Full fit-epoch covariance, including NG cross terms, mapped to the moving
		// minimum and incoming Opik B-plane (equatorial J2000 pole convention).
		// Distance bounds are extrema of the filled 3-sigma ellipse, converted to
		// pericentre distances with nominal v-infinity, not d +/- 3 sigma(d).
		// Linearized estimates, not nonlinear confidence limits or impact odds.
		// Missing for invalid covariance, degenerate minima or interior passages;
		// non-hyperbolic/undefined B-planes retain only local time uncertainty.
		std::optional<double> minimumDistanceAu;
		std::optional<double> maximumDistanceAu;
		std::optional<double> timeUncertaintyDays;
		std::string uncertaintyNote;
	};

	struct CloseApproachSearchResult {
		std::vector<NominalCloseApproach> events;
		bool cancelled{};
		// The osculating pericentre of an impacting nominal may be reported;
		// no later encounters are retained after its first surface entry.
		bool stoppedAtSurface{};
	};

	[[nodiscard]] ORBITAL_CORE_API double defaultCloseApproachLimitAu(Planet planet) noexcept;

	// Searches all eight planets and the Moon using the accepted force/NG model.
	// Input and returned epochs are TDB. Throws on invalid
	// input or failed nominal propagation; cancellation never publishes partial rows.
	[[nodiscard]] ORBITAL_CORE_API CloseApproachSearchResult
	searchCloseApproaches(const FitPipelineResult& solution, const CloseApproachSearchRequest& request,
	                      const FitRunControl& control = {});

}  // namespace orbital_core::hazard
