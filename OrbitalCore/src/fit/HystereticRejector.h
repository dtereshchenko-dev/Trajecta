#pragma once

// Hysteretic outlier rejection with re-include and safeguards.
//
// Post-fit chi squares are evaluated by the fit's own observation model.
// Included observations are excluded for chi2 > rejectChiSquared;
// excluded observations are restored for chi2 <= recoverChiSquared.
// A candidate mask is rejected if it cuts the
// time span or kept-count below the configured floors.
//
// Radar measurements are never rejected automatically: their chi2 does not set
// the batch threshold, an excluded radar row may still be recovered, and a
// manual Exclude lock still removes one. Occultations follow the optical rules.

#include "orbital_core/ObsNode.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <cstdint>
#include <span>
#include <vector>

namespace Calculations {

	struct HystereticRejectConfig {
		// Recovery/rejection ratio (9.21/10), scaled without
		// increasing the user-selected sigma-squared rejection threshold.
		static constexpr double recoveryRatio = 9.21 / 10.0;
		double rejectChiSquared{3.0 * 3.0};
		double recoverChiSquared{rejectChiSquared * recoveryRatio};
		// Reject values at least 0.25 of the largest chi2 first; for small
		// samples or repeated modification counts, only the maximum (and ties).
		bool rejectLargestOnly{false};
		double minKeptFraction{2.0 / 3.0};
		double minSpanFraction{0.5};
		double retryScale{2.25};
		double maxChiSquared{3600.0};
	};

	enum class RejectStatus {
		NoChange,
		Changed,
		SafeguardRollback,
	};

	struct HystereticRejectResult {
		RejectStatus status{RejectStatus::NoChange};
		int nFlippedIn{0};
		int nFlippedOut{0};
		double appliedRejectChiSquared{0.0};
		double appliedRecoverChiSquared{0.0};
		int nKept{0};
		int nKeptBefore{0};
	};

	[[nodiscard]] ORBITAL_CORE_API HystereticRejectResult applyHystereticReject(
	    std::span<const ObsNode> allObs, std::vector<std::uint8_t>& included, std::span<const double> chiSquared,
	    const HystereticRejectConfig& cfg = {}, std::span<const std::uint8_t> locked = {});

}  // namespace Calculations
