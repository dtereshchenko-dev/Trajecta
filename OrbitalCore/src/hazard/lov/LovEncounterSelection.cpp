#include "hazard/lov/LovTraceDetail.h"
#include <optional>
#include <span>
#include <vector>

#include <algorithm>
#include <cmath>

namespace orbital_core::hazard::lov_detail {

	std::optional<std::size_t> selectReturnTargetPlaneMinimum(std::span<const CloseApproachEvent> events,
	                                                          std::span<const std::optional<TraceEvaluation>> traces,
	                                                          std::span<const CloseApproachPassage> passages,
	                                                          double windowStartJd, double windowEndJd, double epochJd,
	                                                          double referenceJd) noexcept
	{
		const double start = std::min(windowStartJd, windowEndJd);
		const double end = std::max(windowStartJd, windowEndJd);
		const bool forward = referenceJd >= epochJd;

		const CloseApproachPassage* terminal = nullptr;
		for (const auto& passage : passages) {
			if (forward) {
				if (passage.lastInsideJd < epochJd || passage.firstInsideJd > end)
					continue;
				if (!terminal || passage.firstInsideJd > terminal->firstInsideJd)
					terminal = &passage;
			} else {
				if (passage.firstInsideJd > epochJd || passage.lastInsideJd < start)
					continue;
				if (!terminal || passage.lastInsideJd < terminal->lastInsideJd)
					terminal = &passage;
			}
		}
		if (!terminal)
			return std::nullopt;

		std::optional<std::size_t> lastMinimum;
		std::optional<std::size_t> selected;
		for (std::size_t i = 0; i < events.size(); ++i) {
			const auto& event = events[i];
			if (event.encounterIndex != terminal->encounterIndex || (forward ? event.jd < epochJd : event.jd > epochJd))
				continue;
			if (!lastMinimum || (forward && event.jd > events[*lastMinimum].jd) ||
			    (!forward && event.jd < events[*lastMinimum].jd))
				lastMinimum = i;
			if (!selected || event.distAu < events[*selected].distAu)
				selected = i;
		}

		if (!lastMinimum || events[*lastMinimum].jd < start || events[*lastMinimum].jd > end)
			return std::nullopt;
		if (!selected || *selected >= traces.size() || !traces[*selected] || !traces[*selected]->trace.valid)
			return std::nullopt;
		return selected;
	}

	std::vector<double> encounterCheckpointJds(const std::vector<CloseApproachEvent>& events, double epochJd,
	                                           double targetJd)
	{
		const bool forward = targetJd >= epochJd;
		std::vector<double> checkpoints;
		checkpoints.reserve(events.size());
		for (const CloseApproachEvent& event : events) {
			if (!std::isfinite(event.jd))
				continue;
			const bool isIntermediate =
			    forward ? event.jd > epochJd && event.jd < targetJd : event.jd < epochJd && event.jd > targetJd;
			if (isIntermediate)
				checkpoints.push_back(event.jd);
		}

		std::sort(checkpoints.begin(), checkpoints.end());
		checkpoints.erase(std::unique(checkpoints.begin(), checkpoints.end()), checkpoints.end());
		if (!forward)
			std::reverse(checkpoints.begin(), checkpoints.end());
		return checkpoints;
	}

}  // namespace orbital_core::hazard::lov_detail
