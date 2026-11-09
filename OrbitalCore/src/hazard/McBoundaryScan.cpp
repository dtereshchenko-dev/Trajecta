#include "McBoundaryScan.h"
#include "BoundaryScan.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core::hazard::mc {
	void scan(McClone& clone, const IntervalMcConfig& config, NBodyPropagator& propagator, IEphemerisSource& ephemeris,
	          std::stop_token cancel)
	{
		const auto result = boundary::search(clone.state, clone.epochMjdTdb, clone.nonGrav, config.target,
		                                     config.startMjdTdb, config.endMjdTdb, config.altitudeKm, propagator,
		                                     ephemeris, cancel);
		clone.entryBeforeInterval = result.beforeInterval;
		clone.entry = result.entry;
	}

	void summarize(IntervalMcResult& result)
	{
		std::vector<double> epochs;
		for (const auto& clone : result.clones) {
			switch (clone.status) {
				case McCloneStatus::Complete: {
					++result.completed;
					result.entriesBeforeInterval += clone.entryBeforeInterval;
					if (!clone.entry)
						break;
					++result.entryClones;
					const double epoch = clone.entry->mjdTdb;
					if (!std::isfinite(epoch))
						throw std::runtime_error("Non-finite MC entry epoch");
					epochs.push_back(epoch);
					break;
				}
				case McCloneStatus::FitFailed:
					++result.failedFit;
					break;
				case McCloneStatus::PropagationFailed:
					++result.failedPropagation;
					break;
				case McCloneStatus::Cancelled:
					++result.cancelledClones;
					break;
				case McCloneStatus::Unstarted:
					++result.unstarted;
					break;
			}
		}
		std::sort(epochs.begin(), epochs.end());
		for (const double epoch : epochs) {
			if (result.events.empty() || epoch - result.events.back().endMjdTdb > eventMaxGapDays)
				result.events.push_back({.startMjdTdb = epoch, .endMjdTdb = epoch});
			auto& row = result.events.back();
			row.endMjdTdb = epoch;
			++row.entryClones;
		}
		for (auto& row : result.events) {
			const double n = result.completed, p = row.entryClones / n, z = 1.959963984540054;
			const double denominator = 1 + z * z / n;
			const double center = (p + z * z / (2 * n)) / denominator;
			const double half = z * std::sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator;
			row.probability = p;
			row.probabilityLow = std::max(0.0, center - half);
			row.probabilityHigh = std::min(1.0, center + half);
		}
	}
}  // namespace orbital_core::hazard::mc
