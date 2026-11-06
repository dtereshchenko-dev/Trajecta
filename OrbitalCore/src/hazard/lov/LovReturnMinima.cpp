#include "hazard/lov/LovTraceDetail.h"
#include <vector>
#include "hazard/HazardAudit.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

namespace orbital_core::hazard::lov_detail {

	std::vector<TraceRecord> addRefinedMinima(
	    const std::vector<TraceRecord>& records, double sigmaMax,
	    const std::function<double(const TraceRecord&, bool)>& endpointLimit,
	    const std::function<LovMinimumResult(const TraceRecord&, double)>& refineEndpoint,
	    const std::function<LovMinimumResult(const TraceRecord&, const TraceRecord&)>& refineBracket)
	{
		constexpr double kNewtonDistanceBodyRadii = 60.0;
		if (records.empty())
			return {};

		std::vector<TraceRecord> refined;
		const auto runEndpoint = [&](const TraceRecord& record, double stepLimit, bool requireConvergence = false) {
			LovMinimumResult result = refineEndpoint(record, stepLimit);
			if (result.usable(requireConvergence))
				refined.push_back(*result.record);
			return result;
		};
		const auto tryHeadOrTail = [&](const TraceRecord& record, bool increasingHead) {
			if (minimumPossibleBodyRadii(record, sigmaMax) < kNewtonDistanceBodyRadii)
				runEndpoint(record, endpointLimit(record, increasingHead));
		};

		if (records.front().trace.dDistance2DSigmaBodyRadii2 > 0.0)
			tryHeadOrTail(records.front(), true);
		if (records.back().trace.dDistance2DSigmaBodyRadii2 < 0.0)
			tryHeadOrTail(records.back(), false);

		for (std::size_t i = 1; i < records.size(); ++i) {
			const TraceRecord& left = records[i - 1];
			const TraceRecord& right = records[i];
			const double minPossible =
			    std::min(minimumPossibleBodyRadii(left, sigmaMax), minimumPossibleBodyRadii(right, sigmaMax));
			if (minPossible >= kNewtonDistanceBodyRadii)
				continue;

			const LovPairType type = classifyLovPair(left, right);
			const double stepLimit = 0.5 * std::abs(right.trace.sigma - left.trace.sigma);
			switch (type) {
				case LovPairType::SimpleMinimum:
				case LovPairType::InterruptedMinimum: {
					LovMinimumResult found = refineBracket(left, right);
					bool reclassified = false;
					if (found.usable()) {
						refined.push_back(*found.record);
					} else {
						const LovMinimumResult fromLeft = runEndpoint(left, stepLimit, true);
						constexpr double kUniformProbabilitySamplingStep = 0.0;
						found = runEndpoint(right, kUniformProbabilitySamplingStep / 2.0, true);
						reclassified = fromLeft.reclassified || found.reclassified;
					}

					if (type != LovPairType::InterruptedMinimum || reclassified || !found.record) {
						break;
					}
					constexpr double kAlmostInterruptedCos = 0.70;
					constexpr double kDistanceTolerance = 2.3e-3;
					const TraceRecord& anchor = *found.record;
					const double leftAngle = angleDifferenceAbs(left.trace.alphaLovRad, anchor.trace.alphaLovRad);
					const double rightAngle = angleDifferenceAbs(right.trace.alphaLovRad, anchor.trace.alphaLovRad);
					const double foundDistance = anchor.trace.distanceBodyRadii;

					audit::Scope auditScope(left.sampleIndex, anchor.trace.sigma);
					audit::write("interrupted-extra-choice", [&](auto& log) {
						log << " rightSample=" << right.sampleIndex << " leftAlpha=" << left.trace.alphaLovRad
						    << " rightAlpha=" << right.trace.alphaLovRad << " anchorAlpha=" << anchor.trace.alphaLovRad
						    << " leftCos=" << std::cos(leftAngle) << " rightCos=" << std::cos(rightAngle)
						    << " threshold=" << kAlmostInterruptedCos << " chosen="
						    << (std::cos(rightAngle) > kAlmostInterruptedCos  ? "left"
						        : std::cos(leftAngle) > kAlmostInterruptedCos ? "right"
						                                                      : "none");
					});

					LovMinimumResult extra;
					if (std::cos(rightAngle) > kAlmostInterruptedCos) {
						const double limit = 0.5 * std::abs(left.trace.sigma - anchor.trace.sigma);
						extra = refineEndpoint(left, limit);
					} else if (std::cos(leftAngle) > kAlmostInterruptedCos) {
						const double limit = 0.5 * std::abs(right.trace.sigma - anchor.trace.sigma);
						extra = refineEndpoint(right, limit);
					}
					if (extra.usable(true) && std::abs(extra.record->trace.distanceBodyRadii - foundDistance) >
					                              kDistanceTolerance * foundDistance) {
						refined.push_back(std::move(*extra.record));
					}
					break;
				}
				case LovPairType::InterruptedUnknown:
					if (left.trace.dDistance2DSigmaBodyRadii2 < 0.0)
						runEndpoint(left, stepLimit);
					else if (right.trace.dDistance2DSigmaBodyRadii2 > 0.0)
						runEndpoint(right, stepLimit);
					break;
				case LovPairType::EntangledDescending:
					runEndpoint(left, stepLimit);
					break;
				case LovPairType::EntangledAscending:
					runEndpoint(right, stepLimit);
					break;
				case LovPairType::NoConcern:
				case LovPairType::SimpleMaximum:
				case LovPairType::InterruptedMaximum:
				case LovPairType::Unknown:
					break;
			}
		}
		return refined;
	}

}  // namespace orbital_core::hazard::lov_detail
