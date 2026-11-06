#include "hazard/lov/LovTraceDetail.h"
#include <utility>

#include "hazard/LovReturnAnalysis.h"
#include "orbital_core/hazard/PlanetData.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace orbital_core::hazard::lov_detail {
	namespace {
		[[nodiscard]] bool isPhysicalImpact(const TraceRecord& record) noexcept
		{
			const double radiusKm = planetRadiusKm(record.event.planet);
			return radiusKm > 0.0 && std::isfinite(record.event.distAu) && record.event.distAu >= 0.0 &&
			       record.event.distAu * kAuKm < radiusKm;
		}
	}  // namespace

	bool isReportableRisk(const TraceRecord& record) noexcept
	{
		return record.riskCandidate && record.probability.valid && record.probability.candidate &&
		       record.probability.probability > kRiskProbabilityFloor;
	}

	int countImpacts(const std::vector<TraceRecord>& records) noexcept
	{
		return static_cast<int>(std::count_if(records.begin(), records.end(), isPhysicalImpact));
	}

	double sampleDeltaSigma(const std::vector<double>& sigmas, int sampleIndex) noexcept
	{
		if (sampleIndex <= 0 || sampleIndex >= static_cast<int>(sigmas.size()))
			return 0.0;

		const double delta =
		    sigmas[static_cast<std::size_t>(sampleIndex)] - sigmas[static_cast<std::size_t>(sampleIndex - 1)];
		return std::isfinite(delta) && delta > 0.0 ? delta : 0.0;
	}

	std::vector<TraceRecord> addRefinedMinima(const std::vector<TraceRecord>& recordsIn,
	                                          const TraceEvaluationContext& ctx)
	{
		const auto endpointLimit = [&](const TraceRecord& record, bool increasingHead) {
			static const std::vector<LovCatalogPoint> emptyCatalog;
			const auto& catalog = ctx.catalog ? *ctx.catalog : emptyCatalog;
			return endpointSigmaStepLimit(catalog, record.sampleIndex, increasingHead, ctx.cfg.sigmaStepMax);
		};
		return addRefinedMinima(
		    recordsIn, ctx.cfg.sigmaMax, endpointLimit,
		    [&ctx](const TraceRecord& record, double stepLimit) {
			    return refineEndpointMinimum(record, stepLimit, ctx);
		    },
		    [&ctx](const TraceRecord& left, const TraceRecord& right) {
			    return refineBracketedMinimum(left, right, ctx);
		    });
	}

	std::vector<std::vector<TraceRecord>> buildReturnRecords(std::vector<TraceRecord> records)
	{
		std::vector<std::vector<TraceRecord>> returns;
		if (records.empty())
			return returns;

		std::vector<LovTraceEntry> entries;
		entries.reserve(records.size());
		for (int i = 0; i < static_cast<int>(records.size()); ++i) {
			const TraceRecord& record = records[static_cast<std::size_t>(i)];
			entries.push_back({i, record.sampleIndex, record.event.jd});
		}

		const std::vector<LovReturnFilament> filaments = buildLovReturnFilaments(entries, LovReturnAnalysisConfig{});

		returns.reserve(filaments.size());
		for (const LovReturnFilament& filament : filaments) {
			if (filament.entryIndices.empty())
				continue;

			std::vector<TraceRecord> ret;
			ret.reserve(filament.entryIndices.size());
			for (const int entryIndex : filament.entryIndices) {
				const LovTraceEntry& entry = entries[static_cast<std::size_t>(entryIndex)];
				ret.push_back(std::move(records[static_cast<std::size_t>(entry.recordIndex)]));
			}

			std::sort(ret.begin(), ret.end(), [](const TraceRecord& a, const TraceRecord& b) {
				if (a.sampleIndex != b.sampleIndex)
					return a.sampleIndex < b.sampleIndex;
				return a.event.jd < b.event.jd;
			});
			returns.push_back(std::move(ret));
		}
		return returns;
	}

	std::vector<TraceRecord> addDenseImpactSegments(const std::vector<TraceRecord>& records,
	                                                const std::vector<double>& sigmas)
	{
		constexpr int kMinimumSegmentPoints = 10;
		if (records.empty())
			return {};

		std::vector<TraceRecord> denseRepresentatives;

		std::size_t begin = 0;
		while (begin < records.size()) {
			while (begin < records.size() && !isPhysicalImpact(records[begin])) {
				++begin;
			}
			if (begin >= records.size())
				break;

			std::size_t end = begin + 1;
			while (end < records.size() && isPhysicalImpact(records[end]) &&
			       records[end].sampleIndex == records[end - 1].sampleIndex + 1) {
				++end;
			}

			const std::size_t count = end - begin;
			if (count >= kMinimumSegmentPoints) {
				const std::size_t extendedBegin = begin > 0 ? begin - 1 : begin;
				const std::size_t extendedEnd = end < records.size() ? end + 1 : end;

				auto minIt = std::min_element(records.begin() + static_cast<std::ptrdiff_t>(extendedBegin),
				                              records.begin() + static_cast<std::ptrdiff_t>(extendedEnd),
				                              [](const TraceRecord& a, const TraceRecord& b) {
					                              return a.trace.distanceBodyRadii < b.trace.distanceBodyRadii;
				                              });

				double p = 0.0;
				for (std::size_t i = begin; i < end; ++i) {
					const TraceRecord& r = records[i];
					const double delta = sampleDeltaSigma(sigmas, r.sampleIndex);
					if (delta <= 0.0)
						continue;
					p += delta * std::exp(-0.5 * r.trace.sigma * r.trace.sigma) * kInvSqrtTwoPi;
				}

				if (p > 0.0) {
					TraceRecord representative = *minIt;
					const std::size_t outerLeft = extendedBegin;
					const std::size_t outerRight = extendedEnd - 1;
					const double sigmaSpan = std::abs(records[outerRight].trace.sigma - records[outerLeft].trace.sigma);
					if (sigmaSpan > 0.0) {
						const double stretchKm =
						    (records[outerRight].trace.coordsKm - records[outerLeft].trace.coordsKm).norm() / sigmaSpan;
						representative.trace.stretchKm = stretchKm;
						representative.trace.stretchBodyRadii = stretchKm / representative.trace.bodyRadiusKm;
					}

					const double widthKm = 0.5 * (records[outerLeft].trace.widthKm + records[outerRight].trace.widthKm);
					representative.trace.widthKm = widthKm;
					representative.trace.widthBodyRadii = widthKm / representative.trace.bodyRadiusKm;
					representative.probability.valid = true;
					representative.probability.candidate = true;
					representative.probability.sigmaImpact = 0.0;
					representative.probability.probability = std::clamp(p, 0.0, 1.0);
					representative.riskCandidate = true;
					denseRepresentatives.push_back(std::move(representative));
				}
			}

			begin = end;
		}

		return denseRepresentatives;
	}

}  // namespace orbital_core::hazard::lov_detail
