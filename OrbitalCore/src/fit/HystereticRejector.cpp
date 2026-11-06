#include "HystereticRejector.h"
#include "OppositionCoverage.h"
#include <span>
#include <utility>
#include <vector>
#include <cstdint>

#include "logger/Logger.h"
#include "orbital_core/Types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Calculations {
	namespace {

		double includedTimeSpanDays(const std::vector<std::uint8_t>& mask, std::span<const ObsNode> allObs)
		{
			double jdMin = std::numeric_limits<double>::infinity();
			double jdMax = -std::numeric_limits<double>::infinity();
			std::size_t nIncluded = 0;
			for (std::size_t i = 0; i < allObs.size(); ++i) {
				if (!mask[i])
					continue;
				jdMin = std::min(jdMin, allObs[i].JD);
				jdMax = std::max(jdMax, allObs[i].JD);
				++nIncluded;
			}
			return (nIncluded > 0) ? (jdMax - jdMin) : 0.0;
		}

		std::size_t countIncluded(const std::vector<std::uint8_t>& mask)
		{
			std::size_t n = 0;
			for (auto v : mask)
				if (v)
					++n;
			return n;
		}

	}  // namespace

	HystereticRejectResult applyHystereticReject(std::span<const ObsNode> allObs, std::vector<std::uint8_t>& included,
	                                             std::span<const double> chiSquared, const HystereticRejectConfig& cfg,
	                                             std::span<const std::uint8_t> locked)
	{
		if (!std::isfinite(cfg.rejectChiSquared) || !std::isfinite(cfg.recoverChiSquared) ||
		    cfg.recoverChiSquared <= 0.0 || cfg.rejectChiSquared < cfg.recoverChiSquared ||
		    !std::isfinite(cfg.retryScale) || cfg.retryScale <= 1.0 || !std::isfinite(cfg.maxChiSquared) ||
		    cfg.maxChiSquared <= cfg.rejectChiSquared)
			throw std::invalid_argument("Invalid chi-squared rejection thresholds");
		HystereticRejectResult res;
		res.nKeptBefore = static_cast<int>(countIncluded(included));

		res.nKept = res.nKeptBefore;
		if (allObs.empty() || included.size() != allObs.size() || chiSquared.size() != allObs.size() ||
		    (!locked.empty() && locked.size() != allObs.size())) {
			LOG_WARN << "[HystereticRejector] invalid post-fit residual field";
			return res;
		}

		std::vector<std::uint8_t> eligible(allObs.size(), 1);
		for (std::size_t i = 0; i < locked.size(); ++i)
			if (locked[i] && !included[i])
				eligible[i] = 0;
		const auto eligibleCount = countIncluded(eligible);
		const double origSpan = includedTimeSpanDays(eligible, allObs);
		const orbital_core::fit::OppositionCoverage originalCoverage(allObs, included);
		double largestChiSquared = 0.0;
		for (std::size_t i = 0; i < allObs.size(); ++i)
			if (included[i] && !allObs[i].radar && (locked.empty() || !locked[i]) && !originalCoverage.isLast(i) &&
			    std::isfinite(chiSquared[i]))
				largestChiSquared = std::max(largestChiSquared, chiSquared[i]);
		const double batchThreshold = largestChiSquared * (cfg.rejectLargestOnly ? 1.0 : 0.25);

		double rejectLimit = cfg.rejectChiSquared;
		double recoverLimit = cfg.recoverChiSquared;
		while (rejectLimit < cfg.maxChiSquared) {

			std::vector<std::uint8_t> origMask = included;
			auto coverage = originalCoverage;
			int flippedIn = 0;
			int flippedOut = 0;
			for (std::size_t i = 0; i < allObs.size(); ++i) {
				if (!locked.empty() && locked[i])
					continue;
				if (!std::isfinite(chiSquared[i]) || chiSquared[i] < 0.0)
					continue;
				const double r2 = chiSquared[i];
				if (!included[i] && r2 <= recoverLimit) {
					included[i] = 1;
					coverage.include(i);
					++flippedIn;
					LOG_DEBUG << "[HystereticRejector] IN  idx=" << i << " jd=" << allObs[i].JD
					          << " date=" << Convertors::jdToPackedDate(allObs[i].JD) << " code=" << allObs[i].obsCode
					          << " postFitChi2=" << r2 << " recoverChi2=" << recoverLimit;
				} else if (included[i] && !allObs[i].radar && r2 > rejectLimit && r2 >= batchThreshold) {
					if (coverage.isLast(i)) {
						LOG_WARN << "[HystereticRejector] keeping last observation in opposition: idx=" << i
						         << " date=" << Convertors::jdToPackedDate(allObs[i].JD)
						         << " code=" << allObs[i].obsCode << " postFitChi2=" << r2;
						continue;
					}
					included[i] = 0;
					coverage.exclude(i);
					++flippedOut;
					LOG_DEBUG << "[HystereticRejector] OUT idx=" << i << " jd=" << allObs[i].JD
					          << " date=" << Convertors::jdToPackedDate(allObs[i].JD) << " code=" << allObs[i].obsCode
					          << " postFitChi2=" << r2 << " rejectChi2=" << rejectLimit;
				}
			}

			const std::size_t nKept = countIncluded(included);
			const double tSpan = includedTimeSpanDays(included, allObs);
			const std::size_t nMin =
			    static_cast<std::size_t>(std::floor(static_cast<double>(eligibleCount) * cfg.minKeptFraction));
			const double tMin = origSpan * cfg.minSpanFraction;
			const bool okN = (nKept >= nMin);
			const bool okT = (tSpan >= tMin);

			if (okN && okT) {
				res.appliedRejectChiSquared = rejectLimit;
				res.appliedRecoverChiSquared = recoverLimit;
				res.nFlippedIn = flippedIn;
				res.nFlippedOut = flippedOut;
				res.nKept = static_cast<int>(nKept);
				res.status = (flippedIn + flippedOut > 0) ? RejectStatus::Changed : RejectStatus::NoChange;

				LOG_INFO << "[HystereticRejector] rejectChi2=" << rejectLimit << " recoverChi2=" << recoverLimit
				         << " batchChi2=" << batchThreshold << ":"
				         << " in=" << flippedIn << " out=" << flippedOut << " kept=" << nKept << "/" << allObs.size()
				         << " span=" << tSpan << "d/" << origSpan << "d"
				         << " status=" << (res.status == RejectStatus::Changed ? "CHANGED" : "NO_CHANGE");
				return res;
			}

			LOG_INFO << "[HystereticRejector] rejectChi2=" << rejectLimit << " safeguard fail: kept=" << nKept
			         << " (min " << nMin << ")"
			         << " span=" << tSpan << "d (min " << tMin << "d)"
			         << " — rolling back, escalating";
			included = std::move(origMask);
			rejectLimit *= cfg.retryScale;
			recoverLimit *= cfg.retryScale;
		}

		res.appliedRejectChiSquared = rejectLimit;
		res.appliedRecoverChiSquared = recoverLimit;
		res.nKept = static_cast<int>(countIncluded(included));
		res.status = RejectStatus::SafeguardRollback;
		LOG_INFO << "[HystereticRejector] escalation exhausted at chi2=" << rejectLimit
		         << " — no rejection passes safeguards, leaving mask unchanged";
		return res;
	}

}  // namespace Calculations
