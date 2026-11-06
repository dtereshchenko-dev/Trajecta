#include "hazard/lov/LovTraceDetail.h"

#include "orbital_core/Constants.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

namespace orbital_core::hazard::lov_detail {

	namespace {
		BPlaneCovariance analyticTargetPlaneCovariance(const TraceRecord& representative)
		{
			const TargetPlaneTracePoint& trace = representative.trace;
			BPlaneCovariance cov;
			cov.meanXi = trace.coordsKm.x();
			cov.meanZeta = trace.coordsKm.y();
			cov.numSamples = 1;
			cov.analytic = true;

			if (!(trace.stretchKm > 0.0) || !(trace.widthKm > 0.0) || !std::isfinite(trace.alphaRad)) {
				return cov;
			}

			const double sinAlpha = std::sin(trace.alphaRad);
			const double cosAlpha = std::cos(trace.alphaRad);
			const double major2 = trace.stretchKm * trace.stretchKm;
			const double minor2 = trace.widthKm * trace.widthKm;
			const double varXi = major2 * sinAlpha * sinAlpha + minor2 * cosAlpha * cosAlpha;
			const double varZeta = major2 * cosAlpha * cosAlpha + minor2 * sinAlpha * sinAlpha;
			const double covXiZeta = (minor2 - major2) * sinAlpha * cosAlpha;

			cov.sigmaXi = std::sqrt(std::max(0.0, varXi));
			cov.sigmaZeta = std::sqrt(std::max(0.0, varZeta));
			const double sigmaProduct = cov.sigmaXi * cov.sigmaZeta;
			cov.correlation = sigmaProduct > 0.0 ? std::clamp(covXiZeta / sigmaProduct, -1.0, 1.0) : 0.0;
			return cov;
		}
	}  // namespace

	static BPlaneCovariance fitScatterCovariance(const std::vector<TraceRecord>& records)
	{
		BPlaneCovariance cov;
		cov.numSamples = static_cast<int>(records.size());
		if (records.size() < 2)
			return cov;

		for (const TraceRecord& r : records) {
			cov.meanXi += r.trace.coordsKm.x();
			cov.meanZeta += r.trace.coordsKm.y();
		}
		cov.meanXi /= static_cast<double>(records.size());
		cov.meanZeta /= static_cast<double>(records.size());

		double varXi = 0.0;
		double varZeta = 0.0;
		double covXiZeta = 0.0;
		for (const TraceRecord& r : records) {
			const double dx = r.trace.coordsKm.x() - cov.meanXi;
			const double dy = r.trace.coordsKm.y() - cov.meanZeta;
			varXi += dx * dx;
			varZeta += dy * dy;
			covXiZeta += dx * dy;
		}

		const double denom = static_cast<double>(records.size() - 1);
		varXi /= denom;
		varZeta /= denom;
		covXiZeta /= denom;

		cov.sigmaXi = std::sqrt(std::max(0.0, varXi));
		cov.sigmaZeta = std::sqrt(std::max(0.0, varZeta));
		const double sigmaProduct = cov.sigmaXi * cov.sigmaZeta;
		cov.correlation = sigmaProduct > 0.0 ? std::clamp(covXiZeta / sigmaProduct, -1.0, 1.0) : 0.0;
		return cov;
	}

	ImpactRiskRow makeRow(std::size_t clusterIndex, const std::vector<TraceRecord>& records,
	                      const TraceRecord& representative, const FitParameterCovariance& covariance)
	{
		ImpactRiskRow row;
		row.clusterIndex = clusterIndex;
		row.representativeSampleIndex = representative.sampleIndex;
		row.representativeEvent = representative.event;
		row.representativeParameters = representative.parameters;
		row.cluster.startJd = records.front().event.jd;
		row.cluster.endJd = records.front().event.jd;
		row.cluster.meanJd = 0.0;
		row.cluster.minDistAu = std::numeric_limits<double>::infinity();
		row.cluster.firstEncounterIndex = 0;
		row.cluster.lastEncounterIndex = static_cast<int>(records.size());
		row.cluster.representativeEncounterIndex = representative.sampleIndex;

		std::unordered_set<int> sampleIds;
		std::vector<double> misses;
		misses.reserve(records.size());

		for (const TraceRecord& r : records) {
			row.cluster.startJd = std::min(row.cluster.startJd, r.event.jd);
			row.cluster.endJd = std::max(row.cluster.endJd, r.event.jd);
			row.cluster.meanJd += r.event.jd;
			row.cluster.minDistAu = std::min(row.cluster.minDistAu, r.event.distAu);
			sampleIds.insert(r.sampleIndex);
			row.result.scatterKm.push_back({r.trace.coordsKm.x(), r.trace.coordsKm.y()});
			misses.push_back(r.trace.distanceKm);
		}

		row.cluster.numEncounters = static_cast<int>(records.size());
		row.cluster.numSamples = static_cast<int>(sampleIds.size());
		row.cluster.meanJd /= static_cast<double>(records.size());
		row.cluster.minDistLd = row.cluster.minDistAu * kAuKm / Constants::LUNAR_DISTANCE_KM;

		std::sort(misses.begin(), misses.end());

		ImpactRiskResult& r = row.result;
		r.status = ImpactRiskStatus::Completed;
		r.probabilitySource = ImpactProbabilitySource::LineOfVariation;
		r.probabilityImpact = representative.probability.probability;
		r.sigmaProbability = representative.probability.sigmaImpact;
		r.probabilityBeforeRmsScaling = r.probabilityImpact;
		r.sigmaBeforeRmsScaling = r.sigmaProbability;
		r.sigmaLov = representative.trace.sigma;
		r.computed = representative.probability.valid;
		r.converged = representative.trace.valid;
		r.gamma = representative.gamma;
		r.vInfKmS = representative.event.vInfKmS;
		r.bPlaneVIKm = representative.trace.coordsKm;
		r.basis = representative.basis;
		r.scatterCov = fitScatterCovariance(records);
		if (!r.scatterCov.isValid())
			r.scatterCov = analyticTargetPlaneCovariance(representative);
		r.sigmaAlongMissKm = representative.trace.stretchKm;
		r.sigmaTransverseKm = representative.trace.widthKm;
		r.minMissKm = representative.trace.distanceKm;
		r.medianMissKm = misses.empty() ? representative.trace.distanceKm : misses[misses.size() / 2];
		r.bestMatchTcaJd = representative.trace.jd;
		r.numImpactHits = static_cast<int>(
		    std::count_if(records.begin(), records.end(), [](const TraceRecord& tr) { return tr.trace.impact; }));
		r.parameterVI = representative.parameters;
		r.parameterStepPerSigmaVI = representative.stepPerSigma;
		r.parameterCovarianceVI = covariance;
		return row;
	}

}  // namespace orbital_core::hazard::lov_detail
