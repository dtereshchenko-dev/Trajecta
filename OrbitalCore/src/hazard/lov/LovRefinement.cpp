#include "hazard/lov/LovTraceDetail.h"

#include "hazard/HazardAudit.h"
#include "hazard/VirtualImpactorRefinement.h"
#include "orbital_core/hazard/PlanetData.h"
#include "logger/Logger.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <optional>
#include <utility>

namespace orbital_core::hazard::lov_detail {

	LovMinimumResult refineBracketedMinimum(const TraceRecord& left, const TraceRecord& right,
	                                        const TraceEvaluationContext& ctx)
	{
		audit::Scope auditScope(left.sampleIndex, left.trace.sigma, true);
		return refineBracketedMinimum(left, right, [&ctx](int sampleIndex, double sigma, double startJd, double endJd) {
			return evaluateSigmaNearEncounter(sampleIndex, sigma, startJd, endJd, ctx);
		});
	}

	OffLovCheckResult refineVirtualImpactorRecord(const TraceRecord& representative, const TraceEvaluationContext& ctx)
	{
		constexpr int kMaxIterations = 12;
		constexpr int kMaxFolds = 10;
		constexpr double kFoldSigma = 10.0;

		const ImpactRiskConfig& refinementCfg = ctx.cfg.impactRefinement;

		if (!refinementCfg.rmsScaling.linearization) {
			return {};
		}
		if (representative.parameters.size() < 6) {
			return {OffLovCheckStatus::Rejected, std::nullopt, "invalid-parameters"};
		}

		const double impactRadiusKm =
		    representative.trace.impactCrossSectionBodyRadii *
		    std::max(representative.trace.bodyRadiusKm, planetRadiusKm(representative.event.planet));
		if (!(impactRadiusKm > 0.0) || !std::isfinite(impactRadiusKm)) {
			return {OffLovCheckStatus::Rejected, std::nullopt, "invalid-impact-section"};
		}

		if (representative.trace.distanceKm <= 1.01 * impactRadiusKm) {
			return {OffLovCheckStatus::AlreadyInside, representative, "already-inside"};
		}

		const double qNom = nominalObservationChiSquare(ctx);
		if (!std::isfinite(qNom)) {
			return {OffLovCheckStatus::Rejected, std::nullopt, "invalid-nominal-objective"};
		}


		TraceRecord current = representative;
		const Eigen::Vector2d originalCoords = current.trace.coordsKm;
		const double originalSigmaLov = current.trace.sigma;
		const double originalSigmaImpact = current.probability.sigmaImpact;
		const double originalProbability = current.probability.probability;
		VecX baseParameters = representative.parameters;
		int folds = 0;

		for (int iter = 0; iter < kMaxIterations; ++iter) {
			const Eigen::Vector2d target = impactBoundaryTarget(current.trace.coordsKm, impactRadiusKm);


			VirtualImpactorRefinementInput input;
			input.cfg = &refinementCfg;
			input.dCoordsDParams = current.dCoordsDParams;
			input.parameterCovariance = current.covariance;
			input.currentCoordsKm = current.trace.coordsKm;
			input.targetCoordsKm = target;
			input.initialParameters = baseParameters;
			const VirtualImpactorRefinementResult refined = refineVirtualImpactorWithObservations(input);
			if (!refined.valid) {
				return {OffLovCheckStatus::Rejected, std::nullopt,
				        refined.badConditioning ? "bad-conditioning" : "semilinear-step-failed"};
			}

			std::optional<TraceRecord> next =
			    evaluateParametersNearEncounter(current, refined.parameters, current.event.jd, ctx, 50.0);
			if (!next) {
				return {OffLovCheckStatus::Rejected, std::nullopt, "no-target-plane"};
			}
			next->covariance = refined.parameterCovariance;

			const double deltaQ = std::max(0.0, refined.observationChiSquare - qNom);
			const double objectiveSigma = std::sqrt(deltaQ);

			const Eigen::Vector2d beforeMove = target - current.trace.coordsKm;
			const Eigen::Vector2d afterMove = next->trace.coordsKm - current.trace.coordsKm;
			const bool folded = beforeMove.dot(afterMove) < 0.0 || objectiveSigma > kFoldSigma;
			if (folded)
				++folds;
			if (folds > kMaxFolds) {
				return {OffLovCheckStatus::Rejected, std::nullopt, "fold"};
			}

			if (next->trace.distanceKm <= 1.01 * impactRadiusKm) {
				const double exponent =
				    -0.5 * (deltaQ - originalSigmaLov * originalSigmaLov - originalSigmaImpact * originalSigmaImpact);
				const double factor = std::exp(exponent);
				if (!std::isfinite(factor) || factor <= 0.0) {
					return {OffLovCheckStatus::Rejected, std::nullopt, "invalid-probability-correction"};
				}

				const double correctedProbability = originalProbability * factor;
				audit::write("off-lov-probability", [&](auto& log) {
					log << " iteration=" << iter + 1 << " jd=" << next->event.jd << " qNom=" << qNom
					    << " qObs=" << refined.observationChiSquare << " deltaQ=" << deltaQ
					    << " sigmaLov=" << originalSigmaLov << " sigmaImpact=" << originalSigmaImpact
					    << " originalP=" << originalProbability << " factor=" << factor
					    << " correctedP=" << correctedProbability;
				});
				if (correctedProbability < kRiskProbabilityFloor || correctedProbability > 1.0) {
					return {OffLovCheckStatus::Rejected, std::nullopt, "corrected-probability-out-of-range"};
				}

				next->probability.valid = true;
				next->probability.candidate = true;
				next->probability.sigmaImpact = originalSigmaImpact;
				next->probability.probability = correctedProbability;
				next->riskCandidate = true;

				const double objScale = (std::isfinite(ctx.cfg.impactRefinement.rmsScaling.objectiveScale) &&
				                         ctx.cfg.impactRefinement.rmsScaling.objectiveScale > 0.0)
				                            ? ctx.cfg.impactRefinement.rmsScaling.objectiveScale
				                            : 1.0;

				LOG_INFO << std::fixed << std::setprecision(6) << "[LineOfVariation] VI refinement accepted"
				         << "  iter=" << iter + 1 << "  jd=" << next->event.jd << "  objectiveScale=" << objScale
				         << "  qNomScaled=" << qNom << "  qNomRaw=" << qNom * objScale
				         << "  qObsScaled=" << refined.observationChiSquare
				         << "  qObsRaw=" << refined.observationChiSquare * objScale << "  factor=" << factor
				         << "  p=" << std::scientific << std::setprecision(12) << next->probability.probability
				         << std::fixed << std::setprecision(6) << "  missKm=" << next->trace.distanceKm
				         << "  sigmaLov=" << originalSigmaLov << "  sigmaObjectiveScaled=" << objectiveSigma
				         << "  sigmaObjectiveRaw=" << std::sqrt(deltaQ * objScale)
				         << "  shiftKm=" << (next->trace.coordsKm - originalCoords).norm();
				return {OffLovCheckStatus::Refined, std::move(next), "refined"};
			}

			if (!folded)
				baseParameters = refined.parameters;
			current = std::move(*next);
		}

		return {OffLovCheckStatus::Rejected, std::nullopt, "iteration-limit"};
	}

	LovMinimumResult refineEndpointMinimum(const TraceRecord& start, double sigmaStepLimit,
	                                       const TraceEvaluationContext& ctx)
	{
		audit::Scope auditScope(start.sampleIndex, start.trace.sigma, true);
		return refineEndpointMinimum(start, sigmaStepLimit,
		                             [&ctx](int sampleIndex, double sigma, double startJd, double endJd) {
			                             return evaluateSigmaNearEncounter(sampleIndex, sigma, startJd, endJd, ctx);
		                             });
	}


}  // namespace orbital_core::hazard::lov_detail
