#include "orbital_core/FitSession.h"
#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "fit/FitCancellation.h"
#include "fit/FitPassRunner.h"
#include "fit/FitSessionRunner.h"
#include "fit/FitObservations.h"
#include "fit/HystereticRejector.h"
#include "fit/ObservationSelection.h"
#include "logger/Logger.h"

namespace orbital_core {
	namespace {
		struct FitRun {
			const FitConfig& cfg;
			const FitRunControl& control;
			double outputEpochJdTdb;
			std::span<const obs::PreparedObservation> suppliedPreparation{};

			[[nodiscard]] FitResult run(const FitRequest& request) const;

			struct KeptObservations {
				std::vector<ObsNode> observations;
				std::vector<std::size_t> indices;
			};
		};
	}  // namespace

	FitSession::FitSession(FitConfig cfg) : config_(std::move(cfg)) {}

	FitResult FitRun::run(const FitRequest& request) const
	try {
		const auto cancel = control.cancel;
		fit::checkCancellation(cancel);
		const auto& initial = request.initialOrbit;
		const auto observations = request.observations;
		const auto initialIncluded = request.initialIncluded;
		const auto* warmStartNg = request.warmStartNonGrav;
		if (!cfg.nonGrav.enabled || !warmStartNg || !warmStartNg->enabled() ||
		    warmStartNg->model() != static_cast<NgModel>(cfg.nonGrav.model))
			warmStartNg = nullptr;
		const auto& options = request.options;
		if (!initial)
			throw std::invalid_argument("FitSession::optimize: initial orbit is null");

		if (observations.size() < 3)
			throw std::invalid_argument("FitSession::optimize: need at least 3 observations");

		if (request.outputEpochJdTdb && (!std::isfinite(*request.outputEpochJdTdb) || *request.outputEpochJdTdb <= 0.0))
			throw std::invalid_argument("FitSession::optimize: output epoch must be a finite positive JD TDB");

		const auto& clipCfg = cfg.sigmaClip;
		if (clipCfg.enabled && (!std::isfinite(clipCfg.thresholdNSigma) || clipCfg.thresholdNSigma <= 0.0 ||
		                        clipCfg.thresholdNSigma >= 60.0))
			throw std::invalid_argument("Invalid sigma rejection threshold");

		auto selection = fit::resolveObservationSelection(observations, initialIncluded, request.observationOverrides);
		auto& included = selection.included;

		LOG_INFO << "[FitSession] START observations=" << observations.size()
		         << " selection=" << (options.selection == FitSelectionMode::Fixed ? "fixed" : "configured")
		         << " evaluateOnly=" << options.evaluateOnly << " clipping=" << clipCfg.enabled
		         << " sigma=" << clipCfg.thresholdNSigma << " nonGrav=" << cfg.nonGrav.enabled
		         << " model=" << static_cast<int>(cfg.nonGrav.model)
		         << " solveMask=" << static_cast<int>(cfg.nonGrav.solveMask)
		         << " warmNonGrav=" << (warmStartNg != nullptr)
		         << " manualOverrides=" << request.observationOverrides.size();

		fit::FitPassRunner passes{cfg, outputEpochJdTdb, control};
		std::vector<obs::PreparedObservation> ownedPreparation;
		auto prepared = suppliedPreparation;
		if (prepared.empty()) {
			ownedPreparation = passes.prepareObservations(observations);
			prepared = ownedPreparation;
		}
		if (prepared.size() != observations.size())
			throw std::invalid_argument("Prepared observation count does not match fit input");

		auto buildKept = [&]() {
			KeptObservations k;
			k.observations.reserve(observations.size());
			k.indices.reserve(observations.size());
			for (std::size_t i = 0; i < observations.size(); ++i) {
				if (included[i]) {
					k.observations.push_back(observations[i]);
					k.indices.push_back(i);
				}
			}
			return k;
		};

		auto kept = buildKept();

		if (options.evaluateOnly || options.selection == FitSelectionMode::Fixed) {
			FitResult result = passes.run({.initial = initial,
			                               .observations = kept.observations,
			                               .prepared = prepared,
			                               .warmStartNonGrav = warmStartNg,
			                               .uncertainty = options.uncertainty,
			                               .evaluateOnly = options.evaluateOnly,
			                               .included = &included});
			result.keptObsIndices = kept.indices;
			result.sigmaClipRejected = static_cast<int>(observations.size() - kept.observations.size());
			return result;
		}

		if (clipCfg.enabled && clipCfg.maxPasses > 0) {
			Calculations::HystereticRejectConfig reject;
			reject.rejectChiSquared = clipCfg.thresholdNSigma * clipCfg.thresholdNSigma;
			reject.recoverChiSquared = reject.rejectChiSquared * Calculations::HystereticRejectConfig::recoveryRatio;
			const int maxChanges = clipCfg.maxPasses;
			int changes = 0;
			int previousModifications = -1;
			int repeatedModifications = 0;
			FitSelectionStatus status = FitSelectionStatus::NotApplied;
			const auto update = [&](std::vector<std::uint8_t>& mask, const std::vector<double>& chi,
			                        int parameterCount) {
				fit::checkCancellation(cancel);
				if (status == FitSelectionStatus::PassLimit || status == FitSelectionStatus::InsufficientObservations ||
				    status == FitSelectionStatus::SafeguardRollback || status == FitSelectionStatus::CorrectionFailed)
					return false;
				if (std::none_of(chi.begin(), chi.end(), [](double value) { return std::isfinite(value); })) {
					status = FitSelectionStatus::CorrectionFailed;
					return false;
				}
				const auto before = mask;
				const auto active = static_cast<std::size_t>(
				    std::count_if(mask.begin(), mask.end(), [](auto value) { return value != 0; }));
				reject.rejectLargestOnly =
				    repeatedModifications >= 4 || active <= static_cast<std::size_t>(6 * ((parameterCount + 1) / 2));
				const auto decision =
				    Calculations::applyHystereticReject(observations, mask, chi, reject, selection.locked);
				if (decision.status == Calculations::RejectStatus::NoChange) {
					status = FitSelectionStatus::Stable;
					LOG_INFO << "[FitSession] selection status=NO_CHANGE";
					return false;
				}
				if (decision.status == Calculations::RejectStatus::SafeguardRollback) {
					status = FitSelectionStatus::SafeguardRollback;
					return false;
				}
				if (changes >= maxChanges) {
					mask = before;
					status = FitSelectionStatus::PassLimit;
					return false;
				}
				int dimensions = 0;
				for (std::size_t i = 0; i < observations.size(); ++i)
					if (mask[i])
						dimensions += observations[i].measurementDimensions();
				if (static_cast<std::size_t>(decision.nKept) < std::max(std::size_t{3}, clipCfg.minObsKept) ||
				    dimensions <= parameterCount) {
					mask = before;
					status = FitSelectionStatus::InsufficientObservations;
					return false;
				}
				const int modifications = decision.nFlippedIn + decision.nFlippedOut;
				repeatedModifications = modifications == previousModifications ? repeatedModifications + 1 : 0;
				previousModifications = modifications;
				++changes;
				status = FitSelectionStatus::NotApplied;
				fit::reportProgress(control, "σ-clip update " + std::to_string(changes) +
				                                 ": chi2=" + std::to_string(decision.appliedRejectChiSquared) +
				                                 ", in=" + std::to_string(decision.nFlippedIn) +
				                                 " out=" + std::to_string(decision.nFlippedOut));
				return true;
			};
			FitResult result = passes.run({.initial = initial,
			                               .observations = kept.observations,
			                               .prepared = prepared,
			                               .warmStartNonGrav = warmStartNg,
			                               .uncertainty = options.uncertainty,
			                               .included = &included,
			                               .updateSelection = update});
			if (!result.converged)
				status = FitSelectionStatus::CorrectionFailed;
			kept = buildKept();
			result.selectionStatus = status;
			result.converged =
			    result.converged && (status == FitSelectionStatus::Stable || status == FitSelectionStatus::NotApplied);
			result.sigmaClipPassesRun = changes;
			result.sigmaClipRejected = static_cast<int>(observations.size() - kept.observations.size());
			result.keptObsIndices = std::move(kept.indices);
			return result;
		}

		FitResult result = passes.run({.initial = initial,
		                               .observations = kept.observations,
		                               .prepared = prepared,
		                               .warmStartNonGrav = warmStartNg,
		                               .uncertainty = options.uncertainty,
		                               .included = &included});
		result.keptObsIndices = std::move(kept.indices);
		return result;
	} catch (const fit::FitCancelled&) {
		return {.cancelled = true};
	}

	FitResult FitSession::optimize(const FitRequest& request, const FitRunControl& control) const
	{
		return FitRun{config_, control, request.outputEpochJdTdb.value_or(0.0)}.run(request);
	}

	FitResult fit::runPreparedFitSession(const FitConfig& config, const FitRequest& request,
	                                     std::span<const obs::PreparedObservation> prepared,
	                                     const FitRunControl& control)
	{
		return FitRun{config, control, request.outputEpochJdTdb.value_or(0.0), prepared}.run(request);
	}

}  // namespace orbital_core
