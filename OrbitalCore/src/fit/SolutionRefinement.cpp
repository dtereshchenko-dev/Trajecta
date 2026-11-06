#include "orbital_core/SolutionRefinement.h"
#include "fit/FitCancellation.h"
#include "fit/ObservationSelection.h"
#include "orbital_core/Types.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core {
	namespace {
		bool changedInitialValues(const NonGravSettings& a, const NonGravSettings& b)
		{
			return a.useInitialValues &&
			       (!b.useInitialValues || a.initialA1 != b.initialA1 || a.initialA2 != b.initialA2 ||
			        a.initialA3 != b.initialA3 || a.initialDT != b.initialDT);
		}
	}  // namespace

	FitPipelineResult refineSolution(const SolutionRefinementRequest& request, const FitRunControl& control)
	try {
		fit::checkCancellation(control.cancel);
		const auto& previous = request.previous;
		if (previous.fit.cancelled || !previous.fit.orbit || previous.observations.all.empty())
			throw std::invalid_argument("Refinement requires an existing solution and its observations");
		if (request.preparation.weightProfile != request.config.weightProfile)
			throw std::invalid_argument("Preparation and fit weighting profiles must match");

		FitPipelineResult result;
		result.observations = previous.observations;
		const auto& prepared = result.observations.preparation;
		if (!prepared || prepared->debiasRequested != request.preparation.debias ||
		    prepared->weightProfile != request.preparation.weightProfile) {
			result.observations = prepareObservationCopy(previous.observations, request.preparation);
		}
		fit::checkCancellation(control.cancel);
		result.effectiveConfig = request.config;
		result.initial = previous.initial;
		result.selection = previous.selection;
		result.observationOverrides.assign(request.observationOverrides.begin(), request.observationOverrides.end());

		std::vector<std::uint8_t> initial(result.observations.all.size(), previous.fit.keptObsIndices.empty() ? 1 : 0);
		for (const auto i : previous.fit.keptObsIndices)
			initial.at(i) = 1;
		const auto selection =
		    fit::resolveObservationSelection(result.observations.all, initial, result.observationOverrides);
		result.fitInputObservations = std::ranges::count(selection.included, std::uint8_t{1});
		NonGravModel target;
		target.setEnabled(request.config.nonGrav.enabled);
		target.setModel(static_cast<NgModel>(request.config.nonGrav.model));
		target.setFitParameterMask(
		    request.config.nonGrav.solveMask
		        ? static_cast<NgParameterMask>(request.config.nonGrav.solveMask)
		        : NonGravModel::defaultFitParameterMaskFor(static_cast<NgModel>(request.config.nonGrav.model)));
		const auto parameters = 6 + (target.enabled() ? target.numFitParams() : 0);
		int dimensions = 0;
		for (std::size_t i = 0; i < selection.included.size(); ++i)
			if (selection.included[i])
				dimensions += result.observations.all[i].measurementDimensions();
		if (dimensions <= parameters)
			throw std::invalid_argument("Too few selected observations for the requested fit parameters");

		const auto* warmStart = changedInitialValues(request.config.nonGrav, previous.effectiveConfig.nonGrav)
		                            ? nullptr
		                            : &previous.fit.finalNonGrav;
		fit::reportProgress(control, "Refining current solution...");
		result.fit = FitSession{result.effectiveConfig}.optimize(
		    {.initialOrbit = previous.fit.fitEpochOrbit ? previous.fit.fitEpochOrbit : previous.fit.orbit,
		     .observations = result.observations.all,
		     .initialIncluded = initial,
		     .warmStartNonGrav = warmStart,
		     .options = {.selection = result.selection},
		     .outputEpochJdTdb = previous.fit.orbit->getEpoch(),
		     .observationOverrides = result.observationOverrides},
		    control);
		if (result.fit.cancelled)
			return {.fit = {.cancelled = true}};
		fit::checkCancellation(control.cancel);
		if (!result.fit.orbit)
			throw std::runtime_error("Refinement returned no orbit");
		return result;
	} catch (const fit::FitCancelled&) {
		return {.fit = {.cancelled = true}};
	}
}  // namespace orbital_core
