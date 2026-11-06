#include "fit/FitPassRunner.h"

#include "OrbitOptimizer.h"
#include "fit/FitCancellation.h"
#include "fit/FitEpochPolicy.h"
#include "fit/FitObservations.h"
#include "fit/OptimizerConfigFactory.h"
#include "logger/Logger.h"
#include <stdexcept>

namespace orbital_core::fit {

	FitPassRunner::FitPassRunner(const FitConfig& config, double outputEpochJdTdb, const FitRunControl& control)
	    : config_(config), outputEpochJdTdb_(outputEpochJdTdb), control_(control)
	{
	}

	OrbitOptimizerConfig FitPassRunner::optimizerConfig() const noexcept
	{
		auto config = makeOrbitOptimizerConfig(config_);
		config.outputEpochTdb = outputEpochJdTdb_;
		return config;
	}

	std::vector<obs::PreparedObservation> FitPassRunner::prepareObservations(
	    std::span<const ObsNode> observations) const
	{
		return prepareFitObservations(observations, optimizerConfig());
	}

	FitResult FitPassRunner::run(const FitPass& pass)
	{
		const auto cancel = control_.cancel;
		checkCancellation(cancel);
		if (!pass.included)
			throw std::invalid_argument("Fit pass requires an observation mask");
		const auto& initial = pass.initial;
		const auto& obs = pass.observations;
		const auto* warmStartNg = pass.warmStartNonGrav;
		OrbitOptimizerConfig oc = optimizerConfig();
		oc.computeUncertainty = pass.uncertainty == FitUncertaintyMode::Compute;
		oc.evaluateOnly = pass.evaluateOnly;

		if (!oc.nonGrav.enabled || !warmStartNg || !warmStartNg->enabled() ||
		    warmStartNg->model() != NonGravModel::canonicalModel(static_cast<NgModel>(oc.nonGrav.model)))
			warmStartNg = nullptr;
		if (warmStartNg) {
			oc.nonGrav.useInitialValues = true;
			oc.nonGrav.initialA1 = warmStartNg->A1();
			oc.nonGrav.initialA2 = warmStartNg->A2();
			oc.nonGrav.initialA3 = warmStartNg->A3();
			oc.nonGrav.initialDT = warmStartNg->DT();
			oc.skipNonGravBootstrap = true;
		}

		auto propagationConfig = config_;
		propagationConfig.nonGrav = oc.nonGrav;
		Orbit_p fitInitial = pass.evaluateOnly
		                         ? std::make_shared<Orbit>(*initial)
		                         : rebaseInitialOrbitToFitEpoch(initial, obs, propagationConfig, nullptr, false);
		checkCancellation(cancel);

		OrbitOptimizer optimizer(fitInitial,
		                         std::vector<obs::PreparedObservation>(pass.prepared.begin(), pass.prepared.end()), oc);
		FitResult result = optimizer.optimizeWithSelection(*pass.included, pass.updateSelection, control_);
		checkCancellation(cancel);
		return result;
	}

}  // namespace orbital_core::fit
