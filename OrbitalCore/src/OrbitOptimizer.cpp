#include "OrbitOptimizer.h"
#include "fit/FitCancellation.h"
#include "fit/OrbitFitModel.h"
#include "fit/OrbitFitOutput.h"
#include "fit/OrbitFitSolver.h"
#include "orbital_core/DE440Source.h"
#include "logger/Logger.h"
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace {
	struct CacheStatsLogger {
		bool enabled;

		explicit CacheStatsLogger(bool log) : enabled(log)
		{
			if (enabled)
				DE440Source::resetCacheStats();
		}

		~CacheStatsLogger()
		{
			if (!enabled)
				return;
			const auto stats = DE440Source::getCacheStats();
			LOG_INFO << "[DE440 cache] hits=" << stats.hits << " misses=" << stats.misses
			         << " evictions=" << stats.evictions << " hitRate=" << std::fixed << std::setprecision(2)
			         << stats.hitRate() * 100.0 << "%"
			         << " total=" << stats.totalLookups();
		}
	};

	std::vector<double> canonicalParameters(std::span<const double> physical)
	{
		std::vector<double> params(physical.begin(), physical.end());
		if (params.size() >= 6) {
			params[3] /= Constants::GAUSS_K;
			params[4] /= Constants::GAUSS_K;
			params[5] /= Constants::GAUSS_K;
		}
		return params;
	}
}  // namespace

struct OrbitOptimizer::Impl {
	OrbitOptimizerConfig config;
	orbital_core::fit::OrbitFitModel model;
	orbital_core::fit::OrbitFitProblem problem;
	orbital_core::fit::OrbitFitSolver solver;
	orbital_core::fit::OrbitFitOutput output;
	orbital_core::FitResult result;

	Impl(Orbit_p initial, std::vector<orbital_core::obs::PreparedObservation> observations,
	     OrbitOptimizerConfig settings)
	    : config(std::move(settings)),
	      model(std::move(initial), std::move(observations), config),
	      problem(model),
	      solver(model, problem, config),
	      output(model, problem, config)
	{
	}
};

OrbitOptimizer::OrbitOptimizer(Orbit_p initial, std::vector<orbital_core::obs::PreparedObservation> observations,
                               OrbitOptimizerConfig config)
    : d(std::make_unique<Impl>(std::move(initial), std::move(observations), std::move(config)))
{
}

OrbitOptimizer::~OrbitOptimizer() = default;

const orbital_core::FitResult& OrbitOptimizer::optimize(const orbital_core::FitRunControl& control)
{
	const auto cancel = control.cancel;
	try {
		orbital_core::fit::checkCancellation(cancel);
		d->problem.setIncluded({});
		const CacheStatsLogger cacheStats(d->config.logCacheStats);
		const auto params = d->solver.solve(control);
		orbital_core::fit::checkCancellation(cancel);
		d->result = d->output.build(params, d->solver.summary(), cancel);
		orbital_core::fit::checkCancellation(cancel);
	} catch (const orbital_core::fit::FitCancelled&) {
		d->result = {.cancelled = true};
	}
	return d->result;
}

const orbital_core::FitResult& OrbitOptimizer::optimizeWithSelection(std::vector<std::uint8_t>& included,
                                                                     const SelectionUpdate& update,
                                                                     const orbital_core::FitRunControl& control)
{
	try {
		orbital_core::fit::checkCancellation(control.cancel);
		const CacheStatsLogger cacheStats(d->config.logCacheStats);
		d->problem.setIncluded(included);
		auto params = d->solver.solve(control);
		orbital_core::fit::checkCancellation(control.cancel);
		auto summary = d->solver.summary();
		const auto select = [&](const std::vector<double>& current) {
			orbital_core::fit::checkCancellation(control.cancel);
			d->problem.saveLinearizationFor(current);
			const auto chi = d->output.selectionChiSquares(current);
			const bool changed = update(included, chi, static_cast<int>(current.size()));
			orbital_core::fit::checkCancellation(control.cancel);
			if (!changed)
				return false;
			d->problem.setIncluded(included);
			return true;
		};
		if (summary.converged && update && select(params)) {
			params = d->solver.refine(params, control, select);
			const auto refinement = d->solver.summary();
			summary.iterations += refinement.iterations;
			summary.cost = refinement.cost;
			summary.converged = refinement.converged;
			summary.cancelled = refinement.cancelled;
		}
		orbital_core::fit::checkCancellation(control.cancel);
		d->model.nonGravProfile().updateFromParams(params, 6);
		d->result = d->output.build(params, summary, control.cancel);
		d->result.allResiduals = d->result.residuals;
		const auto diagnostics = d->output.diagnostics(d->result, {}, included);
		if (diagnostics.valid) {
			d->result.allResiduals = diagnostics.residuals;
			d->result.uncertainty.fitParameterCovarianceFormal = diagnostics.formalCovariance;
		}
		Residuals selected;
		const auto& all = diagnostics.valid ? diagnostics.residuals : d->result.residuals;
		if (all.count() == included.size()) {
			for (std::size_t i = 0; i < included.size(); ++i)
				if (included[i])
					selected.add(all.data()[i]);
			d->result.residuals = std::move(selected);
			d->result.cost = d->result.residuals.rms();
		}
	} catch (const orbital_core::fit::FitCancelled&) {
		d->result = {.cancelled = true};
	}
	return d->result;
}

const orbital_core::FitResult& OrbitOptimizer::finalizeOutput(double outputEpochJdTdb, bool computeUncertainty,
                                                              const orbital_core::FitRunControl& control)
{
	try {
		orbital_core::fit::checkCancellation(control.cancel);
		const auto& snapshot = d->problem.snapshot();
		if (d->result.cancelled || !d->result.fitEpochOrbit || !snapshot.valid)
			throw std::logic_error("Cannot finalize without an accepted fit linearization");
		auto config = d->config;
		config.outputEpochTdb = outputEpochJdTdb;
		config.computeUncertainty = computeUncertainty;
		const orbital_core::fit::OrbitFitOutput output(d->model, d->problem, config);
		auto result = d->result;
		result.fitEpochOrbit = d->model.orbitFromStateVector(snapshot.params);
		result.uncertainty = {};
		output.finalize(result, snapshot.params, snapshot.rawResiduals, control.cancel);
		d->result = std::move(result);
	} catch (const orbital_core::fit::FitCancelled&) {
		d->result = {.cancelled = true};
	}
	return d->result;
}

double OrbitOptimizer::weightedObservationChiSquarePhysical(std::span<const double> physicalParams) const
{
	return d->problem.weightedChiSquare(canonicalParameters(physicalParams));
}

OrbitOptimizer::WeightedObservationLinearization OrbitOptimizer::weightedObservationLinearizationPhysical(
    std::span<const double> physicalParams) const
{
	if (physicalParams.size() < 6)
		return {};
	return d->problem.weightedLinearization(canonicalParameters(physicalParams));
}

OrbitOptimizer::PostFitDiagnostics OrbitOptimizer::postFitDiagnostics(
    const orbital_core::FitParameterCovariance& formalCovariance, std::span<const std::uint8_t> included) const
{
	if (d->result.cancelled)
		return {};
	return d->output.diagnostics(d->result, formalCovariance, included);
}
