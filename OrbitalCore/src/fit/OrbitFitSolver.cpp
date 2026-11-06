#include "fit/OrbitFitSolver.h"
#include "fit/FitCancellation.h"
#include "orbital_core/optimizers/EigenLevenbergMarquardt.h"
#include "logger/Logger.h"
#include <iomanip>

namespace orbital_core::fit {
	namespace {
		struct GravityBootstrap {
			NonGravModel& model;
			NonGravModel saved;

			explicit GravityBootstrap(NonGravModel& target) : model(target), saved(target) { model.setEnabled(false); }

			~GravityBootstrap() { model = saved; }
		};
	}  // namespace

	std::vector<double> OrbitFitSolver::solve(const FitRunControl& control)
	{
		checkCancellation(control.cancel);
		summary_ = {};
		auto initial = OrbitFitModel::stateVectorFromOrbit(model_.initialOrbit());
		if (config_.evaluateOnly) {
			model_.differenceSteps() = {1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8};
			model_.nonGravProfile().extendParams(initial, model_.differenceSteps());
			reportProgress(control, "Evaluating supplied orbit...");
			return runCorrection(initial, control.cancel);
		}
		return model_.propagator().nonGrav().enabled() ? solveNonGrav(initial, control)
		                                               : solveGravity(initial, control);
	}

	void OrbitFitSolver::logState(const char* tag, const std::vector<double>& p, double cost) const
	{
		auto orb = model_.orbitFromStateVector(p);
		LOG_INFO << "[" << tag << "] cost=" << std::scientific << std::setprecision(6) << cost << "  a=" << std::fixed
		         << std::setprecision(9) << orb->getSemiMajorAxis() << "  e=" << orb->getExcentricity()
		         << "  i=" << Convertors::radiansToDegrees(orb->getInclination())
		         << "  O=" << Convertors::radiansToDegrees(orb->getLongitudeOfAscendingNode())
		         << "  w=" << Convertors::radiansToDegrees(orb->getArgumentOfPericentre())
		         << "  M=" << Convertors::radiansToDegrees(orb->getMeanAnomaly());
		if (p.size() == 7) {
			LOG_INFO << "[" << tag << "] A2=" << std::scientific << std::setprecision(6) << p[6] << " (Yarkovsky)";
		}

		else if (p.size() >= 9) {
			LOG_INFO << "[" << tag << "] A1=" << std::scientific << std::setprecision(6) << p[6] << "  A2=" << p[7]
			         << "  A3=" << p[8] << ((p.size() >= 10) ? ("  DT=" + std::to_string(p[9])) : std::string{});
		}
	}

	std::vector<double> OrbitFitSolver::bootstrapNonGrav(const std::vector<double>& initialParams,
	                                                     const FitRunControl& control)
	{
		const bool wantsNBody = config_.nbody.enabled;
		auto& nonGrav = model_.propagator().nonGrav();
		std::vector<double> seedParams = initialParams;
		if (!config_.skipNonGravBootstrap) {
			const GravityBootstrap bootstrap(nonGrav);

			reportProgress(control,
			               wantsNBody ? "N-body gravity-only bootstrap..." : "Numerical gravity-only bootstrap...");

			if (config_.logDiagnostics) {
				LOG_INFO << "=== Non-grav bootstrap: " << (wantsNBody ? "N-body" : "Sun-only numerical")
				         << " 6-param LM before A solve ===";
			}

			seedParams = runCorrection(initialParams, control.cancel);
			const double bootstrapCost = summary_.cost;
			if (config_.logDiagnostics)
				logState("after gravity bootstrap", seedParams, model_.computeCost(seedParams));
			if (config_.logDiagnostics)
				LOG_INFO << "Non-grav bootstrap done: cost=" << bootstrapCost;
		}
		return seedParams;
	}

	std::vector<double> OrbitFitSolver::solveExtended(const std::vector<double>& seedParams, std::stop_token cancel)
	{
		const bool wantsNBody = config_.nbody.enabled;
		auto& nonGrav = model_.propagator().nonGrav();
		const auto& profile = model_.nonGravProfile();
		std::vector<double> finalParams = seedParams;
		double startCost = std::numeric_limits<double>::quiet_NaN();
		int pass = 0;
		bool dtUnlocked = false;

		for (;;) {
			checkCancellation(cancel);
			std::vector<double> params(finalParams.begin(), finalParams.begin() + 6);
			model_.differenceSteps() = {1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8};
			profile.extendParams(params, model_.differenceSteps());

			const bool suppressDt = hasNgParameter(nonGrav.fitParameterMask(), NgParameterMask::DT) &&
			                        params.size() == 10 && nonGrav.A1() == 0.0 && nonGrav.A2() == 0.0 &&
			                        nonGrav.A3() == 0.0;
			if (suppressDt) {
				params.pop_back();
				model_.differenceSteps().pop_back();
			}

			if (pass == 0)
				startCost = model_.computeCost(params);
			checkCancellation(cancel);

			if (config_.logDiagnostics) {
				LOG_INFO << "=== Non-grav extended solve pass " << (pass + 1) << " (" << params.size() << "-param, "
				         << (wantsNBody ? "N-body" : "Sun-only numerical")
				         << " force model, integrator=" << model_.propagator().integratorName()
				         << ", initialStep=" << model_.propagator().stepSize()
				         << "d, fixedAccuracy=" << model_.propagator().fixedAccuracyParameter() << ") ===";
				logState(pass == 0 ? "non-grav start" : "non-grav restart", params,
				         pass == 0 ? startCost : model_.computeCost(params));
			}

			finalParams = runCorrection(params, cancel);
			profile.updateFromParams(finalParams, 6);

			const bool dtCanEnter = suppressDt && !dtUnlocked &&
			                        hasNgParameter(nonGrav.fitParameterMask(), NgParameterMask::DT) &&
			                        (nonGrav.A1() != 0.0 || nonGrav.A2() != 0.0 || nonGrav.A3() != 0.0);
			if (!dtCanEnter)
				break;

			dtUnlocked = true;
			++pass;
		}
		if (config_.logDiagnostics)
			logState("after non-grav", finalParams, model_.computeCost(finalParams));
		logNonGravResult(finalParams, startCost);
		return finalParams;
	}

	void OrbitFitSolver::logNonGravResult(const std::vector<double>& finalParams, double gravCost) const
	{
		const auto& ng = model_.propagator().nonGrav();
		const int n = ng.numFitParams();
		if (ng.model() == NgModel::Yarkovsky && n >= 1) {
			LOG_INFO << "Non-grav done. A2=" << std::scientific << std::setprecision(6) << finalParams[6]
			         << " (Yarkovsky, A1=A3=0 by construction)"
			         << "  Δcost = " << std::scientific << std::setprecision(3) << gravCost << " → " << summary_.cost;
		}

		else if (n == 2) {
			LOG_INFO << "Non-grav done. A1=" << std::scientific << std::setprecision(6) << finalParams[6]
			         << " A2=" << finalParams[7] << "  deltaCost = " << std::scientific << std::setprecision(3)
			         << gravCost << " -> " << summary_.cost;
		}

		else if (n >= 3) {
			LOG_INFO << "Non-grav done. A1=" << std::scientific << std::setprecision(6) << finalParams[6]
			         << " A2=" << finalParams[7] << " A3=" << finalParams[8]
			         << ((n >= 4 && finalParams.size() >= 10) ? (" DT=" + std::to_string(finalParams[9]))
			                                                  : std::string{})
			         << "  Δcost = " << std::scientific << std::setprecision(3) << gravCost << " → " << summary_.cost;
		}

		else {
			LOG_INFO << "Non-grav done. (no fit params)";
		}
	}

	std::vector<double> OrbitFitSolver::solveNonGrav(const std::vector<double>& initialParams,
	                                                 const FitRunControl& control)
	{
		const bool wantsNBody = config_.nbody.enabled;
		if (!wantsNBody && config_.logDiagnostics)
			LOG_INFO << "Non-grav enabled: using numerical Sun + non-grav propagation "
			         << "with planetary perturbers disabled.";
		model_.differenceSteps() = {1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8};
		const auto seed = bootstrapNonGrav(initialParams, control);
		reportProgress(control, wantsNBody ? "N-body + non-grav solve..." : "Numerical non-grav solve...");
		return solveExtended(seed, control.cancel);
	}

	std::vector<double> OrbitFitSolver::solveGravity(const std::vector<double>& initialParams,
	                                                 const FitRunControl& control)
	{
		bool wantsNBody = model_.usesNBody();
		std::vector<double> params;

		if (wantsNBody) {
			reportProgress(control, "N-body LM...");

			if (config_.logDiagnostics) {
				LOG_INFO << "=== N-body LM (integrator=" << model_.propagator().integratorName()
				         << ", initialStep=" << model_.propagator().stepSize()
				         << "d, fixedAccuracy=" << model_.propagator().fixedAccuracyParameter() << ") ===";
			}
			params = runCorrection(initialParams, control.cancel);
			const double nbodyCost = summary_.cost;
			if (config_.logDiagnostics)
				LOG_INFO << "N-body LM done: cost=" << nbodyCost;
		}

		else {
			reportProgress(control, "Kepler LM...");

			if (config_.logDiagnostics)
				LOG_INFO << "=== Kepler LM ===";
			params = runCorrection(initialParams, control.cancel);
			if (config_.logDiagnostics)
				LOG_INFO << "Kepler LM done: cost=" << summary_.cost;
		}

		return params;
	}

	std::vector<double> OrbitFitSolver::refine(const std::vector<double>& params, const FitRunControl& control,
	                                           std::function<bool(const std::vector<double>&)> updateSelection)
	{
		summary_ = {};
		return runCorrection(params, control.cancel, std::move(updateSelection));
	}

	std::vector<double> OrbitFitSolver::runCorrection(const std::vector<double>& initialParams, std::stop_token cancel,
	                                                  std::function<bool(const std::vector<double>&)> updateSelection)
	{
		using namespace orbital_core::optimizers;
		checkCancellation(cancel);
		problem_.beginPass(static_cast<int>(initialParams.size()));
		if (config_.evaluateOnly) {
			problem_.saveLinearizationFor(initialParams);
			checkCancellation(cancel);
			summary_.cost = problem_.displayCost(initialParams);
			summary_.converged = false;
			summary_.cancelled = false;
			return initialParams;
		}
		if (problem_.activeResidualCount() < static_cast<int>(initialParams.size())) {
			summary_.cost = std::numeric_limits<double>::quiet_NaN();
			summary_.converged = false;
			return initialParams;
		}
		const auto result = EigenLevenbergMarquardt().optimize(
		    initialParams, problem_, [cancel] { return cancel.stop_requested(); }, std::move(updateSelection));
		checkCancellation(cancel);

		if (!result.cancelled) {
			problem_.saveLinearizationFor(result.params);
			checkCancellation(cancel);
		}

		summary_.iterations += result.iterations;
		summary_.cost = result.finalCost;
		summary_.converged = result.converged;
		summary_.cancelled = result.cancelled;

		return result.params;
	}

}  // namespace orbital_core::fit
