#pragma once

#include "fit/OrbitFitProblem.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/FitRunControl.h"
#include <stop_token>
#include <string>
#include <functional>

namespace orbital_core::fit {

	struct FitSolveSummary {
		int iterations{0};
		double cost{0.0};
		bool converged{false};
		bool cancelled{false};
	};

	// Gravity bootstrap, optional DT activation, and least-squares correction.
	// It does not prepare observations, build covariance, or change fit epochs.
	class OrbitFitSolver {
	   public:
		OrbitFitSolver(OrbitFitModel& model, OrbitFitProblem& problem, const OrbitOptimizerConfig& config)
		    : model_(model), problem_(problem), config_(config)
		{
		}

		std::vector<double> solve(const FitRunControl& control);
		std::vector<double> refine(const std::vector<double>& params, const FitRunControl& control,
		                           std::function<bool(const std::vector<double>&)> updateSelection);

		[[nodiscard]] const FitSolveSummary& summary() const noexcept { return summary_; }

	   private:
		std::vector<double> runCorrection(const std::vector<double>& params, std::stop_token cancel,
		                                  std::function<bool(const std::vector<double>&)> updateSelection = {});
		std::vector<double> solveGravity(const std::vector<double>& initial, const FitRunControl& control);
		std::vector<double> solveNonGrav(const std::vector<double>& initial, const FitRunControl& control);
		std::vector<double> bootstrapNonGrav(const std::vector<double>& initial, const FitRunControl& control);
		std::vector<double> solveExtended(const std::vector<double>& seed, std::stop_token cancel);
		void logState(const char* tag, const std::vector<double>& params, double cost) const;
		void logNonGravResult(const std::vector<double>& params, double startCost) const;

		OrbitFitModel& model_;
		OrbitFitProblem& problem_;
		const OrbitOptimizerConfig& config_;
		FitSolveSummary summary_;
	};

}  // namespace orbital_core::fit
