#pragma once

#include "fit/OrbitFitUncertainty.h"
#include "fit/OrbitFitSolver.h"
#include <span>

namespace orbital_core::fit {
	struct PhotometricOrbitSigmas {
		std::optional<double> asteroid;
		std::optional<double> cometTotal;
		std::optional<double> cometNuclear;
	};

	// Builds an owning result; diagnostics reuse the accepted linearization.
	class OrbitFitOutput {
	   public:
		OrbitFitOutput(const OrbitFitModel& model, const OrbitFitProblem& problem, const OrbitOptimizerConfig& config)
		    : model_(model), problem_(problem), config_(config), uncertainty_(model, problem, config.outputEpochTdb)
		{
		}

		FitResult build(const std::vector<double>& params, const FitSolveSummary& summary,
		                std::stop_token cancel) const;
		void finalize(FitResult& result, const std::vector<double>& params, const std::vector<double>& rawResiduals,
		              std::stop_token cancel) const;
		// External formal covariance uses physical parameters at the fit epoch.
		PostFitDiagnostics diagnostics(const FitResult& result, const FitParameterCovariance& formalCovariance,
		                               std::span<const std::uint8_t> included) const;
		[[nodiscard]] std::vector<double> selectionChiSquares(const std::vector<double>& params) const;
		// Covariance uses the optimizer's canonical parameter units.
		[[nodiscard]] PhotometricOrbitSigmas photometricOrbitSigmas(const Residuals& nominal,
		                                                            const std::vector<double>& params,
		                                                            const Eigen::MatrixXd& covariance,
		                                                            std::stop_token cancel) const;

	   private:
		// Match Residuals::rms units and summation order without constructing report rows.
		double reportRms(std::span<const double> rawResiduals) const;
		void populateUncertainty(FitResult& result, const std::vector<double>& params,
		                         const std::vector<double>& rawResiduals, std::stop_token cancel) const;
		Residuals populateResidualsFromRaw(const std::vector<double>& rawResiduals,
		                                   const std::vector<double>& params) const;
		Orbit_p propagateToTargetEpoch(const std::vector<double>& params, const Orbit& fitOrbit) const;

		const OrbitFitModel& model_;
		const OrbitFitProblem& problem_;
		const OrbitOptimizerConfig& config_;
		OrbitFitUncertainty uncertainty_;
	};

}  // namespace orbital_core::fit
