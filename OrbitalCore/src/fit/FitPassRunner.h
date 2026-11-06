#pragma once

#include "fit/OrbitFitConfig.h"
#include "obs/PreparedObservation.h"
#include "orbital_core/FitSession.h"

#include <span>
#include <vector>
#include <functional>

namespace orbital_core::fit {

	// Borrows the complete prepared pool and the current selection for one fit.
	struct FitPass {
		Orbit_p initial;
		std::span<const ObsNode> observations;
		std::span<const obs::PreparedObservation> prepared;
		const NonGravModel* warmStartNonGrav{nullptr};
		FitUncertaintyMode uncertainty{FitUncertaintyMode::Compute};
		bool evaluateOnly{false};
		std::vector<std::uint8_t>* included{nullptr};
		std::function<bool(std::vector<std::uint8_t>&, const std::vector<double>&, int)> updateSelection;
	};

	// Internal synchronous executor; config and run control outlive the runner.
	class FitPassRunner {
	   public:
		FitPassRunner(const FitConfig& config, double outputEpochJdTdb, const FitRunControl& control);

		[[nodiscard]] std::vector<obs::PreparedObservation> prepareObservations(
		    std::span<const ObsNode> observations) const;
		[[nodiscard]] FitResult run(const FitPass& pass);

	   private:
		[[nodiscard]] OrbitOptimizerConfig optimizerConfig() const noexcept;

		const FitConfig& config_;
		double outputEpochJdTdb_;
		const FitRunControl& control_;
	};

}  // namespace orbital_core::fit
