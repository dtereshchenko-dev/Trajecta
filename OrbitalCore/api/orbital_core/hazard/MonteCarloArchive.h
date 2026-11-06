#pragma once

#include "orbital_core/Exporter.h"
#include "orbital_core/OrbitSolution.h"
#include "orbital_core/hazard/ObsMonteCarlo.h"

namespace orbital_core::hazard {

	enum class CovarianceSampling { State, Parameters };

	struct MonteCarloArchiveInput {
		OrbitSolution solution;
		ObsMcConfig sampling;
		CovarianceSampling kind{CovarianceSampling::State};
		// Provenance only: solution.covariance is already scaled. Never apply twice.
		double covarianceObjectiveScale{1.0};
		std::optional<double> albedo;
	};

	// Synchronous export into a new work/<object>/MC/run-* directory. Publishes the
	// manifest last; never starts atmospheric propagation or overwrites older runs.
	[[nodiscard]] ORBITAL_CORE_API ExportReport saveMonteCarloArchive(const std::filesystem::path& solutionDirectory,
	                                                                  const MonteCarloArchiveInput& input,
	                                                                  const NBodyPropagator& propagator,
	                                                                  std::span<const CloseApproachEvent> events,
	                                                                  const ObsMcResult& result);

}  // namespace orbital_core::hazard
