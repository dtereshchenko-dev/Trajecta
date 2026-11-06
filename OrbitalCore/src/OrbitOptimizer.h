#pragma once

#include "fit/OrbitFitConfig.h"
#include "fit/OrbitFitLinearization.h"
#include "obs/PreparedObservation.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/FitRunControl.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <memory>
#include <span>
#include <functional>

// Single-sample optimizer. Selection/clipping and epoch policy live in FitSession.
// It consumes fit-ready observations (fit::prepareFitObservations) and never
// prepares them, so a caller fitting several subsets prepares its set once.
class ORBITAL_CORE_API OrbitOptimizer {
   public:
	using WeightedObservationLinearization = orbital_core::fit::WeightedObservationLinearization;

	using PostFitDiagnostics = orbital_core::fit::PostFitDiagnostics;
	using SelectionUpdate = std::function<bool(std::vector<std::uint8_t>&, const std::vector<double>&, int)>;

	OrbitOptimizer(Orbit_p initialOrbit, std::vector<orbital_core::obs::PreparedObservation> observations,
	               OrbitOptimizerConfig config = {});
	~OrbitOptimizer();
	OrbitOptimizer(const OrbitOptimizer&) = delete;
	OrbitOptimizer& operator=(const OrbitOptimizer&) = delete;

	// Borrowed until optimize()/finalizeOutput() or destruction. Copy for ownership.
	// FitResult keeps display-epoch elements separate from fit-epoch state/COV.
	// A cancelled call returns no solution and does not publish partial covariance.
	[[nodiscard]] const orbital_core::FitResult& optimize(const orbital_core::FitRunControl& control = {});
	[[nodiscard]] const orbital_core::FitResult& optimizeWithSelection(std::vector<std::uint8_t>& included,
	                                                                   const SelectionUpdate& update,
	                                                                   const orbital_core::FitRunControl& control = {});

	// Finalizes the accepted linearization without fitting or rebasing it again.
	// Invalidates previously borrowed results, just like optimize().
	[[nodiscard]] const orbital_core::FitResult& finalizeOutput(double outputEpochJdTdb, bool computeUncertainty,
	                                                            const orbital_core::FitRunControl& control = {});

	[[nodiscard]] double weightedObservationChiSquarePhysical(std::span<const double> physicalParams) const;
	[[nodiscard]] WeightedObservationLinearization weightedObservationLinearizationPhysical(
	    std::span<const double> physicalParams) const;
	// Empty covariance uses the accepted fit; external COV is physical at fit epoch.
	[[nodiscard]] PostFitDiagnostics postFitDiagnostics(
	    const orbital_core::FitParameterCovariance& formalCovariance = {},
	    std::span<const std::uint8_t> included = {}) const;

   private:
	struct Impl;
	std::unique_ptr<Impl> d;
};
