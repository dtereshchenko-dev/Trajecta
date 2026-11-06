#pragma once

#include "PreparedObservation.h"
#include "orbital_core/AstrometricWeights.h"

#include <span>
#include <vector>

class ObsNode;

namespace orbital_core::obs {
	// ---------------------------------------------------------------------
	// ObservationPipeline — internal preparation layer between raw MPC
	// observations (`ObsNode`) and the LM hot path (`PreparedObservation`).
	//
	// Lifecycle:
	//   1. Construction warms up the astrometric σ table.  Stations resolve
	//      through the process-wide observatoryCatalog(), so constructing a
	//      pipeline does not read `OBSCODE.dat`.  Any fatal I/O error throws
	//      `std::runtime_error`.
	//   2. `prepare()` runs the per-obs pipeline:
	//        a. Resolve observatory data (longitude / parallax constants).
	//        b. UTC → TDB normalisation (DE440 is indexed by TDB).
	//        c. Earth heliocentric state at jdTdb.
	//        d. Carry through satellite offset for spacecraft observations.
	//        e. Per-observation sigmas from the selected weighting profile
	//           and, when enabled, the observer's reported ADES uncertainty.
	//        f. Sort by jdTdb so light-time iteration walks the arc
	//           monotonically.
	//
	// The class holds no observation state between calls.  Every step is
	// per observation; the set-dependent vfcc17 batch relaxation is applied
	// by fit::prepareFitObservations.  Observation selection belongs to
	// FitSession.
	// ---------------------------------------------------------------------
	class ObservationPipeline {
	   public:
		// `options` selects which uncertainty each observation is weighted by
		// (observer-supplied vs table, and whether the timing term applies);
		// see ObservationWeightingOptions.  Both knobs default off.
		explicit ObservationPipeline(AstrometricWeightProfile weightProfile = AstrometricWeightProfile::Vfcc17,
		                             ObservationWeightingOptions options = {});

		ObservationPipeline(const ObservationPipeline&) = delete;
		ObservationPipeline& operator=(const ObservationPipeline&) = delete;
		ObservationPipeline(ObservationPipeline&&) = default;
		ObservationPipeline& operator=(ObservationPipeline&&) = default;

		// Convert raw observations into the LM-ready prepared form.
		// Output is sorted by TDB epoch.  Throws `std::runtime_error`
		// when an observation references an unknown observatory code.
		// A surface radar return without a body radius gets a 1 km fallback in
		// the prepared copy, with one warning per call; the input is unchanged.
		//
		// `assumeSortedInput` — when true, skip the internal sort by
		// jdTdb.  Caller asserts that `rawObs` is JD-sorted; since
		// UTC→TDB delta is monotonic over our entire ephemeris window
		// (leap-second offsets only grow), JD-sorted input ⇒
		// jdTdb-sorted output without re-sorting.  Default false for
		// backward-compat / safety.
		[[nodiscard]] std::vector<PreparedObservation> prepare(std::span<const ObsNode> rawObs,
		                                                       bool assumeSortedInput = false) const;

	   private:
		AstrometricWeightProfile m_weightProfile{AstrometricWeightProfile::Vfcc17};
		ObservationWeightingOptions m_options;
	};
}  // namespace orbital_core::obs
