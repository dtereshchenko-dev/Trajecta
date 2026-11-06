#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <vector>

namespace Calculations {

	// ---------------------------------------------------------------------------
	// Result of a multi-method IOD attempt — caller uses `state` + `epoch` as
	// LM seed, `method` for display/logging, `costRad2` / `score` for
	// diagnostics.  `state.isValid()` is false ⇒ ALL methods failed.
	// ---------------------------------------------------------------------------
	struct ORBITAL_CORE_API InitialOrbitResult {
		StateVectors state;
		double epoch{0.0};  // JD TDB
		const char* method{nullptr};
		double costRad2{0.0};  // keplerResidualCost on `allObs`
		double score{0.0};     // initial orbit quality score
		std::size_t subarcStart{static_cast<std::size_t>(-1)};
		std::size_t subarcEnd{static_cast<std::size_t>(-1)};
	};

	// ---------------------------------------------------------------------------
	// Multi-method IOD orchestrator.  Pipeline:
	//
	//   arclen = min(span(groundObs), 730 d)
	//
	//   while best.score > kAcceptableScore and arclen > 0.001 d:
	//     (start, end) = lookForBestSubarc(groundObs, arclen)
	//     subarc = groundObs[start..end]
	//
	//     # Classical seeds on a triplet drawn from the subarc:
	//     for each surviving Gauss polynomial root → score → Herget-refine
	//     Herget(triplet)
	//     Vaisala(triplet endpoints)
	//
	//     # Pseudo-r scans on the subarc (only if subSpan < 230 d):
	//     Vaisala-scan: r = 0.1..200 AU × 1.2  (flip-out at score>5e4)
	//     Herget-scan:  ρ = ρ_min..5 AU × 1.2  (no flip-out — near-parabolic
	//                   comets need to scan past intermediate Lambert blow-ups)
	//
	//     arclen *= 0.7
	//
	// The search is bounded by the decreasing arclen and finite scan ranges, not
	// elapsed wall-clock time. Machine load and logging cannot truncate candidates.
	// If the subarc bound is exhausted, the best orbit found so far is returned.
	//
	// Each candidate is scored using weighted RMS/2, a small eccentricity term,
	// hyperbolic and high-inclination penalties, an empirical
	// inclination-vs-a scattergram adjustment for main-belt-like orbits,
	// and the physical sun-surface rejection for bound orbits.
	//
	// The selected winner carries the subarc mask that produced it.  This is
	// used as the starting mask for arc extension.
	//
	// Caller responsibilities:
	//   `groundObs` — JD-sorted astrometric observations used for IOD.
	//                 Subarc selection happens here, NOT upstream — this
	//                 is what lets the dynamic arclen-shrink find a
	//                 Gauss-convergent window for multi-apparition objects.
	//                 ≥ 3 entries.
	//   `allObs`   — JD-sorted full observation set used for scoring.
	//
	// Returns InitialOrbitResult with state.isValid() == false if every
	// method failed.
	[[nodiscard]] ORBITAL_CORE_API InitialOrbitResult bestInitialOrbit(
	    const std::vector<ObsNode>& groundObs, const std::vector<ObsNode>& allObs,
	    orbital_core::AstrometricWeightProfile weightProfile = orbital_core::AstrometricWeightProfile::Vfcc17);
}  // namespace Calculations
