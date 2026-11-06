#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/OrbitalCoreExport.h"

namespace Calculations {

	// ---------------------------------------------------------------------------
	// Väisälä's initial orbit determination.
	//
	// Assumes a (near-)circular orbit (e ≈ 0) and fits a single parameter —
	// heliocentric distance r — by grid-searching for the value that best
	// reproduces the angular separation between two observations under Kepler
	// circular-orbit dynamics.  The associated velocity is recovered from the
	// 2nd-order Lagrange f&g series (same convention as Gauss IOD).
	//
	// Designed for distant / slow-moving objects (TNOs, Centaurs, long-period
	// comets) where Gauss's 8th-degree polynomial becomes ill-conditioned and
	// Herget's Lambert solver tends to oscillate on the near-parabolic branch.
	// Not appropriate for near-Earth objects where e is substantial.
	//
	// Parameters:
	//   obs1   Earliest observation.
	//   obs2   Latest observation.  Must differ from obs1 by ≥ 0.01 days.
	//   rMin   Lower heliocentric-distance search bound (AU).
	//   rMax   Upper heliocentric-distance search bound (AU).
	//
	// Returns StateVectors at obs1.JD; isValid() == false on failure (arc too
	// short, no valid circular solution in [rMin, rMax], etc.).
	//
	// Thread-safe: pure function, no shared state.
	// ---------------------------------------------------------------------------
	[[nodiscard]] ORBITAL_CORE_API StateVectors vaisala(const ObsNode& obs1, const ObsNode& obs2, double rMin = 5.0,
	                                                    double rMax = 200.0);

}  // namespace Calculations
