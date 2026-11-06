#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <vector>

namespace Calculations {

	// ---------------------------------------------------------------------------
	// Herget's method of initial orbit determination.
	//
	// Two-parameter fit: topocentric distances ρ₁ at obs.front() and ρₙ at
	// obs.back().  Endpoints feed a Lambert solve (universal variables) → V₁;
	// the resulting two-body orbit is propagated (Kepler UV) to every
	// intermediate observation and the sky-plane residuals are driven down by
	// Gauss-Newton on (ρ₁, ρₙ).
	//
	// More tolerant than Gauss on short arcs (days–weeks) where the 8th-
	// degree polynomial has multiple/ill-conditioned roots, and uses ALL
	// observations (Gauss uses only 3).
	//
	// Units internally: AU; μ = 1; time = τ = Constants::GAUSS_K·Δt — matches
	// Gauss.h / Vaisala.h output conventions.  Observer→target sign follows the
	// Vaisala convention in this codebase: heliocentric target position is
	//   R = L̂·ρ − obsVector   (where obsVector points observer → Sun).
	//
	// Parameters:
	//   obs       JD-sorted observations, size ≥ 3.
	//   rho1Init  Initial ρ at obs.front() (AU).  Pass ≤ 0 to grid-search
	//             in [0.5, 5.0] AU at 0.25 AU spacing.
	//   rhoNInit  Initial ρ at obs.back()  (AU).  Same grid-search fallback.
	//   maxIter   Upper bound on Gauss-Newton iterations.
	//
	// Returns StateVectors at obs.front().JD; isValid() == false on failure
	// (too few obs, zero arc, grid search found no valid endpoint pair,
	// Lambert/Kepler divergence, line search stalled).
	//
	// Thread-safe: pure function, no shared state.
	//
	// Reference: Herget 1939/65; Bate-Mueller-White (universal Kepler);
	// Vallado 2013 Alg. 8 (Kepler UV) & Alg. 58 (Lambert UV).
	// ---------------------------------------------------------------------------
	[[nodiscard]] ORBITAL_CORE_API StateVectors herget(const std::vector<ObsNode>& obs, double rho1Init = -1.0,
	                                                   double rhoNInit = -1.0, int maxIter = 25);

	// Two-endpoint transfer used to improve the initial orbit.
	// Distances are topocentric AU at obs1/obs2; returned state is at obs1.JD
	// in the same canonical units as herget()/gauss()/vaisala().
	[[nodiscard]] ORBITAL_CORE_API StateVectors transferOrbitFromEndpointDistances(const ObsNode& obs1,
	                                                                               const ObsNode& obs2, double rho1,
	                                                                               double rho2);

}  // namespace Calculations
