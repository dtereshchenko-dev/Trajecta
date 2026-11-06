#pragma once

#include "orbital_core/nongrav/NonGravModel.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <vector>

// ---------------------------------------------------------------------------
// NonGravProfile
//   * `NonGravModel` is **physics** — given (model_type, A1, A2, A3, DT)
//     it computes g(r), acceleration, ∂a/∂A_i, ∂a/∂DT.  Useful anywhere
//     a propagator or N-body integrator needs the force, regardless of
//     whether anyone is fitting against it.
//
//   * `NonGravProfile` packs and unpacks the active model parameters for
//     least-squares fitting. The propagator only needs `NonGravModel`.
//
// ---------------------------------------------------------------------------
class ORBITAL_CORE_API NonGravProfile {
   public:
	explicit NonGravProfile(NonGravModel& model) noexcept;
	~NonGravProfile() = default;

	NonGravProfile(const NonGravProfile&) = delete;
	NonGravProfile& operator=(const NonGravProfile&) = delete;
	NonGravProfile(NonGravProfile&&) noexcept = default;
	NonGravProfile& operator=(NonGravProfile&&) noexcept = default;

	[[nodiscard]] NonGravModel& model() noexcept;
	[[nodiscard]] const NonGravModel& model() const noexcept;

	// ---- LM parameter-vector adapters ----
	//
	// Append the model's fittable parameters (initial guesses) and matching
	// numerical-derivative steps to the optimiser's working vectors.  Adds
	// `model().numFitParams()` entries; no-op when the model is disabled.
	// Reads current A1/A2/A3/DT from the wrapped model — call `setParams`
	// / `setDT` on the model first if a warm-start is needed.
	void extendParams(std::vector<double>& params, std::vector<double>& deltaVals) const;

	// Read fitted values back from `params[offset .. offset + numFitParams - 1]`
	// into the wrapped model (writes A1/A2/A3, plus DT for sublimation models).
	// Yarkovsky variant pins A1 = A3 = 0 by construction.
	void updateFromParams(const std::vector<double>& params, int offset) const;

	// Static overload — write LM-fit params into an arbitrary NonGravModel.
	// Used by OMP-parallel Jacobian columns where each thread snapshots its
	// own local NonGravModel copy and updates it from a perturbed parameter
	// vector without mutating the shared model.
	static void updateFromParams(NonGravModel& model, const std::vector<double>& params, int offset);

   private:
	NonGravModel* p_model;
};
