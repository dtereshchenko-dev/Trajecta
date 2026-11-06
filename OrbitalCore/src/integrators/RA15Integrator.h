#pragma once

#include "IIntegrator.h"
#include <span>

// Everhart's 15th-order Gauss-Radau integrator, ported from reference RA15V.
// One shared corrector updates state, STM and sensitivities; no state polynomial
// is frozen. As in REBOUND IAS15, convergence uses the relative change of the
// highest acceleration coefficient, stopping below 1e-16 or on stagnation after
// three sweeps, with a twelve-sweep safety limit. This is a roundoff-level
// corrector criterion, not RA15's truncation-error parameter kAccuracy.
// Adaptive step selection uses only the physical acceleration coefficients.
// RA15's existing adaptive-step retry and fixed-encounter warning behavior is
// retained if the corrector exhausts its limit; non-finite coefficients fail.
// NG terms in the shared corrector norm retain fixed scales: GM_SUN for d/dAi
// and 1/GAUSS_K for d/dDT. They do not depend on fitted amplitudes or solve masks.
// Equations, outputs and covariance retain physical AU/day units. This changes
// RA15's control policy, not its quadrature, force model or encounter stepping.
// State and elapsed time use compensated accumulation. Stage and dense-output
// polynomials include the low parts; physical/time compensation survives an
// encounter restart, while reset local STM/sensitivity blocks clear theirs.
// Accepted physical steps expose this same dense-output polynomial to the
// synchronous observer; evaluating it does not advance or restart integration.
class RA15Integrator final : public IIntegrator {
   public:
	static constexpr double kAccuracy = 1.0e-12;

	IntegrationReport integrateState(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
	                                 const IntegratedState& initial, double epochMjd,
	                                 std::span<const double> targetMjds, std::span<IntegratedState> outputs,
	                                 const IntegrationOptions& options,
	                                 const AcceptedStateStepObserver& observer = {}) const override;

	IntegrationReport integrateStateWithStm(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
	                                        const IntegratedStateWithStm& initial, double epochMjd,
	                                        std::span<const double> targetMjds,
	                                        std::span<IntegratedStateWithStm> outputs,
	                                        const IntegrationOptions& options) const override;

	IntegrationReport integrateStateWithStmAndSens(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
	                                               const IntegratedStateFull& initial, double epochMjd,
	                                               std::span<const double> targetMjds,
	                                               std::span<IntegratedStateFull> outputs,
	                                               const IntegrationOptions& options) const override;

	[[nodiscard]] const char* name() const noexcept override { return "RA15"; }
};
