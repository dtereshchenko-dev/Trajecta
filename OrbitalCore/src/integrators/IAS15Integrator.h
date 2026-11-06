#pragma once

#include "IIntegrator.h"
#include <span>

// Gauss-Radau integrator with IAS15 corrector and GLOBAL step control.
// Its dimensionless epsilon is distinct from RA15's truncation-control level.
// As in REBOUND, one corrector updates state, STM and sensitivities together;
// GLOBAL selects the step from the physical three-vector only, excluding all
// variational blocks. No physical polynomial is frozen during correction.
// The shared highest-coefficient norm stops below 1e-16 or on stagnation after
// three sweeps, with a twelve-sweep cap. Non-finite corrections fail explicitly.
// Compensated coefficient updates and position/velocity increments follow the
// REBOUND approach; predictors and dense output retain the starting roundoff
// compensation. Endpoint polynomial terms are accumulated separately from small
// to large, multiplying by the step twice instead of a pre-rounded step squared.
// Node differences/transforms are rounded once from high-precision constants.
// Elapsed time is accumulated with compensation, retained across encounter
// restarts and included in stage times, dense output and the terminal step.
// Polynomial prediction uses rational Horner evaluation and the binomial shift,
// retaining the correction to the previous prediction also after the first step.
// A rejected attempt predicts again from the last accepted polynomial, not from
// the rejected endpoint. Excessive step growth clears the predictor as in REBOUND.
// Encounter control prefers the epoch-plus-offset callback, as in RA15.
class IAS15Integrator final : public IIntegrator {
   public:
	// Standard REBOUND IAS15 epsilon, fixed independently of RA15's accuracy.
	static constexpr double kEpsilon = 1.0e-9;
	// Project minimum for adaptive steps.
	// Encounter caps and shortened final steps may be smaller.
	static constexpr double kMinimumAdaptiveStepDays = 0.001;

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

	[[nodiscard]] const char* name() const noexcept override { return "IAS15"; }
};
