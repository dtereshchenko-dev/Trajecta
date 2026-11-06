#pragma once

#include "IIntegrator.h"

// NBodyPropagator keeps heliocentric inputs, outputs, observers and force APIs.
// This adapter integrates barycentric states, avoiding numerical integration
// of the second derivative of the piecewise DE solar ephemeris. Translation
// depends only on epoch, so state/NG variations retain their original units
// and require no Jacobian change. NG forcing still uses heliocentric RTN.
class BarycentricIntegrator final : public IIntegrator {
   public:
	explicit BarycentricIntegrator(IIntegrator_p engine);
	virtual ~BarycentricIntegrator() = default;

	IntegrationReport integrateState(const IIntegrationDynamics&, const NonGravModel&, const IntegratedState&, double,
	                                 std::span<const double>, std::span<IntegratedState>, const IntegrationOptions&,
	                                 const AcceptedStateStepObserver& = {}) const override;

	IntegrationReport integrateStateWithStm(const IIntegrationDynamics&, const NonGravModel&,
	                                        const IntegratedStateWithStm&, double, std::span<const double>,
	                                        std::span<IntegratedStateWithStm>,
	                                        const IntegrationOptions&) const override;

	IntegrationReport integrateStateWithStmAndSens(const IIntegrationDynamics&, const NonGravModel&,
	                                               const IntegratedStateFull&, double, std::span<const double>,
	                                               std::span<IntegratedStateFull>,
	                                               const IntegrationOptions&) const override;

	[[nodiscard]] const char* name() const noexcept override;

   private:
	IIntegrator_p m_engine;
};
