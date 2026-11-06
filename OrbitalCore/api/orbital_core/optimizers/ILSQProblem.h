#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <Eigen/Core>
#include <vector>

namespace orbital_core::optimizers {

	// Problem-specific residual model for least-squares correction.
	//
	// The optimizer never sees observations, light-time iteration, weights, or
	// any astrometry-specific concerns — those are baked into the implementor's
	// `residuals()` / `jacobian()`.  Pre-weighting (W = diag(1/σ²)) is the
	// implementor's responsibility, applied in-line before returning so the
	// returned residuals/Jacobian are in pre-whitened (Fisher-info-ready) form.
	class ORBITAL_CORE_API ILSQProblem {
	   public:
		virtual ~ILSQProblem() = default;

		[[nodiscard]] virtual int numParams() const noexcept = 0;
		[[nodiscard]] virtual int numResiduals() const noexcept = 0;

		// Pre-weighted residuals.  outResiduals.resize(numResiduals()).
		virtual void residuals(const std::vector<double>& params, std::vector<double>& outResiduals) const = 0;

		// Pre-weighted Jacobian.  outJ.resize(numResiduals(), numParams()).
		// The implementor decides analytic-vs-FD internally; the optimizer
		// never sees that distinction.
		virtual void jacobian(const std::vector<double>& params, Eigen::MatrixXd& outJ) const = 0;

		// Optional param-space guard (e.g. orbit-type unchanged).  Returning
		// false rejects the candidate before evaluation; the optimizer can shrink.
		[[nodiscard]] virtual bool acceptParams(const std::vector<double>& /*candidate*/) const { return true; }

		// Final reporting metric, evaluated once after the solve unless cancelled.
		// Does not control step acceptance; a nonfinite value invalidates the result.
		[[nodiscard]] virtual double displayCost(const std::vector<double>& params) const = 0;
	};

}  // namespace orbital_core::optimizers
