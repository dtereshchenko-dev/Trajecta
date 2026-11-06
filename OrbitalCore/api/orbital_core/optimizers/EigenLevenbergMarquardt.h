#pragma once
#include "orbital_core/optimizers/ILSQProblem.h"
#include <functional>
#include <vector>

namespace orbital_core::optimizers {

	class ORBITAL_CORE_API EigenLevenbergMarquardt final {
	   public:
		struct Result {
			std::vector<double> params;
			double finalCost{0.0};
			int iterations{0};
			bool converged{false};
			bool cancelled{false};
		};

		[[nodiscard]] Result optimize(const std::vector<double>& initialParams, ILSQProblem& problem,
		                              std::function<bool()> isCancelled = {});
		// Returning true signals an objective change at an accepted state; LM is reinitialized.
		[[nodiscard]] Result optimize(const std::vector<double>& initialParams, ILSQProblem& problem,
		                              std::function<bool()> isCancelled,
		                              std::function<bool(const std::vector<double>&)> updateSelection);
	};

}  // namespace orbital_core::optimizers
