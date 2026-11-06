#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <vector>

namespace Calculations {

	[[nodiscard]] ORBITAL_CORE_API StateVectors gauss(const ObsNode& obs1, const ObsNode& obs2, const ObsNode& obs3);


	[[nodiscard]] ORBITAL_CORE_API std::vector<StateVectors> gaussAllRoots(const ObsNode& obs1, const ObsNode& obs2,
	                                                                       const ObsNode& obs3);

	namespace detail {

		[[nodiscard]] ORBITAL_CORE_API std::vector<double> findInitialGuesses(double a, double b, double c,
		                                                                      double start, double end, double step);

	}  // namespace detail
}  // namespace Calculations
