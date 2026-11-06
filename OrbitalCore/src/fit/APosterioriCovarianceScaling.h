#pragma once

#include "orbital_core/OrbitalCoreExport.h"

namespace orbital_core::fit {

	[[nodiscard]] ORBITAL_CORE_API double aPosterioriVarianceFactor(double chiSquare, int numUsedResiduals,
	                                                                int numSolvedParameters) noexcept;

}  // namespace orbital_core::fit
