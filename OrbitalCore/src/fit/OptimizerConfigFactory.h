#pragma once

#include "OrbitFitConfig.h"
#include "orbital_core/FitConfig.h"

namespace orbital_core::fit {

	[[nodiscard]] OrbitOptimizerConfig makeOrbitOptimizerConfig(const FitConfig& cfg) noexcept;

}
