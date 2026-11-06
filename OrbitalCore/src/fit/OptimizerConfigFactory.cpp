#include "fit/OptimizerConfigFactory.h"

namespace orbital_core::fit {
	OrbitOptimizerConfig makeOrbitOptimizerConfig(const FitConfig& config) noexcept
	{
		OrbitOptimizerConfig result;
		result.nbody = config.nbody;
		result.nonGrav = config.nonGrav;
		result.useWeights = config.useWeights;
		result.weightProfile = config.weightProfile;
		result.obsWeighting = config.obsWeighting;
		result.gaiaPhotocenterRadiusKm = config.gaiaPhotocenterRadiusKm;
		return result;
	}
}  // namespace orbital_core::fit
