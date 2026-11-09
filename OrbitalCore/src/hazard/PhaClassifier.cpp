#include "orbital_core/hazard/PhaClassifier.h"

#include <cmath>

namespace orbital_core::hazard {
	PhaStatus classifyPha(bool isAsteroid, double absoluteMag, double earthMoidAu) noexcept
	{
		if (!isAsteroid)
			return PhaStatus::NotApplicable;
		if (!std::isfinite(absoluteMag) || !std::isfinite(earthMoidAu) || earthMoidAu < 0.0)
			return PhaStatus::Unknown;
		return absoluteMag <= kPhaAbsoluteMagThreshold && earthMoidAu <= kPhaMoidThresholdAu ? PhaStatus::Pha
		                                                                                     : PhaStatus::NotPha;
	}
}  // namespace orbital_core::hazard
