#pragma once

#include "orbital_core/OrbitalCoreExport.h"

namespace orbital_core::hazard {
	inline constexpr double kPhaMoidThresholdAu{0.05};
	inline constexpr double kPhaAbsoluteMagThreshold{22.0};

	enum class PhaStatus { NotApplicable, Unknown, NotPha, Pha };

	// H and Earth MOID must refer to the asteroid solution; unavailable values
	// are represented by NaN. MOID must be non-negative.
	[[nodiscard]] ORBITAL_CORE_API PhaStatus classifyPha(bool isAsteroid, double absoluteMag,
	                                                     double earthMoidAu) noexcept;
}  // namespace orbital_core::hazard
