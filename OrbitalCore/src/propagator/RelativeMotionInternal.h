#pragma once

#include "NBodyPropagator.h"
#include "orbital_core/FitConfig.h"

namespace orbital_core::detail {
	NBodyPropagator relativeMotionPropagator(const NBodySettings& forces, IEphemerisSource& source);
}
