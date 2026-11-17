#pragma once

#include "orbital_core/hazard/SemilinearContours.h"
#include "orbital_core/hazard/Planet.h"
#include "orbital_core/FitPipeline.h"
#include "NBodyPropagator.h"
#include <stop_token>

namespace orbital_core::hazard {
	SemilinearContours computeSemilinearContours(const std::array<double, 6>& state, double epoch,
	                                             const NonGravModel& nonGrav, const FitParameterCovariance& covariance,
	                                             const std::optional<BoundaryCrossing>& nominal, Planet target,
	                                             double start, double end, double altitude, unsigned workers,
	                                             const NBodyPropagator& propagator, IEphemerisSource& ephemeris,
	                                             std::stop_token cancel);
}
