#pragma once

#include "NBodyPropagator.h"
#include "orbital_core/hazard/BoundaryCrossing.h"
#include "orbital_core/hazard/Planet.h"
#include "orbital_core/FitPipeline.h"

#include <array>
#include <optional>
#include <stop_token>

namespace orbital_core::hazard::boundary {

	struct SearchResult {
		bool beforeInterval{};
		std::optional<BoundaryCrossing> entry;
	};

	SearchResult search(const std::array<double, 6>& state, double epochMjdTdb, const NonGravModel& nonGrav,
	                    Planet target, double startMjdTdb, double endMjdTdb, double altitudeKm,
	                    NBodyPropagator& propagator, IEphemerisSource& ephemeris, std::stop_token cancel);

	std::optional<NominalEarthEntry> nearbyNominalEarthEntry(const FitResult& fit, const NBodySettings& forces,
	                                                         std::stop_token cancel);
}  // namespace orbital_core::hazard::boundary
