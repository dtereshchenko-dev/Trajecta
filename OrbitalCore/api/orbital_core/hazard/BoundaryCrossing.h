#pragma once

#include "orbital_core/Vec3D.h"

namespace orbital_core::hazard {

	struct BoundaryCrossing {
		double mjdTdb{};
		Vec3D relativePositionKm;
		Vec3D relativeVelocityKmS;
		double latitudeDeg{};
		double longitudeDeg{};
	};
}  // namespace orbital_core::hazard
