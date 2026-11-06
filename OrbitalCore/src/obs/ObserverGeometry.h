#pragma once

#include "obs/PreparedObservation.h"
#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/obs/RovingObserver.h"

#include <span>
#include <string>

class ObsNode;

namespace orbital_core::obs {

	struct ObserverState {
		Vec3D offsetEclAu;
		Vec3D velocityEclAuPerDay;
	};

	struct ApparentMotion {
		double raArcsecPerMinute{0.0};
		double decArcsecPerMinute{0.0};
		bool valid{false};
	};

	// WGS84 geodetic coordinates to the same parallax factors used by the
	// terrestrial station catalog. The inertial rotation is performed at epoch.
	[[nodiscard]] ObserverData rovingObserverData(const RovingObserver& position, const std::string& code);

	[[nodiscard]] ObserverState groundObserverStateEclAu(double jdTdb, double mjdUtc, const ObserverData& observer);

	[[nodiscard]] ObserverState groundObserverStateEclAuMjd(double mjdTdb, double mjdUtc, const ObserverData& observer);

	[[nodiscard]] Vec3D groundObserverOffsetEclAu(double jdTdb, double mjdUtc, const ObserverData& observer);

	[[nodiscard]] Vec3D observerOffsetAt(double jdTdb, const PreparedObservation& observation);

	[[nodiscard]] ApparentMotion apparentMotionArcsecPerMinute(const Vec3D& objectHelioEclAu,
	                                                           const Vec3D& objectVelocityEclAuPerDay,
	                                                           const PreparedObservation& observation);

}  // namespace orbital_core::obs
