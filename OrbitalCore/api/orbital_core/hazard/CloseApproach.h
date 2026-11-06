#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/Planet.h"

#include <limits>
#include <vector>

class IEphemerisSource;
class NBodyPropagator;
class NonGravModel;
class Orbit;

// Find local minima of planet-centric distance over [t0, t0+horizon].
// Coarse-grid propagation brackets inbound/outbound radial velocity and refines
// the time of closest approach by solving r_rel · v_rel / |r_rel| = 0.
// Filters events with distance >= thresholdAu.
namespace orbital_core::hazard {

	struct ORBITAL_CORE_API CloseApproachEvent {
		double jd{};                   // TDB Julian Date at TCA
		double distAu{};               // Min distance to planet (AU)
		double vInfKmS{};              // Asymptotic relative speed (km/s)
		Planet planet{Planet::Earth};  // Body whose neighbourhood was scanned
		Vec3D rTarget{};               // Heliocentric AU
		Vec3D vTarget{};               // Heliocentric AU/day
		Vec3D rPlanet{};               // Heliocentric AU
		Vec3D vPlanet{};               // Heliocentric AU/day
		double localMoidAu{std::numeric_limits<double>::quiet_NaN()};
		int encounterIndex{-1};  // Contiguous influence-region passage in this scan
	};

	// Sampled influence-region passages, including those without a minimum.
	struct CloseApproachPassage {
		int encounterIndex{-1};
		double firstInsideJd{};
		double lastInsideJd{};
	};

	struct ORBITAL_CORE_API CloseApproachConfig {
		double horizonYears{100.0};
		double lookbackYears{0.0};
		double thresholdAu{Constants::EARTH_CLOSE_APPROACH_REPORT_AU};  // reference dmeacontr.
		double encounterRadiusAu{Constants::EARTH_TARGET_PLANE_AU};     // reference dmea for NEOs.
		double sampleStepDay{1.0};
		Planet planet{Planet::Earth};
	};

	[[nodiscard]] ORBITAL_CORE_API std::vector<CloseApproachEvent> findCloseApproaches(
	    const Orbit& orbit, const NBodyPropagator& propagator, const CloseApproachConfig& cfg = {},
	    IEphemerisSource* eph = nullptr);

	[[nodiscard]] ORBITAL_CORE_API std::vector<CloseApproachEvent> findCloseApproachesFromState(
	    const Vec3D& r0Au, const Vec3D& v0AuPerDay, double epochTdb, const NBodyPropagator& propagator,
	    const CloseApproachConfig& cfg = {}, IEphemerisSource* eph = nullptr);

	[[nodiscard]] ORBITAL_CORE_API std::vector<CloseApproachEvent> findCloseApproachesFromState(
	    const Vec3D& r0Au, const Vec3D& v0AuPerDay, double epochTdb, const NBodyPropagator& propagator,
	    const NonGravModel& nonGrav, const CloseApproachConfig& cfg = {}, IEphemerisSource* eph = nullptr);

	[[nodiscard]] ORBITAL_CORE_API std::vector<CloseApproachEvent> findCloseApproachesFromState(
	    const Vec3D& r0Au, const Vec3D& v0AuPerDay, double epochTdb, const NBodyPropagator& propagator,
	    const NonGravModel& nonGrav, const CloseApproachConfig& cfg, IEphemerisSource* eph,
	    std::vector<CloseApproachPassage>& passages);

}  // namespace orbital_core::hazard
