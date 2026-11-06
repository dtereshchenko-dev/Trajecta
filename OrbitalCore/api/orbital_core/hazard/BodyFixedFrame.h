#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/Planet.h"
#include <array>
#include <string>
#include <string_view>

class IEphemerisSource;

// J2000 ecliptic -> body-fixed for impact-corridor lat/lon.
// Earth: optical IAU 1976/1980 precession-nutation + UT1/GAST, no polar motion.
// Moon: DE440 principal axes -> ME440 aligned to ME421.
// Other planets: full IAU rotation series in NAIF pck00011, including periodic
// terms. These reference frames do not imply a terrain model or equal accuracy.
namespace orbital_core::hazard {

	struct BodyFixedTransform {
		std::array<Vec3D, 3> equatorialRows{};
		std::string_view frame;
		std::string_view model;
		bool valid{};

		[[nodiscard]] Vec3D applyEquatorial(const Vec3D& vector) const noexcept
		{
			return {equatorialRows[0].dot(vector), equatorialRows[1].dot(vector), equatorialRows[2].dot(vector)};
		}
	};

	// MJD TDB at the public boundary; Earth orientation converts locally to TT/UT1.
	// A supplied source takes precedence over the default. Missing lunar data,
	// unsupported DE frames and non-finite epochs return valid=false, never an
	// unlabeled approximate lunar orientation.
	[[nodiscard]] ORBITAL_CORE_API BodyFixedTransform
	bodyFixedTransformMjd(Planet planet, double mjdTdb, IEphemerisSource* source = nullptr) noexcept;

	// Same transform for an ecliptic vector. Invalid transforms return NaNs.
	[[nodiscard]] ORBITAL_CORE_API Vec3D icrfToBodyFixedMjd(Planet planet, const Vec3D& rPlanetCentricEclKm,
	                                                        double mjdTdb, IEphemerisSource* source = nullptr) noexcept;

	// Includes the rate of the rotating frame; input and output are km/s.
	[[nodiscard]] ORBITAL_CORE_API Vec3D bodyFixedVelocityMjd(Planet planet, const Vec3D& rPlanetCentricEclKm,
	                                                          const Vec3D& vPlanetCentricEclKmS, double mjdTdb,
	                                                          IEphemerisSource* source = nullptr) noexcept;

	// Provenance of the Earth transformation below, distinct from the selected
	// EOP profile used by the dynamical pole. Capture with computed states.
	[[nodiscard]] ORBITAL_CORE_API std::string earthBodyFixedFrameModel();

	struct ORBITAL_CORE_API LatLon {
		double latitudeDeg{};   // [-90, +90]
		double longitudeDeg{};  // [-180, +180]
		bool valid{false};
	};

	// Planet-centric ecliptic-J2000 (km) → planet-centric body-fixed (km).
	[[nodiscard]] ORBITAL_CORE_API Vec3D icrfToBodyFixed(Planet planet, const Vec3D& rPlanetCentricEclKm,
	                                                     double jdTdb) noexcept;

	// Body-fixed cartesian → lat/lon (deg), spherical fallback.
	[[nodiscard]] ORBITAL_CORE_API LatLon latLonFromBodyFixed(const Vec3D& rBodyFixedKm) noexcept;

	// Body-fixed cartesian → geodetic lat/lon (deg) on the planet reference
	// surface when an ellipsoid is defined; spherical fallback otherwise.
	[[nodiscard]] ORBITAL_CORE_API LatLon latLonFromBodyFixed(Planet planet, const Vec3D& rBodyFixedKm) noexcept;
}  // namespace orbital_core::hazard
