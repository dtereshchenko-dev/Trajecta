#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <string_view>

namespace orbital_core {

	// Coordinates are in AU; returned velocity is in AU/day, in J2000 axes.
	ORBITAL_CORE_API bool earthPckState(double mjdTdb, const Vec3D& terrestrial, Vec3D& position, Vec3D& velocity);
	ORBITAL_CORE_API bool earthPckPole(double mjdTdb, Vec3D& pole);
	ORBITAL_CORE_API bool earthPckToTerrestrial(double mjdTdb, const Vec3D& equatorial, Vec3D& terrestrial);
	// Position and velocity must use matching distance units, with velocity per second.
	ORBITAL_CORE_API bool earthPckToTerrestrialState(double mjdTdb, const Vec3D& position, const Vec3D& velocity,
	                                                 Vec3D& terrestrialPosition, Vec3D& terrestrialVelocity);
	ORBITAL_CORE_API bool earthPckLoaded(std::string_view filename);

}  // namespace orbital_core
