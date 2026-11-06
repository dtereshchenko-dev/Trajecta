#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/hazard/Planet.h"

// Per-planet physical constants in km and km^3/s^2. Hazard code works in km
// units at close-encounter scales, separate from the AU/day dynamical model.
namespace orbital_core::hazard {

	struct ORBITAL_CORE_API PlanetData {
		double muKm3PerS2{};  // GM in km³/s²
		double radiusEqKm{};  // equatorial radius (km), used as impact sphere
		double flattening{};  // reference ellipsoid flattening; 0 => spherical
		double j2{};          // zonal harmonic; reserved for J₂ corrections
	};

	[[nodiscard]] ORBITAL_CORE_API const PlanetData& getPlanetData(Planet p) noexcept;

	[[nodiscard]] inline double planetRadiusKm(Planet p) noexcept { return getPlanetData(p).radiusEqKm; }

	[[nodiscard]] inline double planetMuKm3PerS2(Planet p) noexcept { return getPlanetData(p).muKm3PerS2; }

	[[nodiscard]] inline double planetFlattening(Planet p) noexcept { return getPlanetData(p).flattening; }
}  // namespace orbital_core::hazard
