#pragma once

#include "orbital_core/IEphemerisSource.h"

#include <array>
#include <cstdint>
#include <string_view>

// Shared enum for the NEA hazard analysis suite (MOID, impact-risk,
// LineOfVariations, ImpactMonitor).  Deliberately narrower than
// `EphemerisPlanet` (which includes Sun and Moon) — hazard metrics are
// only meaningful for the 8 planets, so the type system enforces that
// at the API boundary instead of relying on runtime checks.
namespace orbital_core::hazard {
	enum class Planet : std::uint8_t {
		Mercury = 0,
		Venus,
		Earth,
		Mars,
		Jupiter,
		Saturn,
		Uranus,
		Neptune,
		Moon,
		Count_
	};

	inline constexpr std::size_t kPlanetCount = static_cast<std::size_t>(Planet::Count_);

	[[nodiscard]] constexpr std::string_view longLabel(Planet p) noexcept
	{
		switch (p) {
			case Planet::Mercury:
				return "Mercury";
			case Planet::Venus:
				return "Venus";
			case Planet::Earth:
				return "Earth";
			case Planet::Mars:
				return "Mars";
			case Planet::Jupiter:
				return "Jupiter";
			case Planet::Saturn:
				return "Saturn";
			case Planet::Uranus:
				return "Uranus";
			case Planet::Neptune:
				return "Neptune";
			case Planet::Moon:
				return "Moon";
			case Planet::Count_:
				break;
		}
		return "??";
	}

	inline constexpr std::array<Planet, kPlanetCount> kAllPlanets = {
	    Planet::Mercury, Planet::Venus,  Planet::Earth,   Planet::Mars, Planet::Jupiter,
	    Planet::Saturn,  Planet::Uranus, Planet::Neptune, Planet::Moon,
	};

	[[nodiscard]] constexpr EphemerisPlanet toEphemerisPlanet(Planet p) noexcept
	{
		switch (p) {
			case Planet::Mercury:
				return EphemerisPlanet::Mercury;
			case Planet::Venus:
				return EphemerisPlanet::Venus;
			case Planet::Earth:
				return EphemerisPlanet::Earth;
			case Planet::Mars:
				return EphemerisPlanet::Mars;
			case Planet::Jupiter:
				return EphemerisPlanet::Jupiter;
			case Planet::Saturn:
				return EphemerisPlanet::Saturn;
			case Planet::Uranus:
				return EphemerisPlanet::Uranus;
			case Planet::Neptune:
				return EphemerisPlanet::Neptune;
			case Planet::Moon:
				return EphemerisPlanet::Moon;
			case Planet::Count_:
				break;
		}
		return EphemerisPlanet::Earth;
	}
}  // namespace orbital_core::hazard
