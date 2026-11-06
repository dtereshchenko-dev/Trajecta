#pragma once

#include "orbital_core/IEphemerisSource.h"
#include <cstdint>

namespace orbital_core {
	[[nodiscard]] constexpr std::uint16_t planetaryBodyMask(EphemerisPlanet body) noexcept
	{
		switch (body) {
			case EphemerisPlanet::Mercury:
				return 1u << 0;
			case EphemerisPlanet::Venus:
				return 1u << 1;
			case EphemerisPlanet::Earth:
				return 1u << 2;
			case EphemerisPlanet::Mars:
				return 1u << 3;
			case EphemerisPlanet::Jupiter:
				return 1u << 4;
			case EphemerisPlanet::Saturn:
				return 1u << 5;
			case EphemerisPlanet::Uranus:
				return 1u << 6;
			case EphemerisPlanet::Neptune:
				return 1u << 7;
			case EphemerisPlanet::Moon:
				return 1u << 11;
			case EphemerisPlanet::Pluto:
				return 1u << 12;
			default:
				return 0;
		}
	}
}  // namespace orbital_core
