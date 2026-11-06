#pragma once

namespace orbital_core::hazard {

	enum class LovCoordinateSystem { Equinoctial, Keplerian, Cometary };

	[[nodiscard]] constexpr const char* lovCoordinateCode(LovCoordinateSystem system) noexcept
	{
		switch (system) {
			case LovCoordinateSystem::Equinoctial:
				return "EQU";
			case LovCoordinateSystem::Keplerian:
				return "KEP";
			case LovCoordinateSystem::Cometary:
				return "COT";
		}
		return "unknown";
	}

	[[nodiscard]] constexpr const char* lovCoordinateName(LovCoordinateSystem system) noexcept
	{
		switch (system) {
			case LovCoordinateSystem::Equinoctial:
				return "equinoctial";
			case LovCoordinateSystem::Keplerian:
				return "keplerian";
			case LovCoordinateSystem::Cometary:
				return "cometary";
		}
		return "unknown";
	}

}  // namespace orbital_core::hazard
