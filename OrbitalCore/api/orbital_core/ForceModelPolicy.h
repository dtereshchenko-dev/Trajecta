#pragma once

#include <cstdint>

namespace orbital_core {

	enum class InactiveInteriorMassPolicy : std::uint8_t {
		ExplicitOnly = 0,
		SmoothCentralMass = 1,
	};

}  // namespace orbital_core
