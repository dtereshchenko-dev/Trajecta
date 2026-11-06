#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <string>

namespace orbital_core {
	using ObservationId = std::string;
	[[nodiscard]] ORBITAL_CORE_API ObservationId newObservationId();
}  // namespace orbital_core
