#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <optional>
#include <stop_token>
#include <string_view>

namespace orbital_core::obs {
	// Radius in kilometres from a positive SBDB diameter, if one is available.
	[[nodiscard]] ORBITAL_CORE_API std::optional<double> parseSbdbRadiusKm(std::string_view json,
	                                                                       int minorPlanetNumber);

	[[nodiscard]] ORBITAL_CORE_API std::optional<double> downloadSbdbRadiusKm(int minorPlanetNumber,
	                                                                          std::stop_token cancel = {});
}  // namespace orbital_core::obs
