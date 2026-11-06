#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <cstddef>
#include <stop_token>
#include <string>
#include <string_view>

namespace orbital_core::obs {

	struct ObservatoryUpdate {
		std::string contents;
		std::size_t updated{};
		std::size_t unmapped{};
		bool changed = true;
	};
	ORBITAL_CORE_API void validateOpticalObservatoryCatalog(std::string_view contents);

	[[nodiscard]] ORBITAL_CORE_API ObservatoryUpdate prepareOpticalObservatoryUpdate(
	    std::string_view source, std::string_view current);
	[[nodiscard]] ORBITAL_CORE_API ObservatoryUpdate prepareRadarObservatoryUpdate(
	    std::string_view source, std::string_view current);
	[[nodiscard]] ORBITAL_CORE_API ObservatoryUpdate downloadOpticalObservatoryUpdate(
	    std::string_view current, std::stop_token cancel = {});
	[[nodiscard]] ORBITAL_CORE_API ObservatoryUpdate downloadRadarObservatoryUpdate(
	    std::string_view current, std::stop_token cancel = {});
}
