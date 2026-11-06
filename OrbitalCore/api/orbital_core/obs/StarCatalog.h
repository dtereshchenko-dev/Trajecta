#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <string_view>

namespace orbital_core::obs {

	struct StarCatalog {
		char mpcCode;
		std::string_view adesName;
		std::string_view name;
	};

	[[nodiscard]] ORBITAL_CORE_API const StarCatalog* findStarCatalog(char mpcCode) noexcept;
	[[nodiscard]] ORBITAL_CORE_API const StarCatalog* findStarCatalog(std::string_view adesName) noexcept;

}  // namespace orbital_core::obs
