#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace orbital_core {

	[[nodiscard]] ORBITAL_CORE_API bool planetaryEphemerisAvailable() noexcept;
	[[nodiscard]] ORBITAL_CORE_API bool asteroidEphemerisAvailable() noexcept;

	struct RuntimeDataEntry {
		std::string id;
		std::filesystem::path path;
		bool valid = false;
		bool active = false;
		std::size_t records = 0;
		std::size_t catalogs = 0;
		double firstMjdUtc = 0;
		double lastMjdUtc = 0;
		std::string error;
	};

	// Inspects the same paths and parsers used by the runtime. File validation is
	// distinct from the state of the readers already open in this process.
	[[nodiscard]] ORBITAL_CORE_API std::vector<RuntimeDataEntry> inspectRuntimeData(
	    const std::filesystem::path& downloadedAsteroidFile = {});

}  // namespace orbital_core
