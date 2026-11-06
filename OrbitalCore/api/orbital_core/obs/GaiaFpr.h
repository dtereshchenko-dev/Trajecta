#pragma once

#include "orbital_core/ObservationSet.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <filesystem>
#include <stop_token>
#include <string>
#include <string_view>
#include <vector>

namespace orbital_core::obs {
	// Gaia Archive TAP JSON, kept separately from transport for reproducibility.
	[[nodiscard]] ORBITAL_CORE_API std::vector<GaiaFprObservation> parseGaiaFprResponse(std::string_view json,
	                                                                                    int minorPlanetNumber);

	enum class GaiaFprSource { Archive, Cache };

	struct GaiaFprDownload {
		std::vector<GaiaFprObservation> observations;
		GaiaFprSource source{GaiaFprSource::Archive};
		std::string notice;
	};

	// Cache is used only after an unavailable network, never for an invalid response.
	[[nodiscard]] ORBITAL_CORE_API GaiaFprDownload downloadGaiaFpr(int minorPlanetNumber,
	                                                               const std::filesystem::path& cacheDirectory = {},
	                                                               std::stop_token cancel = {});
}  // namespace orbital_core::obs
