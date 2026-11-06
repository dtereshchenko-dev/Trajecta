#pragma once

#include "SourceCache.h"
#include <map>

namespace strewn::sources {

	using Bounds = std::array<double, 4>;
	// Windows deployments carry share/gdal and share/proj beside this module.
	// Explicit GDAL_DATA / PROJ_DATA / PROJ_LIB settings retain precedence.
	void initializeRasters();

	struct WeatherFields {
		std::vector<double> latitudes;
		std::vector<double> longitudes;
		std::map<int, std::map<std::string, std::vector<double>>> levels;
		std::int64_t validUnix{};
		std::int64_t cycleUnix{};
	};

	std::vector<double> rasterSamples(const std::filesystem::path& path, const std::vector<double>& latitudes,
	                                  const std::vector<double>& longitudes, bool clampEdges = false);
	void readGrib(const std::filesystem::path& path, const Bounds& bounds, WeatherFields& result);

}  // namespace strewn::sources
