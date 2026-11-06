#pragma once

#include "RasterSource.h"
#include <functional>
#include <span>

namespace strewn::sources {

	struct GfsSelection {
		std::int64_t cycleUnix{};
		std::int64_t validUnix{};
		std::string url;
		bool indexed{};
	};

	GfsSelection selectGfs(double eventMjdUtc, double nowMjdUtc, std::optional<double> cycleMjdUtc = {});
	std::vector<std::filesystem::path> downloadGfsRanges(
	    std::span<const trajecta::http::ByteRange> ranges,
	    const std::function<std::filesystem::path(trajecta::http::ByteRange, std::stop_token)>& download,
	    std::stop_token cancel);
	WeatherFields obtainGfs(const AutomaticEnvironment& options, double eventMjdUtc, const Bounds& bounds,
	                        const std::filesystem::path& directory, std::stop_token cancel, Json& provenance);
	AtmosphereGrid atmosphereFromGfs(const WeatherFields& fields, const std::vector<double>& geoidM);

}  // namespace strewn::sources
