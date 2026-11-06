#include "AutoEnvironment.h"
#include "JobFiles.h"
#include "VectorMath.h"
#include <algorithm>
#include <cmath>
#include <format>
#include <set>
#include <sstream>

namespace strewn::sources {
	namespace {

		void validateBounds(const Bounds& bounds)
		{
			for (double v : bounds)
				if (!std::isfinite(v))
					throw std::invalid_argument("Nonfinite environment bounds");
			if (bounds[0] >= bounds[1] || bounds[2] >= bounds[3] || bounds[0] < -85 || bounds[1] > 85 ||
			    bounds[2] < -179.5 || bounds[3] > 179.5 || bounds[3] - bounds[2] > 30 || bounds[1] - bounds[0] > 30)
				throw std::invalid_argument("Automatic environment requires regional, non-polar, non-dateline bounds");
		}

		std::vector<double> axis(double low, double high, double step)
		{
			const int first = static_cast<int>(std::floor(low / step)), last = static_cast<int>(std::ceil(high / step));
			if (last - first > 10000)
				throw std::invalid_argument("Environment axis is too large");
			std::vector<double> values;
			for (int i = first; i <= last; ++i)
				values.push_back(i * step);
			return values;
		}

		TerrainGrid terrain(const AutomaticEnvironment& options, const Bounds& bounds,
		                    const std::filesystem::path& directory, const std::filesystem::path& geoid,
		                    std::stop_token cancel, Json& provenance)
		{
			TerrainGrid result;
			result.latitudesDeg = axis(bounds[0], bounds[1], options.terrainStepDeg);
			result.longitudesDeg = axis(bounds[2], bounds[3], options.terrainStepDeg);
			const auto count = result.latitudesDeg.size() * result.longitudesDeg.size();
			if (count > 200000)
				throw std::invalid_argument(
				    "DEM request exceeds 200000 samples; increase terrain_step_deg or narrow bounds");
			provenance = {
			    {"provider", "Copernicus DEM GLO-90 / AWS public COG / 2021 release"},
			    {"input_height_datum", "EGM2008"},
			    {"conversion", "h_WGS84 = H_EGM2008 + N_EGM2008"},
			    {"sample_step_deg", options.terrainStepDeg},
			    {"note", "Sampled surface model, not bare-earth DTM; effective grid resolution is the configured step"},
			    {"resampling",
			     "bilinear within each native tile; missing shared edge held constant within one native pixel"},
			    {"raw_sources", Json::array()},
			    {"ocean_tiles", Json::array()}};
			const std::string bucket = "https://copernicus-dem-90m.s3.amazonaws.com/";
			const auto list = fetch(bucket + "tileList.txt", directory, options, cancel);
			std::set<std::string> available;
			std::istringstream lines(readBytes(list, options.maxDownloadBytes));
			for (std::string name; std::getline(lines, name);) {
				while (!name.empty() && (name.back() == '\r' || name.back() == '/'))
					name.pop_back();
				if (!name.starts_with("Copernicus_DSM_COG_30_"))
					throw std::runtime_error("Invalid Copernicus tile list");
				available.insert(name);
			}
			if (available.size() < 10000)
				throw std::runtime_error("Incomplete Copernicus tile inventory");
			provenance["tile_inventory"] = list.filename().string();
			std::map<int, std::vector<std::size_t>> latitudeGroups, longitudeGroups;
			for (std::size_t i = 0; i < result.latitudesDeg.size(); ++i)
				latitudeGroups[static_cast<int>(std::ceil(result.latitudesDeg[i])) - 1].push_back(i);
			for (std::size_t j = 0; j < result.longitudesDeg.size(); ++j)
				longitudeGroups[static_cast<int>(std::floor(result.longitudesDeg[j]))].push_back(j);
			if (latitudeGroups.size() * longitudeGroups.size() > 100)
				throw std::invalid_argument("DEM request exceeds 100 one-degree tiles");
			result.heightsWgs84M.resize(count, 0);
			for (const auto& [lat, rows] : latitudeGroups)
				for (const auto& [lon, cols] : longitudeGroups) {
					if (cancel.stop_requested())
						throw trajecta::http::Cancelled{};
					const auto name =
					    std::format("Copernicus_DSM_COG_30_{}{:02}_00_{}{:03}_00_DEM", lat < 0 ? 'S' : 'N',
					                std::abs(lat), lon < 0 ? 'W' : 'E', std::abs(lon));
					if (!available.contains(name)) {
						provenance["ocean_tiles"].push_back(name);
						continue;
					}
					const auto file = fetch(bucket + name + "/" + name + ".tif", directory, options, cancel);
					std::vector<double> lats, lons;
					for (auto i : rows)
						lats.push_back(result.latitudesDeg[i]);
					for (auto j : cols)
						lons.push_back(result.longitudesDeg[j]);
					const auto values = rasterSamples(file, lats, lons, true);
					for (std::size_t i = 0; i < rows.size(); ++i)
						for (std::size_t j = 0; j < cols.size(); ++j) {
							const auto h = values[i * cols.size() + j];
							if (h < -500 || h > 9000)
								throw std::runtime_error("Invalid DEM elevation");
							result.heightsWgs84M[rows[i] * result.longitudesDeg.size() + cols[j]] = h;
						}
					provenance["raw_sources"].push_back(file.filename().string());
				}
			const auto undulations = rasterSamples(geoid, result.latitudesDeg, result.longitudesDeg);
			for (std::size_t k = 0; k < count; ++k)
				result.heightsWgs84M[k] += undulations[k];
			return result;
		}

	}  // namespace

	void validateAutomatic(const AutomaticEnvironment& options)
	{
		if (!std::isfinite(options.marginKm) || options.marginKm < 1 || options.marginKm > 500 ||
		    !std::isfinite(options.terrainStepDeg) || options.terrainStepDeg < 0.001 || options.terrainStepDeg > 0.1 ||
		    options.timeoutSeconds < 1 || options.timeoutSeconds > 3600 || options.maxDownloadBytes == 0 ||
		    options.maxDownloadBytes > 1024ULL * 1024 * 1024)
			throw std::invalid_argument("Invalid automatic environment bounds, resolution or download limits");
		if (options.boundsDeg)
			validateBounds(*options.boundsDeg);
		if (options.cycleMjdUtc && !options.gfsFile.empty())
			throw std::invalid_argument("gfs_file and cycle_mjd_utc are mutually exclusive");
	}

	Bounds environmentBounds(const AutomaticEnvironment& options, std::span<const EntryState> entries)
	{
		validateAutomatic(options);
		if (entries.empty())
			throw std::invalid_argument("Cannot prepare environment without entries");
		Bounds bounds{90, -90, 180, -180};
		const auto add = [&](const Vector3& point) {
			const auto pos = geodetic(point);
			const double lat = pos.latitudeRad / detail::deg, lon = pos.longitudeRad / detail::deg;
			bounds[0] = std::min(bounds[0], lat);
			bounds[1] = std::max(bounds[1], lat);
			bounds[2] = std::min(bounds[2], lon);
			bounds[3] = std::max(bounds[3], lon);
		};
		for (const auto& entry : entries) {
			if (!detail::finite(entry.velocityEcefMps) || !std::isfinite(entry.epochMjdUtc) ||
			    !std::isfinite(entry.weight) || entry.weight < 0)
				throw std::invalid_argument("Invalid entry for automatic environment");
			add(entry.positionEcefM);
			if (options.boundsDeg)
				continue;
			const auto& p = entry.positionEcefM;
			const auto& v = entry.velocityEcefMps;
			double a = 0, b = 0, c = -1;
			for (std::size_t k = 0; k < 3; ++k) {
				const double radius = k == 2 ? detail::earthA * (1 - detail::earthF) : detail::earthA;
				const double x = p[k] / radius, u = v[k] / radius;
				a += u * u;
				b += 2 * x * u;
				c += x * x;
			}
			const double discriminant = b * b - 4 * a * c;
			if (a <= 0 || b >= 0 || discriminant < 0)
				throw std::invalid_argument(
				    "Entry has no downward straight-line ground intersection; supply bounds_deg");
			const double time = std::max(0.0, (-b - std::sqrt(discriminant)) / (2 * a));
			Vector3 point;
			for (std::size_t k = 0; k < 3; ++k)
				point[k] = p[k] + v[k] * time;
			add(point);
		}
		if (options.boundsDeg) {
			const auto& requested = *options.boundsDeg;
			if (bounds[0] < requested[0] || bounds[1] > requested[1] || bounds[2] < requested[2] ||
			    bounds[3] > requested[3])
				throw std::invalid_argument("Configured environment bounds exclude entry states");
			return requested;
		}
		const double latPad = options.marginKm / 111.195;
		const double lonPad = latPad / std::cos(std::max(std::abs(bounds[0]), std::abs(bounds[1])) * detail::deg);
		bounds[0] -= latPad;
		bounds[1] += latPad;
		bounds[2] -= lonPad;
		bounds[3] += lonPad;
		validateBounds(bounds);
		return bounds;
	}

	Json automaticEnvironment(const AutomaticEnvironment& options, std::span<const EntryState> entries,
	                          const std::filesystem::path& directory, std::stop_token cancel)
	{
		const auto bounds = environmentBounds(options, entries);
		const auto [first, last] = std::minmax_element(
		    entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.epochMjdUtc < b.epochMjdUtc; });
		if (last->epochMjdUtc - first->epochMjdUtc > 1.0 / 24)
			throw std::invalid_argument("Entry times span over one hour; split events into separate weather jobs");
		const double event = (first->epochMjdUtc + last->epochMjdUtc) / 2;
		try {
			Json weatherInfo, terrainInfo;
			const auto weather = obtainGfs(options, event, bounds, directory / "weather", cancel, weatherInfo);
			const std::string geoidUrl = "https://cdn.proj.org/us_nga_egm08_25.tif";
			const auto geoid = fetch(geoidUrl, directory / "terrain", options, cancel);
			const auto undulations = rasterSamples(geoid, weather.latitudes, weather.longitudes);
			const auto atmosphere = atmosphereFromGfs(weather, undulations);
			auto ground = terrain(options, bounds, directory / "terrain", geoid, cancel, terrainInfo);
			auto snapshot = jobfiles::environmentJson(Environment(atmosphere, std::move(ground)));
			snapshot["source"] = "automatic NOAA GFS + Copernicus DEM + EGM2008";
			snapshot["idealized_fallback"] = false;
			snapshot["preparation"] = {
			    {"bounds_deg", bounds},
			    {"entry_epoch_range_mjd_utc", {first->epochMjdUtc, last->epochMjdUtc}},
			    {"weather", weatherInfo},
			    {"terrain", terrainInfo},
			    {"geoid", {{"model", "EGM2008"}, {"source_url", geoidUrl}, {"raw_source", geoid.filename().string()}}}};
			return snapshot;
		} catch (const Unavailable& error) {
			if (!options.idealizedFallback)
				throw;
			auto snapshot = jobfiles::environmentJson(Environment::idealized(0));
			snapshot["idealized_fallback"] = true;
			snapshot["source"] = "US76 / zero wind / flat WGS84 surface; automatic sources unavailable";
			snapshot["fallback_reason"] = error.what();
			return snapshot;
		}
	}

}  // namespace strewn::sources
