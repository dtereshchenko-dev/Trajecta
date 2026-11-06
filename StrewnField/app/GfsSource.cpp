#include "GfsSource.h"
#include "strewn_core/Environment.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <format>
#include <mutex>
#include <set>
#include <sstream>
#include <thread>

namespace strewn::sources {
	namespace {

		std::string date(std::int64_t seconds)
		{
			const auto day =
			    std::chrono::floor<std::chrono::days>(std::chrono::sys_seconds{std::chrono::seconds{seconds}});
			const std::chrono::year_month_day ymd{day};
			return std::format("{:04}{:02}{:02}", static_cast<int>(ymd.year()), static_cast<unsigned>(ymd.month()),
			                   static_cast<unsigned>(ymd.day()));
		}

		GfsSelection selection(std::int64_t cycle, std::int64_t valid)
		{
			const auto day = date(cycle);
			const int hour = static_cast<int>((cycle % 86400) / 3600);
			const auto step = (valid - cycle) / 3600;
			if (step < 0 || step > 384 || cycle % 21600 != 0)
				throw Unavailable("GFS cycle must be 00/06/12/18 UTC and forecast lead within 0..384 h");
			if (day < "20210226") {
				const auto url =
				    "https://www.ncei.noaa.gov/oa/prod-model/global-forecast-system/access/historical/analysis/" +
				    day.substr(0, 6) + "/" + day + std::format("/gfsanl_4_{}_{:02}00_{:03}.grb2", day, hour, step);
				return {cycle, valid, url, false};
			}
			const auto url = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs." + day +
			                 std::format("/{:02}/atmos/gfs.t{:02}z.pgrb2.0p25.f{:03}", hour, hour, step);
			return {cycle, valid, url, true};
		}

		struct InventoryRow {
			std::uint64_t offset;
			std::string field;
			std::string level;
		};

		std::vector<trajecta::http::ByteRange> ranges(std::string_view text)
		{
			std::vector<InventoryRow> rows;
			std::istringstream lines{std::string(text)};
			for (std::string line; std::getline(lines, line);) {
				std::istringstream parts(line);
				std::vector<std::string> values;
				for (std::string value; std::getline(parts, value, ':');)
					values.push_back(value);
				if (values.size() < 6)
					throw std::runtime_error("Malformed GFS inventory");
				rows.push_back({std::stoull(values[1]), values[3], values[4]});
				if (rows.size() > 20000 || (rows.size() > 1 && rows.back().offset < rows[rows.size() - 2].offset))
					throw std::runtime_error("Invalid GFS inventory offsets");
			}
			std::vector<trajecta::http::ByteRange> result;
			std::set<std::uint64_t> added;
			for (std::size_t i = 0; i < rows.size(); ++i) {
				const auto& row = rows[i];
				if (row.field != "HGT" && row.field != "TMP" && row.field != "UGRD" && row.field != "VGRD")
					continue;
				if (!row.level.ends_with(" mb"))
					continue;
				const double p = std::stod(row.level);
				if (p < 1 || p > 1000)
					continue;
				if (!added.insert(row.offset).second)
					continue;
				auto next = i + 1;
				while (next < rows.size() && rows[next].offset == row.offset)
					++next;
				if (next == rows.size())
					throw std::runtime_error("GFS inventory lacks selected message end");
				const trajecta::http::ByteRange range{row.offset, rows[next].offset - 1};
				if (!result.empty() && result.back().last + 1 == range.first)
					result.back().last = range.last;
				else
					result.push_back(range);
			}
			if (result.empty())
				throw std::runtime_error("No pressure-level weather fields in GFS inventory");
			return result;
		}

	}  // namespace

	std::vector<std::filesystem::path> downloadGfsRanges(
	    std::span<const trajecta::http::ByteRange> ranges,
	    const std::function<std::filesystem::path(trajecta::http::ByteRange, std::stop_token)>& download,
	    std::stop_token cancel)
	{
		constexpr std::size_t maxConcurrentDownloads = 4;
		std::vector<std::filesystem::path> files(ranges.size());
		std::atomic<std::size_t> next{0};
		std::stop_source stop;
		std::stop_callback propagateStop(cancel, [&] { stop.request_stop(); });
		std::mutex mutex;
		std::exception_ptr failure;
		std::vector<std::jthread> workers;
		workers.reserve(std::min(maxConcurrentDownloads, ranges.size()));
		try {
			for (std::size_t i = 0; i < std::min(maxConcurrentDownloads, ranges.size()); ++i) {
				workers.emplace_back([&] {
					try {
						while (!stop.stop_requested()) {
							const auto index = next.fetch_add(1, std::memory_order_relaxed);
							if (index >= ranges.size())
								return;
							files[index] = download(ranges[index], stop.get_token());
						}
					} catch (...) {
						{
							std::lock_guard lock(mutex);
							if (!failure)
								failure = std::current_exception();
						}
						stop.request_stop();
					}
				});
			}
		} catch (...) {
			stop.request_stop();
			throw;
		}
		for (auto& worker : workers)
			worker.join();
		if (failure)
			std::rethrow_exception(failure);
		if (cancel.stop_requested())
			throw trajecta::http::Cancelled{};
		return files;
	}

	GfsSelection selectGfs(double eventMjdUtc, double nowMjdUtc, std::optional<double> cycleMjdUtc)
	{
		const auto supportedDate = [](double mjd) { return std::isfinite(mjd) && mjd >= 53000 && mjd <= 88000; };
		if (!supportedDate(eventMjdUtc) || !supportedDate(nowMjdUtc))
			throw std::invalid_argument("Automatic GFS supports event dates 2004..2099");
		const auto valid = static_cast<std::int64_t>(std::floor((eventMjdUtc - 40587) * 8 + 0.5)) * 10800;
		const auto now = static_cast<std::int64_t>((nowMjdUtc - 40587) * 86400);
		auto cycle = valid / 21600 * 21600;
		if (valid > now - 43200)
			cycle = std::min(cycle, now / 21600 * 21600);
		if (cycleMjdUtc) {
			if (!supportedDate(*cycleMjdUtc))
				throw std::invalid_argument("GFS cycle date is outside the supported range");
			const double seconds = (*cycleMjdUtc - 40587) * 86400;
			cycle = static_cast<std::int64_t>(std::llround(seconds));
			if (std::abs(seconds - static_cast<double>(cycle)) > 1e-4 || cycle > now)
				throw std::invalid_argument("Invalid or future GFS cycle");
		}
		return selection(cycle, valid);
	}

	WeatherFields obtainGfs(const AutomaticEnvironment& options, double eventMjdUtc, const Bounds& bounds,
	                        const std::filesystem::path& directory, std::stop_token cancel, Json& provenance)
	{
		WeatherFields result;
		provenance = {{"provider", "NOAA GFS"}, {"requested_mjd_utc", eventMjdUtc}, {"raw_sources", Json::array()}};
		if (!options.gfsFile.empty()) {
			const auto hash = fingerprint(readBytes(options.gfsFile, options.maxDownloadBytes));
			readGrib(options.gfsFile, bounds, result);
			provenance["local_file"] = options.gfsFile.generic_string();
			provenance["fnv1a64"] = hash;
		} else {
			const double now =
			    40587 +
			    std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count() / 86400;
			auto chosen = selectGfs(eventMjdUtc, now, options.cycleMjdUtc);
			if (chosen.indexed) {
				std::filesystem::path index;
				const int attempts = !options.cycleMjdUtc && eventMjdUtc > now - 0.5 ? 8 : 1;
				for (int i = 0; i < attempts; ++i) {
					try {
						index = fetch(chosen.url + ".idx", directory, options, cancel);
						break;
					} catch (const Unavailable&) {
						if (i + 1 == attempts)
							throw;
						chosen = selection(chosen.cycleUnix - 21600, chosen.validUnix);
					}
				}
				const auto selected = ranges(readBytes(index, options.maxDownloadBytes));
				const auto files = downloadGfsRanges(selected, [&](auto range, auto token) {
					return fetch(chosen.url, directory, options, token, range);
				}, cancel);
				for (const auto& file : files) {
					if (cancel.stop_requested()) throw trajecta::http::Cancelled{};
					readGrib(file, bounds, result);
					provenance["raw_sources"].push_back(file.filename().string());
				}
			} else {
				const auto file = fetch(chosen.url, directory, options, cancel);
				readGrib(file, bounds, result);
				provenance["raw_sources"].push_back(file.filename().string());
			}
			if (result.validUnix != chosen.validUnix || result.cycleUnix != chosen.cycleUnix)
				throw std::runtime_error("Downloaded GFS time does not match selected cycle/lead");
			provenance["source_url"] = chosen.url;
		}
		if (std::abs(result.validUnix - (eventMjdUtc - 40587) * 86400) > 5400.001)
			throw std::runtime_error("GFS valid time is more than 90 minutes from the event");
		for (auto it = result.levels.begin(); it != result.levels.end();) {
			if (it->second.size() != 4)
				it = result.levels.erase(it);
			else
				++it;
		}
		if (result.levels.size() < 20 || result.levels.begin()->first > 10 || result.levels.rbegin()->first < 1000)
			throw std::runtime_error("Incomplete GFS pressure-level column");
		provenance["cycle_mjd_utc"] = 40587 + result.cycleUnix / 86400.0;
		provenance["valid_mjd_utc"] = 40587 + result.validUnix / 86400.0;
		provenance["valid_time_offset_s"] = result.validUnix - (eventMjdUtc - 40587) * 86400;
		provenance["forecast_lead_hours"] = (result.validUnix - result.cycleUnix) / 3600;
		provenance["temporal_interpolation"] = "none; nearest 3-hour valid snapshot";
		provenance["density"] = "dry ideal gas p/(287.0531*T)";
		provenance["vertical_wind"] = "zero; no vertical wind field requested";
		provenance["height_conversion"] = "geopotential -> geometric using R=6356766 m, plus EGM2008 undulation";
		provenance["upper_extension"] = "US76, top ratios and horizontal wind taper to standard over 10 km";
		provenance["lower_extension"] = "constant values below lowest pressure level";
		provenance["levels_hpa"] = Json::array();
		provenance["grid_step_deg"] = {result.latitudes.at(1) - result.latitudes.at(0),
		                               result.longitudes.at(1) - result.longitudes.at(0)};
		for (const auto& [level, fields] : result.levels)
			provenance["levels_hpa"].push_back(level);
		return result;
	}

	AtmosphereGrid atmosphereFromGfs(const WeatherFields& fields, const std::vector<double>& geoidM)
	{
		AtmosphereGrid result;
		result.latitudesDeg = fields.latitudes;
		result.longitudesDeg = fields.longitudes;
		result.firstHeightM = -1000;
		result.heightStepM = 50;
		result.heightCount = 2821;
		const auto columns = fields.latitudes.size() * fields.longitudes.size();
		if (columns == 0 || geoidM.size() != columns)
			throw std::invalid_argument("Geoid/weather dimensions mismatch");
		result.samples.reserve(columns * result.heightCount);

		struct Level {
			double h;
			AirSample air;
		};

		for (std::size_t column = 0; column < columns; ++column) {
			std::vector<Level> profile;
			for (const auto& [p, values] : fields.levels) {
				const double h = values.at("HGT").at(column), t = values.at("TMP").at(column);
				if (!std::isfinite(h) || h < -2000 || h > 100000 || !std::isfinite(t) || t < 100 || t > 400)
					throw std::runtime_error("Invalid GFS pressure-level height or temperature");
				profile.push_back(
				    {6356766 * h / (6356766 - h) + geoidM[column],
				     {p * 100 / (287.0531 * t), t, {values.at("UGRD").at(column), values.at("VGRD").at(column), 0}}});
			}
			std::sort(profile.begin(), profile.end(), [](const auto& a, const auto& b) { return a.h < b.h; });
			if (profile.size() < 2)
				throw std::runtime_error("Incomplete GFS profile");
			for (std::size_t i = 1; i < profile.size(); ++i)
				if (profile[i].h <= profile[i - 1].h)
					throw std::runtime_error("Non-monotonic GFS heights");
			const auto& top = profile.back();
			const auto topStd = standardAtmosphere(top.h);
			std::size_t lo = 0;
			for (std::size_t k = 0; k < result.heightCount; ++k) {
				const double h = result.firstHeightM + k * result.heightStepM;
				AirSample air;
				if (h <= profile.front().h)
					air = profile.front().air;
				else if (h > top.h) {
					air = standardAtmosphere(h);
					const double fade = std::max(0.0, 1 - (h - top.h) / 10000);
					air.densityKgM3 *= 1 + (top.air.densityKgM3 / topStd.densityKgM3 - 1) * fade;
					air.temperatureK *= 1 + (top.air.temperatureK / topStd.temperatureK - 1) * fade;
					air.windEnuMps = {top.air.windEnuMps[0] * fade, top.air.windEnuMps[1] * fade, 0};
				} else {
					while (lo + 2 < profile.size() && profile[lo + 1].h < h)
						++lo;
					const auto& a = profile[lo];
					const auto& b = profile[lo + 1];
					const double f = (h - a.h) / (b.h - a.h);
					air.temperatureK = std::lerp(a.air.temperatureK, b.air.temperatureK, f);
					air.densityKgM3 = std::exp(std::lerp(std::log(a.air.densityKgM3), std::log(b.air.densityKgM3), f));
					for (std::size_t j = 0; j < 3; ++j)
						air.windEnuMps[j] = std::lerp(a.air.windEnuMps[j], b.air.windEnuMps[j], f);
				}
				result.samples.push_back(air);
			}
		}
		return result;
	}

}  // namespace strewn::sources
