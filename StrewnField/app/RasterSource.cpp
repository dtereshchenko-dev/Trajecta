#include "RasterSource.h"
#include <gdal_priv.h>
#include <ogr_srs_api.h>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#include <algorithm>
#include <cmath>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace strewn::sources {
	void initializeRasters()
	{
		static std::once_flag initialized;
		std::call_once(initialized, [] {
#ifdef _WIN32
			HMODULE module{};
			std::array<wchar_t, 32768> modulePath{};
			if (GetModuleHandleExW(
			        GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
			        reinterpret_cast<LPCWSTR>(&initializeRasters), &module)) {
				const auto length =
				    GetModuleFileNameW(module, modulePath.data(), static_cast<DWORD>(modulePath.size()));
				if (length > 0 && length < modulePath.size()) {
					const auto directory = std::filesystem::path(modulePath.data()).parent_path() / "share";
					if (!CPLGetConfigOption("GDAL_DATA", nullptr) &&
					    std::filesystem::is_regular_file(directory / "gdal" / "grib2_table_4_2_0_0.csv")) {
						const auto path = (directory / "gdal").u8string();
						CPLSetConfigOption("GDAL_DATA", reinterpret_cast<const char*>(path.c_str()));
					}
					if (!CPLGetConfigOption("PROJ_DATA", nullptr) && !CPLGetConfigOption("PROJ_LIB", nullptr) &&
					    std::filesystem::is_regular_file(directory / "proj" / "proj.db")) {
						const auto path = (directory / "proj").u8string();
						const char* paths[] = {reinterpret_cast<const char*>(path.c_str()), nullptr};
						OSRSetPROJSearchPaths(paths);
					}
				}
			}
#endif
			GDALAllRegister();
		});
	}

	namespace {

		using Dataset = std::unique_ptr<GDALDataset, decltype(&GDALClose)>;

		Dataset open(const std::filesystem::path& path, const char* driver)
		{
			initializeRasters();
			const auto utf8 = path.u8string();
			const char* allowed[] = {driver, nullptr};
			auto* raw = GDALOpenEx(reinterpret_cast<const char*>(utf8.c_str()), GDAL_OF_RASTER | GDAL_OF_READONLY,
			                       allowed, nullptr, nullptr);
			if (!raw)
				throw std::runtime_error("Cannot decode raster source: " + path.string());
			return {static_cast<GDALDataset*>(raw), &GDALClose};
		}

		std::array<double, 6> transform(GDALDataset& dataset)
		{
			std::array<double, 6> result;
			if (dataset.GetGeoTransform(result.data()) != CE_None || result[1] <= 0 || result[5] >= 0 ||
			    result[2] != 0 || result[4] != 0 || !dataset.GetSpatialRef() ||
			    !dataset.GetSpatialRef()->IsGeographic())
				throw std::runtime_error("Expected a north-up geographic raster");
			return result;
		}

		std::vector<double> read(GDALRasterBand& band, int x, int y, int width, int height)
		{
			if (width <= 0 || height <= 0 || x < 0 || y < 0 || x + width > band.GetXSize() ||
			    y + height > band.GetYSize())
				throw std::runtime_error("Requested region is outside the raster");
			std::vector<double> result(static_cast<std::size_t>(width) * height);
			if (band.RasterIO(GF_Read, x, y, width, height, result.data(), width, height, GDT_Float64, 0, 0) != CE_None)
				throw std::runtime_error("Raster decode failed");
			int hasMissing = 0;
			const double missing = band.GetNoDataValue(&hasMissing);
			for (double value : result)
				if (!std::isfinite(value) || (hasMissing && value == missing))
					throw std::runtime_error("Missing or invalid raster values in requested region");
			return result;
		}

		std::string metadata(GDALRasterBand& band, const char* key)
		{
			const auto* value = band.GetMetadataItem(key);
			return value ? value : "";
		}

	}  // namespace

	std::vector<double> rasterSamples(const std::filesystem::path& path, const std::vector<double>& latitudes,
	                                  const std::vector<double>& longitudes, bool clampEdges)
	{
		auto dataset = open(path, "GTiff");
		if (dataset->GetRasterCount() != 1 || latitudes.empty() || longitudes.empty())
			throw std::runtime_error("Expected a single-band geoid raster and nonempty axes");
		const auto t = transform(*dataset);
		const auto coordinate = [clampEdges](double pixel, int count) {
			if (clampEdges) {
				if (pixel < -1.001 || pixel > count + 0.001)
					throw std::runtime_error("Point outside DEM tile");
				return std::clamp(pixel, 0.0, static_cast<double>(count - 1));
			}
			if (pixel < 0 || pixel > count - 1)
				throw std::runtime_error("Point outside geoid raster");
			return pixel;
		};
		const auto pixelX = [&](double lon) {
			return coordinate((lon - t[0]) / t[1] - 0.5, dataset->GetRasterXSize());
		};
		const auto pixelY = [&](double lat) {
			return coordinate((lat - t[3]) / t[5] - 0.5, dataset->GetRasterYSize());
		};
		const int x0 =
		    std::min(static_cast<int>(std::floor(pixelX(longitudes.front()))), dataset->GetRasterXSize() - 2);
		const int y0 = std::min(static_cast<int>(std::floor(pixelY(latitudes.back()))), dataset->GetRasterYSize() - 2);
		const int width =
		    std::min(static_cast<int>(std::floor(pixelX(longitudes.back()))) + 2, dataset->GetRasterXSize()) - x0;
		const int height =
		    std::min(static_cast<int>(std::floor(pixelY(latitudes.front()))) + 2, dataset->GetRasterYSize()) - y0;
		const auto grid = read(*dataset->GetRasterBand(1), x0, y0, width, height);
		std::vector<double> result;
		for (double lat : latitudes)
			for (double lon : longitudes) {
				const double x = pixelX(lon) - x0, y = pixelY(lat) - y0;
				const int i = std::min(static_cast<int>(std::floor(x)), width - 2),
				          j = std::min(static_cast<int>(std::floor(y)), height - 2);
				const double a = x - i, b = y - j;
				const auto k = static_cast<std::size_t>(j) * width + i;
				result.push_back((1 - b) * ((1 - a) * grid[k] + a * grid[k + 1]) +
				                 b * ((1 - a) * grid[k + width] + a * grid[k + width + 1]));
			}
		return result;
	}

	void readGrib(const std::filesystem::path& path, const Bounds& bounds, WeatherFields& result)
	{
		auto dataset = open(path, "GRIB");
		const auto t = transform(*dataset);
		const double west = t[0] + t[1] / 2, north = t[3] + t[5] / 2;
		const int x0 = static_cast<int>(std::floor((bounds[2] - west) / t[1] + 1e-8));
		const int x1 = static_cast<int>(std::ceil((bounds[3] - west) / t[1] - 1e-8));
		const int y0 = static_cast<int>(std::floor((bounds[1] - north) / t[5] + 1e-8));
		const int y1 = static_cast<int>(std::ceil((bounds[0] - north) / t[5] - 1e-8));
		const int width = x1 - x0 + 1, height = y1 - y0 + 1;
		if (width < 2 || height < 2 || static_cast<std::int64_t>(width) * height > 512)
			throw std::runtime_error("GFS region must contain 2..512 grid columns");
		std::vector<double> lats, lons;
		for (int y = y1; y >= y0; --y)
			lats.push_back(north + y * t[5]);
		for (int x = x0; x <= x1; ++x)
			lons.push_back(west + x * t[1]);
		if (result.latitudes.empty()) {
			result.latitudes = lats;
			result.longitudes = lons;
		} else if (result.latitudes != lats || result.longitudes != lons)
			throw std::runtime_error("GFS fields have inconsistent grids");
		for (int i = 1; i <= dataset->GetRasterCount(); ++i) {
			auto& band = *dataset->GetRasterBand(i);
			const auto name = metadata(band, "GRIB_ELEMENT");
			if (name != "HGT" && name != "TMP" && name != "UGRD" && name != "VGRD")
				continue;
			const auto shortName = metadata(band, "GRIB_SHORT_NAME");
			if (!shortName.ends_with("-ISBL"))
				continue;
			const double pressure = std::stod(shortName);
			if (pressure < 100 || pressure > 100000)
				continue;
			const int level = static_cast<int>(std::lround(pressure / 100));
			const auto valid = std::stoll(metadata(band, "GRIB_VALID_TIME"));
			const auto cycle = std::stoll(metadata(band, "GRIB_REF_TIME"));
			if (result.validUnix && (result.validUnix != valid || result.cycleUnix != cycle))
				throw std::runtime_error("GFS fields refer to different times");
			result.validUnix = valid;
			result.cycleUnix = cycle;
			const auto unit = metadata(band, "GRIB_UNIT");
			if ((name == "TMP" && unit != "[C]" && unit != "[K]") || (name == "HGT" && unit != "[gpm]") ||
			    ((name == "UGRD" || name == "VGRD") && unit != "[m/s]"))
				throw std::runtime_error("Unexpected GFS field units: " + name + " " + unit);
			auto values = read(band, x0, y0, width, height);
			std::vector<double> ordered;
			ordered.reserve(values.size());
			for (int y = height - 1; y >= 0; --y)
				for (int x = 0; x < width; ++x)
					ordered.push_back(values[static_cast<std::size_t>(y) * width + x] +
					                  (name == "TMP" && unit == "[C]" ? 273.15 : 0));
			auto& fields = result.levels[level];
			if (fields.contains(name) && fields.at(name) != ordered)
				throw std::runtime_error("Conflicting duplicate GFS field");
			fields[name] = std::move(ordered);
		}
	}

}  // namespace strewn::sources
