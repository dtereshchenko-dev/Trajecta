#include "FieldOutputInternal.h"
#include "../app/JobFiles.h"
#include "../app/SourceCache.h"
#include "../app/RasterSource.h"
#include <gdal_priv.h>
#include <gdal_alg.h>
#include <ogrsf_frmts.h>
#include <charconv>
#include <cmath>
#include <memory>

namespace strewn::output {
	namespace {
		using Dataset = std::unique_ptr<GDALDataset, decltype(&GDALClose)>;
		using Feature = std::unique_ptr<OGRFeature, decltype(&OGRFeature::DestroyFeature)>;

		double number(OGRFeature& feature, const char* key)
		{
			const int index = feature.GetFieldIndex(key);
			if (index < 0 || !feature.IsFieldSetAndNotNull(index))
				throw std::runtime_error(std::string("Missing impact field: ") + key);
			const std::string value = feature.GetFieldAsString(index);
			double result = 0;
			const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
			if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size() || !std::isfinite(result))
				throw std::runtime_error(std::string("Invalid impact number: ") + key);
			return result;
		}

		int progress(double, const char*, void* data) { return !static_cast<std::stop_token*>(data)->stop_requested(); }

		std::uint64_t fragmentId(OGRFeature& feature)
		{
			const std::string value = feature.GetFieldAsString("fragment_id");
			std::uint64_t result = 0;
			const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
			if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size())
				throw std::runtime_error("Invalid fragment_id in impact archive");
			return result;
		}

	}  // namespace

	ArchivedField readArchive(const std::filesystem::path& directory, std::stop_token cancel)
	{
		sources::initializeRasters();
		ArchivedField result;
		result.summary = jobfiles::readJson(directory / "summary.json");
		result.job = jobfiles::readJson(directory / "job.json");
		const auto manifest = jobfiles::readJson(directory / "manifest.json");
		if (manifest.at("status") != "complete" || result.summary.at("failed") != 0 ||
		    result.summary.at("truncated") != 0 || result.summary.at("cancelled") != false ||
		    result.summary.at("completed") == 0)
			throw std::runtime_error("Map export requires a complete, nontruncated numerical ensemble");
		const auto path = directory / "impacts.csv";
		result.impactFingerprint = sources::fingerprint(sources::readBytes(path, 512 * 1024 * 1024));
		const auto source = "CSV:" + jobfiles::pathText(path);
		const char* drivers[] = {"CSV", nullptr};
		Dataset data(static_cast<GDALDataset*>(
		                 GDALOpenEx(source.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, drivers, nullptr, nullptr)),
		             &GDALClose);
		if (!data || data->GetLayerCount() != 1)
			throw std::runtime_error("Cannot read archived impact table");
		auto* layer = data->GetLayer(0);
		while (Feature f{layer->GetNextFeature(), &OGRFeature::DestroyFeature}) {
			if (cancel.stop_requested())
				throw trajecta::http::Cancelled{};
			if (std::string(f->GetFieldAsString("run_status")) != "complete")
				throw std::runtime_error("Incomplete run in impact table");
			const double id = number(*f, "run_id");
			if (id < 0 || id >= result.summary.at("completed").get<double>() || std::floor(id) != id)
				throw std::runtime_error("Impact run_id outside the accepted prefix");
			result.points.push_back({number(*f, "latitude_deg"), number(*f, "longitude_deg"), number(*f, "mass_kg"),
			                         number(*f, "multiplicity"), number(*f, "run_weight")});
			result.impactDetails.push_back({static_cast<std::uint64_t>(id), fragmentId(*f),
			                                number(*f, "time_since_entry_s"), number(*f, "speed_mps")});
		}
		if (result.points.size() != result.summary.at("tracked_impacts").get<std::size_t>())
			throw std::runtime_error("Impact table size differs from run summary");
		return result;
	}

	Json exportFields(const std::filesystem::path& directory, const DensityFields& grid, std::stop_token cancel)
	{
		sources::initializeRasters();
		auto* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
		if (!driver)
			throw std::runtime_error("GDAL GeoTIFF driver unavailable");
		const char* options[] = {"COMPRESS=DEFLATE", "PREDICTOR=3", nullptr};
		Dataset raster(
		    driver->Create(jobfiles::pathText(directory / "density.tif").c_str(), static_cast<int>(grid.columns),
		                   static_cast<int>(grid.rows), 2, GDT_Float64, const_cast<char**>(options)),
		    &GDALClose);
		if (!raster)
			throw std::runtime_error("Cannot create density GeoTIFF");
		double transform[] = {grid.firstLongitudeDeg - grid.longitudeStepDeg / 2,
		                      grid.longitudeStepDeg,
		                      0,
		                      grid.firstLatitudeDeg + (grid.rows - .5) * grid.latitudeStepDeg,
		                      0,
		                      -grid.latitudeStepDeg};
		OGRSpatialReference reference;
		if (reference.SetWellKnownGeogCS("WGS84") != OGRERR_NONE || raster->SetGeoTransform(transform) != CE_None ||
		    raster->SetSpatialRef(&reference) != CE_None)
			throw std::runtime_error("Cannot set field spatial reference");
		raster->SetMetadataItem(
		    "NORMALIZATION",
		    "each band sums to 1; conditional weighted fragment distribution, not at-least-one impact probability");
		Json geo = {{"type", "FeatureCollection"}, {"features", Json::array()}};
		Json metadata = {
		    {"schema_version", 1},
		    {"frame", "WGS84 longitude/latitude"},
		    {"cell_centers", true},
		    {"row_order", "north_to_south in GeoTIFF; south_to_north in API"},
		    {"rows", grid.rows},
		    {"columns", grid.columns},
		    {"cell_m", grid.cellM},
		    {"first_latitude_deg", grid.firstLatitudeDeg},
		    {"first_longitude_deg", grid.firstLongitudeDeg},
		    {"latitude_step_deg", grid.latitudeStepDeg},
		    {"longitude_step_deg", grid.longitudeStepDeg},
		    {"local_metric_radius_m", 6367000},
		    {"conditioning",
		     "provided atmospheric entries; normalized number/mass distribution, not Earth-impact probability"}};
		for (int channel = 0; channel < 2; ++channel) {
			if (cancel.stop_requested())
				throw trajecta::http::Cancelled{};
			const auto& field = channel == 0 ? grid.count : grid.mass;
			const std::string name = channel == 0 ? "count" : "mass";
			auto* band = raster->GetRasterBand(channel + 1);
			band->SetDescription(name.c_str());
			band->SetUnitType("normalized weight per cell");
			std::vector<double> northFirst(field.probability.size());
			for (std::size_t y = 0; y < grid.rows; ++y)
				std::copy_n(field.probability.data() + y * grid.columns, grid.columns,
				            northFirst.data() + (grid.rows - y - 1) * grid.columns);
			if (band->RasterIO(GF_Write, 0, 0, static_cast<int>(grid.columns), static_cast<int>(grid.rows),
			                   northFirst.data(), static_cast<int>(grid.columns), static_cast<int>(grid.rows),
			                   GDT_Float64, 0, 0) != CE_None)
				throw std::runtime_error("Cannot write density band");
			auto* memoryDriver = GetGDALDriverManager()->GetDriverByName("MEM");
			if (!memoryDriver || !memoryDriver->GetMetadataItem(GDAL_DCAP_VECTOR))
				memoryDriver = GetGDALDriverManager()->GetDriverByName("Memory");
			if (!memoryDriver)
				throw std::runtime_error("GDAL Memory vector driver unavailable");
			Dataset contours(memoryDriver->Create("", 0, 0, 0, GDT_Unknown, nullptr), &GDALClose);
			if (!contours)
				throw std::runtime_error("Cannot allocate contour layer");
			auto* layer = contours->CreateLayer("contours", &reference, wkbLineString);
			if (!layer)
				throw std::runtime_error("Cannot create contour layer");
			OGRFieldDefn level("threshold", OFTReal);
			if (layer->CreateField(&level) != OGRERR_NONE)
				throw std::runtime_error("Cannot create contour attribute");
			auto thresholds = field.hpdThresholds;
			auto ordered = thresholds;
			const double contourScale = *std::min_element(thresholds.begin(), thresholds.end());
			for (auto& value : northFirst)
				value /= contourScale;
			for (auto& value : ordered)
				value /= contourScale;
			std::sort(ordered.begin(), ordered.end());
			const int levelCount = static_cast<int>(std::unique(ordered.begin(), ordered.end()) - ordered.begin());
			auto* rasterMemory = GetGDALDriverManager()->GetDriverByName("MEM");
			if (!rasterMemory)
				throw std::runtime_error("GDAL MEM raster driver unavailable");
			Dataset normalized(rasterMemory->Create("", static_cast<int>(grid.columns), static_cast<int>(grid.rows), 1,
			                                        GDT_Float64, nullptr),
			                   &GDALClose);
			if (!normalized || normalized->SetGeoTransform(transform) != CE_None)
				throw std::runtime_error("Cannot allocate normalized contour raster");
			auto* contourBand = normalized->GetRasterBand(1);
			if (contourBand->RasterIO(GF_Write, 0, 0, static_cast<int>(grid.columns), static_cast<int>(grid.rows),
			                          northFirst.data(), static_cast<int>(grid.columns), static_cast<int>(grid.rows),
			                          GDT_Float64, 0, 0) != CE_None)
				throw std::runtime_error("Cannot prepare normalized contour raster");
			if (GDALContourGenerate(contourBand, 0, 0, levelCount, ordered.data(), false, 0, layer, -1, 0, progress,
			                        &cancel) != CE_None)
				throw std::runtime_error("Contour extraction failed or cancelled");
			const std::array<double, 3> hpd{.5, .8, .99};
			while (Feature f{layer->GetNextFeature(), &OGRFeature::DestroyFeature}) {
				const auto* line = f->GetGeometryRef()->toLineString();
				const double value = f->GetFieldAsDouble(0) * contourScale;
				std::size_t index = 0;
				for (std::size_t j = 1; j < 3; ++j)
					if (std::abs(value - thresholds[j]) < std::abs(value - thresholds[index]))
						index = j;
				Json coordinates = Json::array();
				for (int i = 0; i < line->getNumPoints(); ++i)
					coordinates.push_back({line->getX(i), line->getY(i)});
				geo["features"].push_back(
				    {{"type", "Feature"},
				     {"properties", {{"quantity", name}, {"hpd", hpd[index]}, {"threshold", value}}},
				     {"geometry", {{"type", "LineString"}, {"coordinates", coordinates}}}});
			}
			metadata[name] = {{"total_weight", field.totalWeight},
			                  {"effective_samples", field.effectiveSamples},
			                  {"kernel_covariance_m2", field.covarianceM2},
			                  {"fixed_kernel_fallback", field.fixedKernelFallback},
			                  {"hpd_levels", hpd},
			                  {"hpd_thresholds", thresholds}};
		}
		raster.reset();
		jobfiles::writeJson(directory / "density.json", metadata);
		jobfiles::writeJson(directory / "contours.geojson", geo);
		return geo;
	}

}  // namespace strewn::output
