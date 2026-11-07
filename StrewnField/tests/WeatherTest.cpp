#include "GfsSource.h"
#include "strewn_core/Environment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <gdal_priv.h>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

using namespace strewn;
using namespace strewn::sources;

TEST(WeatherTest, GfsRangeDownloadsAreBoundedAndKeepInventoryOrder)
{
	std::vector<trajecta::http::ByteRange> ranges;
	for (std::uint64_t i = 0; i < 9; ++i) ranges.push_back({i, i});
	std::mutex mutex;
	std::condition_variable ready;
	std::size_t active = 0, peak = 0, started = 0;
	auto files = downloadGfsRanges(ranges, [&](auto range, auto) {
		std::unique_lock lock(mutex);
		++active;
		peak = std::max(peak, active);
		if (++started <= 4) {
			ready.notify_all();
			if (!ready.wait_for(lock, std::chrono::seconds(5), [&] { return started >= 4; }))
				throw std::runtime_error("GFS downloads did not run concurrently");
		}
		--active;
		return std::filesystem::path(std::to_string(range.first));
	}, {});
	EXPECT_EQ(peak, 4u);
	ASSERT_EQ(files.size(), ranges.size());
	for (std::size_t i = 0; i < files.size(); ++i) EXPECT_EQ(files[i].string(), std::to_string(i));
}

TEST(WeatherTest, GfsRangeFailureStopsOtherDownloads)
{
	const std::vector<trajecta::http::ByteRange> ranges{{0, 0}, {1, 1}};
	std::atomic<bool> secondStarted{false}, sawStop{false};
	EXPECT_THROW((void)downloadGfsRanges(ranges, [&](auto range, auto stop) {
		const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
		if (range.first == 0) {
			while (!secondStarted.load() && std::chrono::steady_clock::now() < deadline) std::this_thread::yield();
			throw std::runtime_error("failed range");
		}
		secondStarted = true;
		while (!stop.stop_requested() && std::chrono::steady_clock::now() < deadline) std::this_thread::yield();
		sawStop = stop.stop_requested();
		return std::filesystem::path{};
	}, {}), std::runtime_error);
	EXPECT_TRUE(sawStop);
	std::stop_source cancel;
	cancel.request_stop();
	EXPECT_THROW((void)downloadGfsRanges(ranges, [](auto, auto) { return std::filesystem::path{}; }, cancel.get_token()),
	             trajecta::http::Cancelled);
}

TEST(WeatherTest, SelectsHistoricalCycleAndNearestThreeHourValidTime)
{
	const auto value = selectGfs(59988.124, 62000);
	EXPECT_EQ(value.url, "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.20230213/00/atmos/gfs.t00z.pgrb2.0p25.f003");
	EXPECT_EQ(value.validUnix - value.cycleUnix, 10800);
	EXPECT_TRUE(value.indexed);
}

TEST(WeatherTest, RollsOverMidnightAndSelectsMigratedHistoricalArchive)
{
	EXPECT_NE(selectGfs(59988.99, 62000).url.find("gfs.20230214/00/"), std::string::npos);
	const auto la = selectGfs(58271.75, 62000);
	EXPECT_EQ(
	    la.url,
	    "https://www.ncei.noaa.gov/oa/prod-model/global-forecast-system/access/historical/analysis/201806/20180602/"
	    "gfsanl_4_20180602_1800_000.grb2");
	EXPECT_FALSE(la.indexed);
}

TEST(WeatherTest, ForecastUsesIssuedCycleAndRejectsBeyondHorizon)
{
	const auto value = selectGfs(61000.375, 61000.2);
	EXPECT_LE(value.cycleUnix, static_cast<std::int64_t>((61000.2 - 40587) * 86400));
	EXPECT_EQ(value.validUnix - value.cycleUnix, 9 * 3600);
	EXPECT_THROW(selectGfs(61017, 61000), Unavailable);
	EXPECT_THROW(selectGfs(61001, 61000, 61000.25), std::invalid_argument);
	EXPECT_THROW(selectGfs(61000, 62000, 60999.1), Unavailable);
	EXPECT_THROW(selectGfs(61000, 1e100), std::invalid_argument);
	EXPECT_THROW(selectGfs(61000, 62000, -1e100), std::invalid_argument);
}

TEST(WeatherTest, ExplicitCycleKeepsSameValidTime)
{
	const auto value = selectGfs(59988.124, 62000, 59987.75);
	EXPECT_NE(value.url.find("gfs.20230212/18/"), std::string::npos);
	EXPECT_TRUE(value.url.ends_with("f009"));
}

TEST(WeatherTest, InterpolatesPressureColumnsAndAppliesGeoidWithoutVerticalWind)
{
	WeatherFields fields;
	fields.latitudes = {49.5, 50};
	fields.longitudes = {0, 0.5};
	fields.levels[1000] = {
	    {"HGT", {0, 0, 0, 0}}, {"TMP", {280, 280, 280, 280}}, {"UGRD", {10, 10, 10, 10}}, {"VGRD", {-5, -5, -5, -5}}};
	fields.levels[500] = {{"HGT", {5000, 5000, 5000, 5000}},
	                      {"TMP", {250, 250, 250, 250}},
	                      {"UGRD", {30, 30, 30, 30}},
	                      {"VGRD", {5, 5, 5, 5}}};
	auto atmosphere = atmosphereFromGfs(fields, {50, 50, 50, 50});
	const auto low = atmosphere.samples[21];
	EXPECT_DOUBLE_EQ(low.temperatureK, 280);
	EXPECT_DOUBLE_EQ(low.densityKgM3, 100000 / (287.0531 * 280));
	EXPECT_DOUBLE_EQ(low.windEnuMps[2], 0);
	const auto high = atmosphere.samples[821];
	EXPECT_DOUBLE_EQ(high.temperatureK, standardAtmosphere(40050).temperatureK);
	EXPECT_DOUBLE_EQ(high.windEnuMps[0], 0);
	for (const auto& sample : atmosphere.samples)
		EXPECT_DOUBLE_EQ(sample.windEnuMps[2], 0);
	EXPECT_THROW(atmosphereFromGfs(fields, {}), std::invalid_argument);
	fields.levels[500]["HGT"] = {0, 0, 0, 0};
	EXPECT_THROW(atmosphereFromGfs(fields, {0, 0, 0, 0}), std::runtime_error);
}

TEST(WeatherTest, RequiredRasterDriversAreAvailable)
{
	initializeRasters();
	for (const auto* name : {"GRIB", "GTiff", "JP2OpenJPEG", "CSV", "MEM"}) {
		EXPECT_NE(GetGDALDriverManager()->GetDriverByName(name), nullptr) << name;
	}
}

TEST(WeatherTest, RuntimeRasterResourcesAreAvailable)
{
	initializeRasters();
	EXPECT_NE(CPLFindFile("gdal", "grib2_table_4_2_0_0.csv"), nullptr);
	OGRSpatialReference reference;
	EXPECT_EQ(reference.importFromEPSG(4326), OGRERR_NONE);
}

TEST(WeatherTest, RasterSamplingPreservesAxesAndPixelCenterConvention)
{
	GDALAllRegister();
	static std::atomic<unsigned> id{0};
	const auto path = std::filesystem::temp_directory_path() /
	                  ("strewn-raster-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
	                   "-" + std::to_string(id++) + ".tif");

	struct Cleanup {
		std::filesystem::path file;

		~Cleanup()
		{
			std::error_code error;
			std::filesystem::remove(file, error);
		}
	} cleanup{path};

	const auto utf8 = path.u8string();
	auto* dataset = GetGDALDriverManager()->GetDriverByName("GTiff")->Create(
	    reinterpret_cast<const char*>(utf8.c_str()), 4, 3, 1, GDT_Float32, nullptr);
	ASSERT_NE(dataset, nullptr);
	double transform[] = {-1, .5, 0, 51, 0, -.5};
	ASSERT_EQ(dataset->SetGeoTransform(transform), CE_None);
	OGRSpatialReference reference;
	ASSERT_EQ(reference.SetWellKnownGeogCS("WGS84"), OGRERR_NONE);
	ASSERT_EQ(dataset->SetSpatialRef(&reference), CE_None);
	std::array<float, 12> values;
	for (int y = 0; y < 3; ++y)
		for (int x = 0; x < 4; ++x)
			values[y * 4 + x] = static_cast<float>(100 + 5 * (50.75 - y * .5) + 3 * (-.75 + x * .5));
	const auto status =
	    dataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, 4, 3, values.data(), 4, 3, GDT_Float32, 0, 0);
	GDALClose(dataset);
	ASSERT_EQ(status, CE_None);
	const std::vector<double> lats{50, 50.5}, lons{-.5, 0, .5};
	const auto sampled = rasterSamples(path, lats, lons);
	for (std::size_t i = 0; i < 2; ++i)
		for (std::size_t j = 0; j < 3; ++j)
			EXPECT_NEAR(sampled[i * 3 + j], 100 + 5 * lats[i] + 3 * lons[j], 1e-10);
	EXPECT_THROW(rasterSamples(path, {45}, {0}), std::runtime_error);
	const auto edge = rasterSamples(path, {49.5}, {.9}, true);
	EXPECT_NEAR(edge[0], 100 + 5 * 49.75 + 3 * .75, 1e-10);
}
