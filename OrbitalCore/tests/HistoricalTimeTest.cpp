#include <gtest/gtest.h>

#include "earth_orientation/EarthOrientation.h"
#include "AdesParser.h"
#include "obs/Mpc80ColumnReader.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/ObservationSet.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/Types.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>

namespace {
	using namespace orbital_core;

	double mjd(int y, unsigned m, unsigned d)
	{
		using namespace std::chrono;
		return static_cast<double>(sys_days(year{y} / month{m} / day{d}).time_since_epoch().count()) + 40587;
	}

	std::filesystem::path dataDirectory()
	{
		return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	}

	class HistoricalTimeTest : public ::testing::Test {
	   protected:
		std::filesystem::path directory;

		void SetUp() override
		{
			directory =
			    std::filesystem::temp_directory_path() /
			    ("orbital-deltat-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
			ASSERT_TRUE(std::filesystem::create_directories(directory / "earth_orientation"));
		}

		void TearDown() override
		{
			std::filesystem::remove(directory / "earth_orientation/ET-UT.dat");
			std::filesystem::remove(directory / "earth_orientation");
			std::filesystem::remove(directory);
		}

		void write(const std::string& records)
		{
			std::ofstream file(directory / "earth_orientation/ET-UT.dat");
			ASSERT_TRUE(file.is_open());
			file << "TT - UT1 seconds\n--------------------\n" << records;
		}
	};

	TEST_F(HistoricalTimeTest, InterpolatesUt1AndInvertsAcrossNegativeDeltaT)
	{
		write("1 1 1900 -2.0\n1 1 1901 3.0\n");
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		const double start = mjd(1900, 1, 1), end = mjd(1901, 1, 1);
		for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
			const double ut1 = start + fraction * (end - start);
			double seconds = 0, recovered = 0;
			ASSERT_TRUE(table.ut1ToTtSeconds(ut1, seconds));
			EXPECT_DOUBLE_EQ(seconds, -2.0 + 5.0 * fraction);
			ASSERT_TRUE(table.ttToUt1Mjd(ut1 + seconds / 86400., recovered));
			EXPECT_NEAR(recovered, ut1, 2e-12);
		}
	}

	TEST_F(HistoricalTimeTest, MissingAndOutOfRangeTablesNeverInventAnOffset)
	{
		double value = 123.0;
		const EarthOrientation missing(directory, EarthOrientationProfile::IAU2000);
		EXPECT_FALSE(missing.ut1ToTtSeconds(mjd(1930, 1, 1), value));
		EXPECT_DOUBLE_EQ(value, 123.0);
		write("1 1 1930 24.02\n1 1 1931 24.00\n");
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		for (double t : {mjd(1929, 12, 31), mjd(1931, 1, 2), std::numeric_limits<double>::quiet_NaN()}) {
			EXPECT_FALSE(table.ut1ToTtSeconds(t, value));
			EXPECT_FALSE(table.ttToUt1Mjd(t, value));
			EXPECT_DOUBLE_EQ(value, 123.0);
		}
	}

	TEST_F(HistoricalTimeTest, RejectsMalformedDuplicateAndReversedDates)
	{
		for (const std::string rows : {"1 1 1930 24.02\n1 1 1930 24.01\n", "1 1 1931 24.00\n1 1 1930 24.02\n",
		                               "1 1 1930 24.02\n31 2 1931 24.00\n", "1 1 1930 24.02\n1 1 1931 24.00 extra\n",
		                               "1 1 1930 24.02\ninvalid\n"}) {
			write(rows);
			const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
			double seconds = 123.0;
			EXPECT_FALSE(table.ut1ToTtSeconds(mjd(1930, 1, 1), seconds)) << rows;
			EXPECT_DOUBLE_EQ(seconds, 123.0);
		}
	}

	TEST_F(HistoricalTimeTest, PackagedTableCoversApolloAndBridgesToDailyFinals)
	{
		const EarthOrientation table(dataDirectory(), EarthOrientationProfile::IAU2000);
		const double epoch = 2426323.692360 - time::kJdAtMjdZero;
		const double fraction = (epoch - mjd(1930, 1, 1)) / 365.;
		double seconds = 0, recovered = 0;
		ASSERT_TRUE(table.ut1ToTtSeconds(epoch, seconds));
		EXPECT_NEAR(seconds, 24.02 + fraction * (24.00 - 24.02), 1e-13);
		ASSERT_TRUE(table.ttToUt1Mjd(epoch + seconds / 86400., recovered));
		EXPECT_NEAR(recovered, epoch, 4e-12);
		for (const double t : {mjd(1800, 1, 1), mjd(1932, 5, 1), mjd(1972, 12, 31), mjd(1973, 1, 1), mjd(1973, 1, 2),
		                       mjd(2020, 11, 16)}) {
			ASSERT_TRUE(table.ut1ToTtSeconds(t, seconds));
			EXPECT_TRUE(std::isfinite(seconds));
		}
		EXPECT_FALSE(table.ut1ToTtSeconds(mjd(1799, 12, 31), seconds));
		EXPECT_FALSE(table.ut1ToTtSeconds(mjd(2200, 1, 1), seconds));
	}

	TEST_F(HistoricalTimeTest, UtcDriftIsDistinctFromHistoricalUt1)
	{
		EXPECT_NEAR(Convertors::utcToTdtSeconds(time::jdFromMjd(37300)), 33.606818, 1e-12);
		EXPECT_NEAR(Convertors::utcToTdtSeconds(time::jdFromMjd(38942)),
		            32.184 + 3.7401300 + (38942 - 38761) * .0012960, 1e-12);
		EXPECT_NEAR(Convertors::utcToTdtSeconds(time::jdFromMjd(39887)),
		            32.184 + 4.2131700 + (39887 - 39126) * .0025920, 1e-12);
		EXPECT_DOUBLE_EQ(Convertors::utcToTdtSeconds(2441317.5), 42.184);
		EXPECT_DOUBLE_EQ(Convertors::utcToTdtSeconds(2459169.5), 69.184);
		EXPECT_TRUE(std::isnan(Convertors::utcToTdtSeconds(2426323.692360)));
		for (const double t : {37300.25, 38942.25, 39887.25, 59169.25}) {
			const auto tt = time::toMjdTt(t, time::Scale::UTC);
			EXPECT_NEAR(time::fromMjdTt(tt, time::Scale::UTC), t, 1e-11);
		}
	}

	TEST_F(HistoricalTimeTest, FormatDeterminesHistoricalTimeScale)
	{
		std::string error;
		const auto old = obs::Mpc80ColumnReader{}.read(
		    "01862          1930 12 13.19236 02 45 33.68 +30 58 37.7                 19628690\n", error);
		ASSERT_EQ(old.size(), 1u) << error;
		EXPECT_EQ(old[0].timeScale, time::Scale::UT1);
		EXPECT_DOUBLE_EQ(old[0].jd, 2426323.692360);
		const auto ades = ades::AdesReader{}.read(
		    "<ades version=\"2022\"><optical><permID>1862</permID><stn>690</stn>"
		    "<mode>CCD</mode><obsTime>1965-07-01T00:00:00Z</obsTime>"
		    "<ra>40.0</ra><dec>30.0</dec></optical></ades>",
		    error);
		ASSERT_EQ(ades.size(), 1u) << error;
		EXPECT_EQ(ades[0].timeScale, time::Scale::UTC);
		EXPECT_EQ(time::legacyObservationScale(2441317.5), time::Scale::UTC);
	}

	TEST_F(HistoricalTimeTest, LegacyAdesXmlAndPsvRetainTheSubmittedTimeScale)
	{
		for (const std::string epoch :
		     {"1951-09-30T10:55:32.736Z", "1960-01-01T00:00:00Z", "1965-07-01T00:00:00Z", "1971-12-31T23:59:59.999Z",
		      "1972-01-01T00:00:00Z", "2020-11-12T04:38:36.067Z"}) {
			for (const std::string format : {"PRE", "M47", "M92", "A17", "A22", "", "Z99"}) {
				const std::string xml =
				    "<ades version=\"2022\"><optical><permID>24P</permID>"
				    "<stn>672</stn><mode>PHO</mode><obsTime>" +
				    epoch + "</obsTime><ra>79.117</ra><dec>13.12422</dec>" +
				    (format.empty() ? "" : "<subFmt>" + format + "</subFmt>") + "</optical></ades>";
				const std::string psv =
				    "permID|stn|mode|obsTime|ra|dec|subFmt\n24P|672|PHO|" + epoch + "|79.117|13.12422|" + format + "\n";
				for (const auto& text : {xml, psv}) {
					SCOPED_TRACE(epoch + " " + format);
					std::string error;
					const auto rows = ades::AdesReader{}.read(text, error);
					ASSERT_EQ(rows.size(), 1u) << error;
					const auto& row = rows.front();
					const bool legacy = format == "PRE" || format == "M47" || format == "M92";
					EXPECT_EQ(row.timeScale, legacy ? time::legacyObservationScale(row.jd) : time::Scale::UTC);
					EXPECT_TRUE(row.hasSplitUtc);
					EXPECT_EQ(row.sourceFormat, format.empty() ? "ADES" : "ADES " + format);
					if (epoch.starts_with("1951")) {
						EXPECT_EQ(row.utcMjdDay, 33919);
						EXPECT_DOUBLE_EQ(row.utcSecondsOfDay, 39332.736);
						if (!legacy)
							EXPECT_TRUE(std::isnan(time::ttOffsetSeconds(time::mjdFromJd(row.jd), row.timeScale)));
					}
				}
			}
		}
	}

	TEST_F(HistoricalTimeTest, SchaumasseLegacyAdesAndMpcHaveTheSamePhysicalEpoch)
	{
		getExecutableDir() = dataDirectory();
		DE440Source source;
		ASSERT_TRUE(source.isAvailable());
		setDefaultEphemerisSource(&source);

		struct Reset {
			~Reset() { setDefaultEphemerisSource(nullptr); }
		} reset;

		std::string error;
		const auto mpc = obs::Mpc80ColumnReader{}.read(
		    "0024P         A1951 09 30.45524 05 16 28.08 +13 07 27.2                 H1148672\n", error);
		ASSERT_EQ(mpc.size(), 1u) << error;
		const auto xml = ades::AdesReader{}.read(
		    "<ades version=\"2022\"><optical><permID>24P</permID><mode>UNK</mode><stn>672</stn>"
		    "<obsTime>1951-09-30T10:55:32.736Z</obsTime><ra>79.11700</ra><dec>13.12422</dec>"
		    "<subFrm>B1950.0</subFrm><subFmt>M92</subFmt><precTime>10</precTime>"
		    "</optical></ades>",
		    error);
		ASSERT_EQ(xml.size(), 1u) << error;
		ASSERT_EQ(xml.front().timeScale, time::Scale::UT1);
		const ObservatoryRecord station{.code = "672", .longitude = 10., .cosPhi = .8, .sinPhi = .6};
		const ObsNode legacy(mpc.front(), station);
		const ObsNode converted(xml.front(), station);
		EXPECT_DOUBLE_EQ(converted.JD, legacy.JD);
		EXPECT_NEAR(converted.ttOffsetSeconds(), legacy.ttOffsetSeconds(), 1e-10);
		EXPECT_NEAR(converted.LST, legacy.LST, 1e-9);
		EXPECT_NEAR(converted.obsVector.x, legacy.obsVector.x, 1e-12);
		EXPECT_NEAR(converted.obsVector.y, legacy.obsVector.y, 1e-12);
		EXPECT_NEAR(converted.obsVector.z, legacy.obsVector.z, 1e-12);
		EXPECT_GT(converted.ttOffsetSeconds(), 29.);
		EXPECT_LT(converted.ttOffsetSeconds(), 30.);
	}

	TEST_F(HistoricalTimeTest, ApolloLoadsAllObservationsWithHistoricalClock)
	{
		getExecutableDir() = dataDirectory();
		DE440Source source;
		ASSERT_TRUE(source.isAvailable());
		setDefaultEphemerisSource(&source);

		struct Reset {
			~Reset() { setDefaultEphemerisSource(nullptr); }
		} reset;

		const auto observations = loadObservations(dataDirectory() / "examples/NEO/apollo_1862.txt");
		ASSERT_EQ(observations.all.size(), 3479u);
		const auto& first = observations.all.front();
		EXPECT_DOUBLE_EQ(first.JD, 2426323.692360);
		EXPECT_EQ(first.timeScale, time::Scale::UT1);
		EXPECT_NEAR(first.ttOffsetSeconds(), 24.00103055561644, 1e-10);
		EXPECT_TRUE(std::isfinite(first.LST));
		for (const auto& obs : observations.all)
			EXPECT_TRUE(std::isfinite(obs.ttOffsetSeconds())) << obs.JD;
	}
}  // namespace
