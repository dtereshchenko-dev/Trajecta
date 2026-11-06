#include "orbital_core/obs/ObservatoryUpdates.h"
#include "obs/ObservatoryCatalog.h"
#include <gtest/gtest.h>
#include <iomanip>
#include <sstream>

using namespace orbital_core::obs;

namespace {
	std::string opticalResponse()
	{
		std::ostringstream out;
		out << "<pre>\nCode  Long.   cos      sin    Name\n";
		for (int i = 0; i < 1000; ++i) {
			out << std::setfill('0') << std::setw(3) << i << std::setfill(' ');
			if (i == 258)
				out << std::string(27, ' ') << "Gaia\n";
			else
				out << ' ' << std::fixed << std::setprecision(4) << std::setw(9) << 12.0
				    << std::setprecision(6) << std::setw(8) << 0.7 << std::showpos << std::setw(9) << 0.6
				    << std::noshowpos << "Station " << i << '\n';
		}
		out << "</pre>\n";
		return out.str();
	}

	std::string radarResponse()
	{
		std::ostringstream out;
		out << R"({"signature":{"version":"1.1"},"coords":{"-1":{"longitude":"293.25","d_xy":404849.13,"d_z":133339.19,"name":"Arecibo"})";
		for (int i = 990; i < 999; ++i)
			out << ",\"-" << i << R"(":{"longitude":"10","d_xy":1000,"d_z":1000,"name":"New radar"})";
		out << "}}";
		return out.str();
	}
}

TEST(ObservatoryUpdates, OpticalPreservesHistoricalCodesAndSpacecraft)
{
	const auto result = prepareOpticalObservatoryUpdate(opticalResponse(), "Code Long cos sin Name\nZ99 1 0.7 0.6 Old station\n");
	EXPECT_EQ(result.updated, 1000);
	EXPECT_NE(result.contents.find("Z99 1 0.7 0.6 Old station"), std::string::npos);
	EXPECT_NE(result.contents.find("258 Gaia"), std::string::npos);
	EXPECT_NE(result.contents.find("500 12.0000 0.700000 +0.600000 Station 500"), std::string::npos);
}

TEST(ObservatoryUpdates, RejectsIncompleteOpticalResponse)
{
	EXPECT_THROW((void)prepareOpticalObservatoryUpdate("<html>error</html>", ""), std::runtime_error);
	EXPECT_THROW((void)prepareOpticalObservatoryUpdate("<pre>\n500 0 0 0\n</pre>", ""), std::runtime_error);
	EXPECT_THROW(validateOpticalObservatoryCatalog("Code Long. cos sin Name\n500 0 0 0 Geocentric\n"), std::runtime_error);
}

TEST(ObservatoryUpdates, UnchangedOpticalDataIgnoreFormatting)
{
	const auto first = prepareOpticalObservatoryUpdate(opticalResponse(), "");
	EXPECT_TRUE(first.changed);
	const auto second = prepareOpticalObservatoryUpdate(opticalResponse(), first.contents);
	EXPECT_FALSE(second.changed);
	std::string formatted = first.contents;
	const auto header = formatted.find("Code Long. cos sin Name");
	ASSERT_NE(header, std::string::npos);
	formatted.replace(header, 23, "Code  Long.  cos  sin  Name");
	EXPECT_FALSE(prepareOpticalObservatoryUpdate(opticalResponse(), formatted).changed);
	EXPECT_TRUE(prepareOpticalObservatoryUpdate(opticalResponse(), formatted + "500 12 0.7 0.6 Station 500\n").changed);
	auto edited = opticalResponse();
	const auto name = edited.find("Station 999");
	ASSERT_NE(name, std::string::npos);
	edited.replace(name, 11, "New station name");
	EXPECT_TRUE(prepareOpticalObservatoryUpdate(edited, formatted).changed);
}

TEST(ObservatoryUpdates, RadarUpdatesKnownStationAndPreservesAliases)
{
	const std::string existing = "! RADCODE format=2\n251 -1 251 293.0 400000 130000 Arecibo\n"
	                             "257 -25 257 243.0 348000 245000 DSS25\n";
	const auto result = prepareRadarObservatoryUpdate(radarResponse(), existing);
	EXPECT_EQ(result.updated, 10);
	EXPECT_EQ(result.unmapped, 9);
	std::istringstream input(result.contents);
	const auto catalog = readRadarCatalog(input);
	ASSERT_EQ(catalog.stations.size(), 11);
	EXPECT_TRUE(catalog.stations.contains("-990"));
	EXPECT_EQ(catalog.aliases.at("-1"), "251");
	EXPECT_EQ(catalog.aliases.at("-25"), "257");
	EXPECT_DOUBLE_EQ(catalog.stations.at("251").longitude, 293.25);
	EXPECT_DOUBLE_EQ(catalog.stations.at("257").longitude, 243.0);
}

TEST(ObservatoryUpdates, RadarCanBootstrapWithoutExistingCatalog)
{
	const auto result = prepareRadarObservatoryUpdate(radarResponse(), "");
	EXPECT_EQ(result.updated, 10);
	EXPECT_EQ(result.unmapped, 10);
	std::istringstream input(result.contents);
	const auto catalog = readRadarCatalog(input);
	EXPECT_EQ(catalog.stations.size(), 10);
	EXPECT_TRUE(catalog.stations.contains("-1"));
}

TEST(ObservatoryUpdates, UnchangedRadarDataDoNotDependOnFormatting)
{
	const auto first = prepareRadarObservatoryUpdate(radarResponse(), "");
	EXPECT_TRUE(first.changed);
	const auto second = prepareRadarObservatoryUpdate(radarResponse(), first.contents);
	EXPECT_FALSE(second.changed);
	auto edited = radarResponse();
	const auto longitude = edited.find("293.25");
	ASSERT_NE(longitude, std::string::npos);
	edited.replace(longitude, 6, "293.26");
	EXPECT_TRUE(prepareRadarObservatoryUpdate(edited, first.contents).changed);
}

TEST(ObservatoryUpdates, RejectsInvalidOrUnmatchedRadarResponse)
{
	const std::string existing = "! RADCODE format=2\n251 -1 251 293 400000 130000 Arecibo\n";
	EXPECT_THROW((void)prepareRadarObservatoryUpdate("{}", existing), std::runtime_error);
	EXPECT_THROW((void)prepareRadarObservatoryUpdate(R"({"signature":{"version":"1.1"},"coords":{"-2":{}}})", existing),
	             std::runtime_error);
}

TEST(ObservatoryUpdates, CancelledDownloadsDoNotStart)
{
	std::stop_source cancellation;
	cancellation.request_stop();
	EXPECT_THROW((void)downloadOpticalObservatoryUpdate("", cancellation.get_token()), std::runtime_error);
	EXPECT_THROW((void)downloadRadarObservatoryUpdate("", cancellation.get_token()), std::runtime_error);
}
