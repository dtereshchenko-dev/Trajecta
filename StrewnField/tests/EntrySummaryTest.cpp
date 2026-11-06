#include "../output/EntrySummary.h"
#include "VectorMath.h"
#include <gtest/gtest.h>
#include <limits>

namespace {
	using namespace strewn;
	using output::summarizeEntries;

	EntryState entry(double latitude = 40, double longitude = 10, Vector3 velocity = {1000, 0, -1000},
	                 double weight = 1)
	{
		EntryState e;
		e.positionEcefM = ecefPosition({latitude * detail::deg, longitude * detail::deg, 100000});
		e.velocityEcefMps = enuToEcef(velocity, latitude * detail::deg, longitude * detail::deg);
		e.epochMjdUtc = 59988.;
		e.weight = weight;
		return e;
	}

	TEST(EntrySummary, SingleStateHasNoEstimatedUncertainty)
	{
		const auto s = summarizeEntries(std::vector{entry()});
		EXPECT_NEAR(s["height_wgs84_km"]["mean"].get<double>(), 100., 1e-8);
		EXPECT_NEAR(s["ecef_speed_km_s"]["mean"].get<double>(), std::sqrt(2.), 1e-12);
		EXPECT_NEAR(s["flight_path_deg"]["mean"].get<double>(), -45., 1e-8);
		EXPECT_NEAR(s["heading_deg"]["mean"].get<double>(), 90., 1e-8);
		EXPECT_EQ(s["latitude_deg"]["sigma"], nullptr);
		EXPECT_EQ(s["time_sigma_s"], nullptr);
		EXPECT_EQ(s["ellipse"]["semi_major_m"], nullptr);
	}

	TEST(EntrySummary, WeightedPopulationMomentsAndStableEpoch)
	{
		auto a = entry(40, 10, {0, 1000, 0}, 1), b = entry(40.01, 10, {0, 2000, 0}, 3);
		b.epochMjdUtc += 2. / 86400;
		auto zero = entry();
		zero.weight = 0;
		zero.positionEcefM[0] = std::numeric_limits<double>::quiet_NaN();
		const auto s = summarizeEntries(std::vector{zero, a, b});
		EXPECT_EQ(s["positive_weight_entries"], 2);
		EXPECT_EQ(s["total_entries"], 3);
		EXPECT_NEAR(s["effective_entries"].get<double>(), 1.6, 1e-12);
		EXPECT_NEAR(s["ecef_speed_km_s"]["mean"].get<double>(), 1.75, 1e-12);
		EXPECT_NEAR(s["ecef_speed_km_s"]["sigma"].get<double>(), std::sqrt(3.) / 4, 1e-12);
		EXPECT_NEAR(s["latitude_deg"]["mean"].get<double>(), 40.0075, 1e-8);
		EXPECT_NEAR(s["latitude_deg"]["sigma"].get<double>(), .01 * std::sqrt(3.) / 4, 1e-8);
		EXPECT_NEAR((s["epoch_mjd_utc"].get<double>() - a.epochMjdUtc) * 86400, 1.5, 1e-6);
		EXPECT_NEAR(s["time_sigma_s"].get<double>(), std::sqrt(3.) / 2, 1e-6);
	}

	TEST(EntrySummary, DatelineAndNorthWrapping)
	{
		const double angle = detail::deg;
		auto a = entry(0, 179.999, {-1000 * std::sin(angle), 1000 * std::cos(angle), 0});
		auto b = entry(0, -179.999, {1000 * std::sin(angle), 1000 * std::cos(angle), 0});
		const auto s = summarizeEntries(std::vector{a, b});
		EXPECT_NEAR(std::abs(s["longitude_deg"]["mean"].get<double>()), 180., 1e-9);
		EXPECT_NEAR(s["longitude_deg"]["sigma"].get<double>(), .001, 1e-9);
		EXPECT_NEAR(std::remainder(s["heading_deg"]["mean"].get<double>(), 360.), 0., 1e-9);
		EXPECT_NEAR(s["heading_deg"]["sigma"].get<double>(), 1., 1e-9);
		EXPECT_NEAR(s["ellipse"]["semi_major_m"].get<double>(),
		            (detail::earthA + 100000) * std::sin(.001 * detail::deg), 1e-6);
		EXPECT_NEAR(s["ellipse"]["azimuth_deg"].get<double>(), 90., 1e-9);
	}

	TEST(EntrySummary, EllipseUsesOneSigmaAxesAndNorthAzimuth)
	{
		const auto center = entry();
		std::vector<EntryState> entries;
		const double angle = 35 * detail::deg;
		for (auto offset : std::vector<Vector3>{{60 * std::sqrt(2.), 0, 0},
		                                        {-60 * std::sqrt(2.), 0, 0},
		                                        {0, 10 * std::sqrt(2.), 0},
		                                        {0, -10 * std::sqrt(2.), 0}}) {
			auto e = center;
			const auto d = enuToEcef({offset[0] * std::sin(angle) + offset[1] * std::cos(angle),
			                          offset[0] * std::cos(angle) - offset[1] * std::sin(angle), 0},
			                         40 * detail::deg, 10 * detail::deg);
			for (std::size_t i = 0; i < 3; ++i)
				e.positionEcefM[i] += d[i];
			entries.push_back(e);
		}
		const auto s = summarizeEntries(entries);
		EXPECT_NEAR(s["ellipse"]["semi_major_m"].get<double>(), 60., 1e-6);
		EXPECT_NEAR(s["ellipse"]["semi_minor_m"].get<double>(), 10., 1e-6);
		EXPECT_NEAR(s["ellipse"]["azimuth_deg"].get<double>(), 35., 1e-6);
	}

	TEST(EntrySummary, DegenerateAndUndefinedDirections)
	{
		auto vertical = entry(40, 10, {0, 0, -1000});
		const auto s = summarizeEntries(std::vector{vertical, vertical});
		EXPECT_EQ(s["heading_deg"]["mean"], nullptr);
		EXPECT_NEAR(s["flight_path_deg"]["mean"].get<double>(), -90., 1e-8);
		EXPECT_EQ(s["ellipse"]["semi_major_m"], 0.);
		EXPECT_EQ(s["ellipse"]["azimuth_deg"], nullptr);
		const auto stationary = summarizeEntries(std::vector{entry(40, 10, {0, 0, 0})});
		EXPECT_EQ(stationary["flight_path_deg"]["mean"], nullptr);
		const auto global = summarizeEntries(std::vector{entry(0, 0), entry(0, 180)});
		EXPECT_EQ(global["longitude_deg"]["mean"], nullptr);
		EXPECT_EQ(global["ellipse"]["semi_major_m"], nullptr);
		const auto wide = summarizeEntries(std::vector{entry(0, -10), entry(0, 10)});
		EXPECT_FALSE(wide["warnings"].empty());
		EXPECT_THROW(summarizeEntries(std::vector<EntryState>{}), std::invalid_argument);
		EXPECT_THROW(summarizeEntries(std::vector{entry(0, 0, {0, 0, 0}, -1)}), std::invalid_argument);
	}
}  // namespace
