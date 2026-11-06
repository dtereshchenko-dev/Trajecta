#include "orbital_core/ObservationDiagnostics.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

namespace {
	ObsNode observation(const char* id, double seconds, double ra, double dec)
	{
		ObsNode o;
		o.id = id;
		o.obsCode = "248";
		o.JD = 2440000.0;
		o.hasSplitUtc = true;
		o.utcMjdDay = 40000;
		o.utcSecondsOfDay = seconds;
		o.RA = ra;
		o.DEC = dec;
		return o;
	}
}  // namespace

TEST(ObservationDiagnostics, SameStationAndTimeUsesSignedAngularSeparation)
{
	std::vector<ObsNode> observations{observation("a", 100, 1e-8, 0.4),
	                                  observation("b", 100, Constants::TWOPI - 1e-8, 0.4 - 1e-8),
	                                  observation("c", 100.02, 1e-8, 0.4), observation("d", 100, 1e-8, 0.4)};
	observations.back().obsCode = "083";
	const auto result = orbital_core::findCoincidentObservations(observations);
	EXPECT_TRUE(result.empty());
	EXPECT_EQ(observations[0].id, "a");
	EXPECT_DOUBLE_EQ(observations[0].RA, 1e-8);
}

TEST(ObservationDiagnostics, LargeTimeGroupDoesNotCreateQuadraticNumberOfPairs)
{
	std::vector<ObsNode> observations;
	for (int i = 0; i < 1000; ++i) {
		auto o = observation("", 0, 0.5, 0.2);
		o.id = std::to_string(i);
		observations.push_back(std::move(o));
	}
	const auto result = orbital_core::findCoincidentObservations(observations);
	ASSERT_EQ(result.size(), 1u);
	EXPECT_EQ(result.front().observations.size(), observations.size());
	ASSERT_EQ(result.front().issues.size(), 1u);
	EXPECT_EQ(result.front().issues.front(), orbital_core::ObservationIssue::CoincidentPosition);
	EXPECT_DOUBLE_EQ(result.front().maxSeparationArcsec, 0.0);
}

TEST(ObservationDiagnostics, DoesNotCompareSpacecraftOrDifferentTimeScales)
{
	std::vector<ObsNode> observations{observation("a", 100, 1, 0.4), observation("b", 100, 1, 0.4),
	                                  observation("c", 100, 1, 0.4)};
	observations[1].isSatellite = true;
	observations[2].timeScale = orbital_core::time::Scale::UT1;
	EXPECT_TRUE(orbital_core::findCoincidentObservations(observations).empty());
}

TEST(ObservationDiagnostics, DifferentPositionsAreReportedWithoutExclusion)
{
	std::vector<ObsNode> observations{observation("a", 100, 1, 0.4), observation("b", 102.5, 1, 0.401)};
	const auto result = orbital_core::findCoincidentObservations(observations);
	ASSERT_EQ(result.size(), 1u);
	ASSERT_EQ(result[0].observations.size(), 2u);
	EXPECT_EQ(result[0].issues.front(), orbital_core::ObservationIssue::NearbyObservationConflict);
	EXPECT_GT(result[0].maxSeparationArcsec, 200.0);
	EXPECT_EQ(observations.size(), 2u);
}

TEST(ObservationDiagnostics, GeometryFlagsGroundObservationsButNotSpacecraft)
{
	auto o = observation("ground", 43200, Constants::PI, 0);
	o.utcMjdDay = 51544;
	o.JD = 2451545.0;
	o.timeScale = orbital_core::time::Scale::TT;
	o.cosPhi = 1.0;
	o.sinPhi = 0.0;
	o.LST = 0.0;
	o.obsVector = Helpers::ecvToEcl(Vec3D(1, 0, 0));
	auto spacecraft = o;
	spacecraft.id = "spacecraft";
	spacecraft.isSatellite = true;
	const auto result = orbital_core::diagnoseObservationVisibility(std::vector<ObsNode>{o, spacecraft});
	ASSERT_EQ(result.size(), 1u);
	EXPECT_EQ(result[0].observation, "ground");
	ASSERT_TRUE(result[0].objectAltitudeDeg);
	ASSERT_TRUE(result[0].sunAltitudeDeg);
	EXPECT_LT(*result[0].objectAltitudeDeg, -89.0);
	EXPECT_GT(*result[0].sunAltitudeDeg, 89.0);
	EXPECT_NE(std::ranges::find(result[0].issues, orbital_core::ObservationIssue::BelowHorizon),
	          result[0].issues.end());
	EXPECT_NE(std::ranges::find(result[0].issues, orbital_core::ObservationIssue::Daylight), result[0].issues.end());
}

TEST(ObservationDiagnostics, HistoricalMidnightAloneIsNotAProblem)
{
	std::vector<ObsNode> observations;
	for (int i = 0; i < 3; ++i) {
		auto o = observation("", 0, 1, 0.4);
		o.id = std::to_string(i);
		o.timeScale = orbital_core::time::Scale::UT1;
		o.utcMjdDay += i;
		observations.push_back(std::move(o));
	}
	const auto result = orbital_core::diagnoseObservationVisibility(observations);
	EXPECT_TRUE(result.empty());
	EXPECT_EQ(observations.size(), 3u);
}
