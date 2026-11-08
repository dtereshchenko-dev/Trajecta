#include <gtest/gtest.h>

#include "orbital_core/hazard/PhaClassifier.h"
#include "orbital_core/Constants.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/Moid.h"

#include <cmath>
#include <filesystem>
#include <limits>

using namespace orbital_core::hazard;

TEST(PhaClassifierTest, IncludesBothThresholds)
{
	EXPECT_EQ(classifyPha(true, 22.0, 0.05), PhaStatus::Pha);
	EXPECT_EQ(classifyPha(true, std::nextafter(22.0, 23.0), 0.05), PhaStatus::NotPha);
	EXPECT_EQ(classifyPha(true, 22.0, std::nextafter(0.05, 1.0)), PhaStatus::NotPha);
	EXPECT_EQ(classifyPha(true, -1.0, 0.0), PhaStatus::Pha);
}

TEST(PhaClassifierTest, MissingOrInvalidInputsAreUnknown)
{
	const auto nan = std::numeric_limits<double>::quiet_NaN();
	const auto inf = std::numeric_limits<double>::infinity();
	EXPECT_EQ(classifyPha(true, nan, 0.01), PhaStatus::Unknown);
	EXPECT_EQ(classifyPha(true, 19.0, nan), PhaStatus::Unknown);
	EXPECT_EQ(classifyPha(true, inf, 0.01), PhaStatus::Unknown);
	EXPECT_EQ(classifyPha(true, 19.0, inf), PhaStatus::Unknown);
	EXPECT_EQ(classifyPha(true, 19.0, -0.01), PhaStatus::Unknown);
}

TEST(PhaClassifierTest, CometsAreNotApplicable)
{
	const auto nan = std::numeric_limits<double>::quiet_NaN();
	EXPECT_EQ(classifyPha(false, 19.0, 0.01), PhaStatus::NotApplicable);
	EXPECT_EQ(classifyPha(false, nan, nan), PhaStatus::NotApplicable);
}

TEST(PhaClassifierTest, ClassifiesComputedEarthMoidsForKnownOrbits)
{
	const auto originalDirectory = getExecutableDir();
	const auto data = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	if (!std::filesystem::exists(data / "linux_p1550p2650.440"))
		GTEST_SKIP() << "DE440 kernel not present";
	getExecutableDir() = data;
	DE440Source ephemeris;
	getExecutableDir() = originalDirectory;
	if (!ephemeris.isAvailable())
		GTEST_SKIP() << "DE440 kernel could not be opened";

	const auto makeOrbit = [](double a, double e, double i, double node, double arg) {
		Orbit orbit;
		orbit.setType(OrbitType::Elliptical);
		orbit.setSemiMajorAxis(a);
		orbit.setExcentricity(e);
		orbit.setInclination(i * Constants::DEG_TO_RAD);
		orbit.setLongitudeOfAscendingNode(node * Constants::DEG_TO_RAD);
		orbit.setArgumentOfPericentre(arg * Constants::DEG_TO_RAD);
		return orbit;
	};
	const auto apophis = makeOrbit(0.9224, 0.1914, 3.339, 204.024, 126.728);
	const auto eros = makeOrbit(1.458, 0.2227, 10.83, 304.30, 178.83);
	const auto apophisMoid = computeMoidsToAllPlanets(apophis, 2459266.348, &ephemeris)[Planet::Earth];
	const auto erosMoid = computeMoidsToAllPlanets(eros, 2459266.348, &ephemeris)[Planet::Earth];
	ASSERT_TRUE(std::isfinite(apophisMoid));
	ASSERT_TRUE(std::isfinite(erosMoid));
	EXPECT_EQ(classifyPha(true, 19.7, apophisMoid), PhaStatus::Pha);
	EXPECT_EQ(classifyPha(true, 11.16, erosMoid), PhaStatus::NotPha);
}
