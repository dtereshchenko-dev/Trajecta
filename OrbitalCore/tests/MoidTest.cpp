#include <gtest/gtest.h>

#include "orbital_core/Constants.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/Moid.h"
#include "orbital_core/hazard/Planet.h"

#include <cmath>
#include <filesystem>
#include <memory>

using namespace orbital_core::hazard;

namespace {
	Orbit_p makeKepler(double aAu, double e, double iDeg, double OmDeg, double wDeg)
	{
		auto o = std::make_shared<Orbit>();
		o->setType(OrbitType::Elliptical);
		o->setSemiMajorAxis(aAu);
		o->setExcentricity(e);
		o->setInclination(iDeg * Constants::DEG_TO_RAD);
		o->setLongitudeOfAscendingNode(OmDeg * Constants::DEG_TO_RAD);
		o->setArgumentOfPericentre(wDeg * Constants::DEG_TO_RAD);
		return o;
	}
}  // namespace

TEST(MoidTest, CoplanarCircles)
{
	auto inner = makeKepler(1.0, 0.0, 0.0, 0.0, 0.0);
	auto outer = makeKepler(2.5, 0.0, 0.0, 0.0, 0.0);

	const auto r = computeMoid(*inner, *outer);
	EXPECT_TRUE(r.converged);
	EXPECT_TRUE(r.isValid());
	EXPECT_NEAR(r.moidAu, 1.5, 1.0e-9);
}

TEST(MoidTest, IdenticalOrbits)
{
	auto a = makeKepler(2.0, 0.3, 15.0, 30.0, 60.0);

	const auto r = computeMoid(*a, *a);
	EXPECT_LT(r.moidAu, 1.0e-8);
}

TEST(MoidTest, ApophisVsEarth_PhaScale)
{
	auto apophis = makeKepler(0.9224, 0.1914, 3.339, 204.024, 126.728);
	auto earth = makeKepler(1.0000, 0.0167, 0.000, 0.000, 102.937);

	const auto r = computeMoid(*apophis, *earth);
	EXPECT_TRUE(r.converged);
	EXPECT_LT(r.moidAu, 0.001);
	EXPECT_GT(r.moidAu, 0.0);
}

TEST(MoidTest, HyperbolicReturnsNan)
{
	auto borisov = makeKepler(-1.0, 3.36, 44.0, 308.2, 209.1);
	borisov->setType(OrbitType::Hyperbolic);

	auto earth = makeKepler(1.0, 0.0167, 0.0, 0.0, 102.937);

	const auto r = computeMoid(*borisov, *earth);
	EXPECT_FALSE(r.converged);
	EXPECT_FALSE(r.isValid());
	EXPECT_TRUE(std::isnan(r.moidAu));
}

class MoidMultiPlanetTest : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (getExecutableDir().empty()) {
			const auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
			const std::filesystem::path candidates[] = {
			    srcDir.parent_path() / "data",
			    srcDir / "data",
			    srcDir / "release",
			    srcDir / "debug",
			};
			for (const auto& d : candidates) {
				if (std::filesystem::exists(d / "linux_p1550p2650.440")) {
					getExecutableDir() = d;
					break;
				}
			}
		}

		static DE440Source de440;
		if (de440.isAvailable() && !getDefaultEphemerisSource())
			setDefaultEphemerisSource(&de440);

		hasDE440 = getDefaultEphemerisSource() && getDefaultEphemerisSource()->isAvailable();
	}

	bool hasDE440{false};
};

TEST_F(MoidMultiPlanetTest, CeresAllPlanets)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 not available";

	auto ceres = makeKepler(2.7658533811, 0.07826903, 10.587999, 80.269653, 73.735909);

	const auto t = computeMoidsToAllPlanets(*ceres, 2459266.348);
	EXPECT_TRUE(t.any());

	struct Expected {
		Planet planet;
		double moidAu;
		double tolAu;
	};

	const Expected ref[] = {
	    {Planet::Mercury, 2.175611, 5.0e-3}, {Planet::Venus, 1.837678, 5.0e-3},    {Planet::Earth, 1.585189, 5.0e-3},
	    {Planet::Mars, 0.931129, 5.0e-3},    {Planet::Jupiter, 2.092539, 5.0e-2},  {Planet::Saturn, 6.325274, 1.0e-1},
	    {Planet::Uranus, 15.764471, 1.5e-1}, {Planet::Neptune, 26.951993, 2.0e-1},
	};

	for (const auto& e : ref) {
		EXPECT_NEAR(t[e.planet], e.moidAu, e.tolAu) << "planet=" << longLabel(e.planet) << " got=" << t[e.planet];
	}
}

TEST(MoidMultiPlanetTest_NoEph, NullEphemerisReturnsAllNan)
{
	auto* saved = getDefaultEphemerisSource();
	setDefaultEphemerisSource(nullptr);

	auto target = makeKepler(2.5, 0.1, 5.0, 30.0, 60.0);
	const auto t = computeMoidsToAllPlanets(*target, 2459266.348, nullptr);

	EXPECT_FALSE(t.any());
	for (auto p : kAllPlanets)
		EXPECT_TRUE(std::isnan(t[p])) << "planet=" << longLabel(p);

	setDefaultEphemerisSource(saved);
}
