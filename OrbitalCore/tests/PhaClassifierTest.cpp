#include <gtest/gtest.h>

#include "orbital_core/Constants.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/Moid.h"
#include "orbital_core/hazard/PhaClassifier.h"

#include <cmath>
#include <filesystem>
#include <memory>

using namespace orbital_core::hazard;

class PhaClassifierTest : public ::testing::Test {
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

		m_de440 = std::make_unique<DE440Source>();
		if (m_de440->isAvailable())
			setDefaultEphemerisSource(m_de440.get());

		hasDE440 = m_de440->isAvailable();
	}

	void TearDown() override { setDefaultEphemerisSource(nullptr); }

	static Orbit_p makeKepler(double aAu, double e, double iDeg, double OmDeg, double wDeg)
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

	std::unique_ptr<DE440Source> m_de440;
	bool hasDE440{false};
};

TEST_F(PhaClassifierTest, ApophisIsPha)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	auto apophis = makeKepler(0.9224, 0.1914, 3.339, 204.024, 126.728);
	const auto moids = computeMoidsToAllPlanets(*apophis, 2459266.348);
	const auto c = classifyPha(*apophis, 19.7, moids);

	EXPECT_TRUE(c.perihelionOk);
	EXPECT_TRUE(c.moidOk) << "Apophis MOID=" << c.moidEarthAu;
	EXPECT_TRUE(c.absoluteMagOk);
	EXPECT_TRUE(c.isPha);
	EXPECT_NEAR(c.semiMajorAxisAu, 0.9224, 1.0e-9);
	EXPECT_LT(c.moidEarthAu, kPhaMoidThresholdAu);
}

TEST_F(PhaClassifierTest, AmorClass_2007VK184_IsPha)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	auto vk184 = makeKepler(1.726, 0.570, 1.221, 253.545, 73.681);
	const auto moids = computeMoidsToAllPlanets(*vk184, 2459266.348);
	const auto c = classifyPha(*vk184, 17.6, moids);

	EXPECT_TRUE(c.perihelionOk) << "q=" << c.perihelionAu;
	EXPECT_TRUE(c.absoluteMagOk);
	EXPECT_NEAR(c.perihelionAu, 0.742, 1.0e-3);
	EXPECT_GT(c.semiMajorAxisAu, kPhaPerihelionThresholdAu);
}

TEST_F(PhaClassifierTest, CeresIsNotPha)
{
	auto ceres = makeKepler(2.7658533811, 0.07826903, 10.587999, 80.269653, 73.735909);
	const auto c = classifyPha(*ceres, 3.34, PlanetMoidTable{});

	EXPECT_FALSE(c.perihelionOk);
	EXPECT_FALSE(c.moidOk);
	EXPECT_FALSE(c.isPha);
	EXPECT_TRUE(std::isnan(c.moidEarthAu));
}

TEST_F(PhaClassifierTest, ErosIsNeaButNotPha)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	auto eros = makeKepler(1.458, 0.2227, 10.83, 304.30, 178.83);
	const auto moids = computeMoidsToAllPlanets(*eros, 2459266.348);
	const auto c = classifyPha(*eros, 11.16, moids);

	EXPECT_TRUE(c.perihelionOk);
	EXPECT_TRUE(c.absoluteMagOk);
	EXPECT_FALSE(c.moidOk) << "Eros MOID=" << c.moidEarthAu;
	EXPECT_FALSE(c.isPha);
	EXPECT_NEAR(c.perihelionAu, 1.133, 1.0e-2);
}

TEST_F(PhaClassifierTest, NanMagnitude_GatesIndependent)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	auto apophis = makeKepler(0.9224, 0.1914, 3.339, 204.024, 126.728);
	const auto moids = computeMoidsToAllPlanets(*apophis, 2459266.348);
	const auto c = classifyPha(*apophis, std::numeric_limits<double>::quiet_NaN(), moids);

	EXPECT_TRUE(c.perihelionOk);
	EXPECT_TRUE(c.moidOk);
	EXPECT_FALSE(c.absoluteMagOk);
	EXPECT_FALSE(c.isPha);
}

TEST(PhaClassifierCachedTest, UsesSuppliedEarthMoid)
{
	Orbit orbit;
	orbit.setType(OrbitType::Elliptical);
	orbit.setSemiMajorAxis(1.0);
	orbit.setExcentricity(0.2);
	PlanetMoidTable moids;
	moids[Planet::Earth] = 0.01;
	moids[Planet::Mars] = 1.0;

	const auto c = classifyPha(orbit, 19.7, moids);
	EXPECT_DOUBLE_EQ(c.moidEarthAu, moids[Planet::Earth]);
	EXPECT_TRUE(c.isPha);

	moids[Planet::Earth] = 0.1;
	moids[Planet::Mars] = 0.001;
	EXPECT_FALSE(classifyPha(orbit, 19.7, moids).isPha);
}

TEST(PhaClassifierCachedTest, MissingEarthMoidDoesNotPass)
{
	Orbit orbit;
	orbit.setType(OrbitType::Elliptical);
	orbit.setSemiMajorAxis(1.0);
	orbit.setExcentricity(0.2);
	PlanetMoidTable moids;
	moids[Planet::Mars] = 0.001;

	const auto c = classifyPha(orbit, 19.7, moids);
	EXPECT_TRUE(c.perihelionOk);
	EXPECT_TRUE(c.absoluteMagOk);
	EXPECT_TRUE(std::isnan(c.moidEarthAu));
	EXPECT_FALSE(c.moidOk);
	EXPECT_FALSE(c.isPha);

	moids[Planet::Earth] = std::numeric_limits<double>::infinity();
	EXPECT_FALSE(classifyPha(orbit, 19.7, moids).moidOk);
}

TEST(PhaClassifierCachedTest, StrictThresholdsArePreserved)
{
	Orbit orbit;
	orbit.setType(OrbitType::Elliptical);
	orbit.setSemiMajorAxis(1.0);
	orbit.setExcentricity(0.0);
	PlanetMoidTable moids;
	moids[Planet::Earth] = std::nextafter(kPhaMoidThresholdAu, 0.0);
	const double bright = std::nextafter(kPhaAbsoluteMagThreshold, 0.0);
	EXPECT_TRUE(classifyPha(orbit, bright, moids).isPha);
	EXPECT_FALSE(classifyPha(orbit, kPhaAbsoluteMagThreshold, moids).isPha);

	moids[Planet::Earth] = kPhaMoidThresholdAu;
	EXPECT_FALSE(classifyPha(orbit, bright, moids).isPha);

	moids[Planet::Earth] = 0.01;
	orbit.setSemiMajorAxis(std::nextafter(kPhaPerihelionThresholdAu, 0.0));
	EXPECT_TRUE(classifyPha(orbit, bright, moids).isPha);
	orbit.setSemiMajorAxis(kPhaPerihelionThresholdAu);
	const auto c = classifyPha(orbit, bright, moids);
	EXPECT_FALSE(c.perihelionOk);
	EXPECT_FALSE(c.isPha);
	EXPECT_TRUE(std::isnan(c.moidEarthAu));
}
