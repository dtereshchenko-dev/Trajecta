#include <gtest/gtest.h>
#include "EphemeridesContainer.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Types.h"
#include <filesystem>

class EphemeridesContainerTest : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (getExecutableDir().empty()) {
			auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
			std::filesystem::path candidates[] = {
			    srcDir.parent_path() / "data", srcDir / "release", srcDir / "debug", srcDir / "data",
			    srcDir / "build/Release",
			};
			for (const auto& dir : candidates) {
				if (std::filesystem::exists(dir / "sun2019.txt") ||
				    std::filesystem::exists(dir / "linux_p1550p2650.440")) {
					getExecutableDir() = dir;
					break;
				}
			}
		}

		auto dir = getExecutableDir();
		hasSunData = std::filesystem::exists(dir / "sun2019.txt");
		hasEarthData = std::filesystem::exists(dir / "earth2019.txt");

		static DE440Source de440;
		if (de440.isAvailable() && !getDefaultEphemerisSource())
			setDefaultEphemerisSource(&de440);

		hasDE440 = getDefaultEphemerisSource() && getDefaultEphemerisSource()->isAvailable();
	}

	bool hasSunData = false;
	bool hasEarthData = false;
	bool hasDE440 = false;
};

TEST_F(EphemeridesContainerTest, SunEphemReturnsData)
{
	if (!hasSunData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	auto data = eph.get(2458634.5);

	EXPECT_EQ(eph.error(), EphemeridesContainer::Error::NoError);
	EXPECT_GT(data.R.magnitude(), 0.0);
}

TEST_F(EphemeridesContainerTest, SunPositionMagnitude)
{
	if (!hasSunData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	auto data = eph.get(2458634.5);

	double dist = data.R.magnitude();
	EXPECT_GT(dist, 0.9);
	EXPECT_LT(dist, 1.1);
}

TEST_F(EphemeridesContainerTest, SunVelocityNonZero)
{
	if (!hasSunData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	auto data = eph.get(2458634.5);

	if (hasDE440)
		EXPECT_GT(data.V.magnitude(), 0.0);
}

TEST_F(EphemeridesContainerTest, SunPositionContinuity)
{
	if (!hasSunData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	auto data1 = eph.get(2458634.5);
	auto data2 = eph.get(2458634.6);

	Vec3D diff = data2.R - data1.R;
	EXPECT_LT(diff.magnitude(), 0.01);
}

TEST_F(EphemeridesContainerTest, SunEphemDifferentDates)
{
	if (!hasSunData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	auto data1 = eph.get(2458634.5);
	auto data2 = eph.get(2458664.5);

	Vec3D diff = data2.R - data1.R;
	EXPECT_GT(diff.magnitude(), 0.01);
}

TEST_F(EphemeridesContainerTest, EarthEphemReturnsData)
{
	if (!hasEarthData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Earth);
	auto data = eph.get(2458634.5);

	EXPECT_EQ(eph.error(), EphemeridesContainer::Error::NoError);
	EXPECT_GT(data.R.magnitude(), 0.0);
}

TEST_F(EphemeridesContainerTest, EarthPositionMagnitude)
{
	if (!hasEarthData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Earth);
	auto data = eph.get(2458634.5);

	double dist = data.R.magnitude();
	EXPECT_GT(dist, 0.98);
	EXPECT_LT(dist, 1.02);
}

TEST_F(EphemeridesContainerTest, EarthOrbitalVelocity)
{
	if (!hasEarthData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Earth);
	auto data = eph.get(2458634.5);

	double v = data.V.magnitude();
	EXPECT_GT(v, 0.015);
	EXPECT_LT(v, 0.020);
}

TEST(EphemeridesContainerErrors, MissingSourceReturnsError)
{
	auto* saved = getDefaultEphemerisSource();
	setDefaultEphemerisSource(nullptr);

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	eph.get(2458634.5);
	EXPECT_EQ(eph.error(), EphemeridesContainer::Error::IOError);

	setDefaultEphemerisSource(saved);
}

TEST(EphemeridesContainerErrors, ErrorStringNotEmpty)
{
	auto* saved = getDefaultEphemerisSource();
	setDefaultEphemerisSource(nullptr);

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	eph.get(2458634.5);

	const std::string errStr = eph.errorString();
	EXPECT_FALSE(errStr.empty());
	EXPECT_NE(errStr, "No error");

	setDefaultEphemerisSource(saved);
}

TEST(EphemeridesContainerErrors, NoErrorStringForSuccess)
{
	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);
	EXPECT_EQ(eph.errorString(), "No error");
}

TEST_F(EphemeridesContainerTest, InterpolationSmoothness)
{
	if (!hasSunData && !hasDE440)
		GTEST_SKIP() << "No ephemeris data available";

	EphemeridesContainer eph(EphemeridesContainer::Planet::Sun);

	double jd0 = 2458634.5;
	auto d1 = eph.get(jd0);
	auto d2 = eph.get(jd0 + 0.5);
	auto d3 = eph.get(jd0 + 1.0);

	double avg_x = (d1.R.x + d3.R.x) / 2.0;
	EXPECT_NEAR(d2.R.x, avg_x, 0.01);
}
