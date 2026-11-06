#include <gtest/gtest.h>
#include "orbital_core/Constants.h"
#include "Vaisala.h"
#include "orbital_core/ObsNode.h"

#include <cmath>

TEST(Vaisala, ZeroArcReturnEmpty)
{
	ObsNode obs1;
	ObsNode obs2;
	obs1.JD = 2458634.0;
	obs2.JD = 2458634.0;

	obs1.vector = Vec3D(1.0, 0.0, 0.0);
	obs2.vector = Vec3D(1.0, 0.0, 0.0);
	obs1.obsVector = Vec3D(1.0, 0.0, 0.0);
	obs2.obsVector = Vec3D(1.0, 0.0, 0.0);

	auto result = Calculations::vaisala(obs1, obs2);
	EXPECT_FALSE(result.isValid());
}

TEST(Vaisala, SubDayArcRejected)
{
	ObsNode obs1;
	ObsNode obs2;
	obs1.JD = 2458634.0;
	obs2.JD = 2458634.005;

	obs1.vector = Vec3D(1.0, 0.0, 0.0);
	obs2.vector = Vec3D(std::cos(0.001), std::sin(0.001), 0.0);
	obs1.obsVector = Vec3D(-1.0, 0.0, 0.0);
	obs2.obsVector = Vec3D(-1.0, 0.0, 0.0);

	auto result = Calculations::vaisala(obs1, obs2);
	EXPECT_FALSE(result.isValid());
}

TEST(Vaisala, EmptySearchRangeReturnEmpty)
{
	ObsNode obs1;
	ObsNode obs2;
	obs1.JD = 2458634.0;
	obs2.JD = 2458664.0;

	obs1.vector = Vec3D(0.0, 0.0, 0.0);
	obs2.vector = Vec3D(0.0, 0.0, 0.0);
	obs1.obsVector = Vec3D(0.0, 0.0, 0.0);
	obs2.obsVector = Vec3D(0.0, 0.0, 0.0);

	auto result = Calculations::vaisala(obs1, obs2);
	EXPECT_FALSE(result.isValid());
}

TEST(Vaisala, SyntheticCircularTNOProducesSaneState)
{
	constexpr double a_target = 10.0;
	const double n_target = Constants::GAUSS_K / (a_target * std::sqrt(a_target));
	constexpr double epoch = 2458634.0;
	constexpr double dtDays = 30.0;

	auto makeObs = [](double jd, const Vec3D& targetHelio, const Vec3D& obsHelio) {
		ObsNode o;
		o.JD = jd;
		const Vec3D los = targetHelio - obsHelio;
		o.vector = los * (1.0 / los.magnitude());
		o.obsVector = obsHelio * (-1.0);
		return o;
	};

	auto obsAt = [&](double jd) {
		const double theta_obj = n_target * (jd - epoch);
		const double theta_earth = (Constants::TWOPI / 365.25) * (jd - epoch);

		const Vec3D target(a_target * std::cos(theta_obj), a_target * std::sin(theta_obj), 0.0);
		const Vec3D earth(std::cos(theta_earth), std::sin(theta_earth), 0.0);
		return makeObs(jd, target, earth);
	};

	const ObsNode o1 = obsAt(epoch);
	const ObsNode o2 = obsAt(epoch + dtDays);

	auto result = Calculations::vaisala(o1, o2, 5.0, 50.0);
	ASSERT_TRUE(result.isValid()) << "Väisälä failed on clean synthetic TNO";

	const double Rmag = result.first.magnitude();
	const double Vmag = result.second.magnitude();

	EXPECT_GT(Rmag, 7.0) << "Heliocentric distance unrealistically small";
	EXPECT_LT(Rmag, 13.0) << "Heliocentric distance unrealistically large";

	const double vExpected = 1.0 / std::sqrt(Rmag);
	EXPECT_GT(Vmag, 0.5 * vExpected);
	EXPECT_LT(Vmag, 2.0 * vExpected);

	EXPECT_TRUE(std::fabs(result.first.x) < 1e10);
	EXPECT_TRUE(std::fabs(result.first.y) < 1e10);
	EXPECT_TRUE(std::fabs(result.first.z) < 1e10);
	EXPECT_TRUE(std::fabs(result.second.x) < 1e10);
	EXPECT_TRUE(std::fabs(result.second.y) < 1e10);
	EXPECT_TRUE(std::fabs(result.second.z) < 1e10);
}

TEST(Vaisala, RestrictiveRangeStaysInsideOrFailsCleanly)
{
	constexpr double a_target = 30.0;
	const double n_target = Constants::GAUSS_K / (a_target * std::sqrt(a_target));
	constexpr double epoch = 2458634.0;
	constexpr double dtDays = 30.0;

	auto obsAt = [&](double jd) {
		const double theta_obj = n_target * (jd - epoch);
		const double theta_earth = (Constants::TWOPI / 365.25) * (jd - epoch);

		const Vec3D target(a_target * std::cos(theta_obj), a_target * std::sin(theta_obj), 0.0);
		const Vec3D earth(std::cos(theta_earth), std::sin(theta_earth), 0.0);

		ObsNode o;
		o.JD = jd;
		const Vec3D los = target - earth;
		o.vector = los * (1.0 / los.magnitude());
		o.obsVector = earth * (-1.0);
		return o;
	};

	const ObsNode o1 = obsAt(epoch);
	const ObsNode o2 = obsAt(epoch + dtDays);

	auto result = Calculations::vaisala(o1, o2, 5.0, 15.0);

	if (result.isValid()) {
		const double Rmag = result.first.magnitude();
		EXPECT_LT(Rmag, 18.0);
	}
}
