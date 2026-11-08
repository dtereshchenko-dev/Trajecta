#include <gtest/gtest.h>
#include "Herget.h"
#include "IodScoring.h"
#include "TwoBodyPropagation.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Types.h"
#include <array>
#include <cmath>
#include <vector>

TEST(IodScoring, KeplerResidualUsesObserverAndRetainsFailureSentinel)
{
	ObsNode observation;
	observation.JD = 2458634.0;
	observation.vector = Vec3D(1.0, 0.0, 0.0);
	observation.obsVector = Vec3D(0.0, 1.0, 0.0);
	const std::vector<ObsNode> observations{observation};
	const double epoch = observation.epochJdTdb();

	EXPECT_NEAR(Calculations::keplerResidualCost(Vec3D(2.0, 0.0, 0.0), Vec3D(), epoch, observations),
	            2.0 - 4.0 / std::sqrt(5.0), 1e-15);
	EXPECT_DOUBLE_EQ(Calculations::keplerResidualCost(Vec3D(), Vec3D(), epoch, observations), 1e30);
}

TEST(Herget, InsufficientObservationsReturnEmpty)
{
	std::vector<ObsNode> obs(2);
	obs[0].JD = 2458634.0;
	obs[1].JD = 2458635.0;

	auto result = Calculations::herget(obs);
	EXPECT_FALSE(result.isValid());
}

TEST(Herget, ZeroDurationReturnEmpty)
{
	std::vector<ObsNode> obs(4);
	for (auto& o : obs) {
		o.JD = 2458634.0;
		o.vector = Vec3D(1.0, 0.0, 0.0);
		o.obsVector = Vec3D(1.0, 0.0, 0.0);
	}

	auto result = Calculations::herget(obs);
	EXPECT_FALSE(result.isValid());
}

TEST(Herget, InvalidObservationsReturnEmpty)
{
	std::vector<ObsNode> obs(3);
	obs[0].JD = 2458634.0;
	obs[1].JD = 2458635.0;
	obs[2].JD = 2458636.0;
	for (auto& o : obs) {
		o.vector = Vec3D(0.0, 0.0, 0.0);
		o.obsVector = Vec3D(0.0, 0.0, 0.0);
	}

	auto result = Calculations::herget(obs);
	EXPECT_FALSE(result.isValid());
}

TEST(Herget, SyntheticCircularOrbitProducesSaneState)
{
	constexpr double a = 2.0;
	const double n = Constants::GAUSS_K / (a * std::sqrt(a));
	constexpr double epoch = 2458634.0;

	auto makeObs = [](double jd, const Vec3D& targetHelio, const Vec3D& obsHelio) {
		ObsNode o;
		o.JD = jd;
		const Vec3D los = (targetHelio - obsHelio);
		o.vector = los * (1.0 / los.magnitude());
		o.obsVector = obsHelio * (-1.0);
		return o;
	};

	std::vector<ObsNode> obs;
	obs.reserve(5);
	for (int i = 0; i < 5; ++i) {
		const double jd = epoch + i * 5.0;
		const double theta_obj = n * (jd - epoch);
		const double theta_earth = (Constants::TWOPI / 365.25) * (jd - epoch);

		const Vec3D target(a * std::cos(theta_obj), a * std::sin(theta_obj), 0.0);
		const Vec3D earth(std::cos(theta_earth), std::sin(theta_earth), 0.0);
		obs.push_back(makeObs(jd, target, earth));
	}

	auto result = Calculations::herget(obs);
	ASSERT_TRUE(result.isValid()) << "Herget failed to converge on clean synthetic data";

	const double Rmag = result.first.magnitude();
	const double Vmag = result.second.magnitude();

	EXPECT_GT(Rmag, 1.0) << "Heliocentric distance unrealistically small";
	EXPECT_LT(Rmag, 4.0) << "Heliocentric distance unrealistically large";

	const double vExpected = 1.0 / std::sqrt(Rmag);
	EXPECT_GT(Vmag, 0.3 * vExpected);
	EXPECT_LT(Vmag, 3.0 * vExpected);

	EXPECT_TRUE(std::isfinite(result.first.x) && std::isfinite(result.first.y) && std::isfinite(result.first.z));
	EXPECT_TRUE(std::isfinite(result.second.x) && std::isfinite(result.second.y) && std::isfinite(result.second.z));
}

TEST(Herget, HyperbolicPropagationMatchesKeplerReferenceInBothTimeDirections)
{
	struct Reference {
		double days;
		Vec3D position;
		Vec3D canonicalVelocity;
	};

	const std::array<Reference, 2> references{{{365.0, Vec3D(-9.675519494550278, 11.150728069623561, 0.0),
	                                            Vec3D(-1.5106031508274609, 1.689245209541408, 0.0)},
	                                           {1000.0, Vec3D(-26.0875279772301, 29.501305749684647, 0.0),
	                                            Vec3D(-1.4982396095725208, 1.675130922578393, 0.0)}}};
	const StateVectors initial(Vec3D(0.1, 0.0, 0.0), Vec3D(0.0, 5.0, 0.0));
	for (const auto& reference : references) {
		for (const double sign : {-1.0, 1.0}) {
			const double days = sign * reference.days;
			SCOPED_TRACE(days);
			const auto propagated = Calculations::propagateKeplerState(initial, days);
			ASSERT_TRUE(propagated.isValid());
			const Vec3D expectedR(reference.position.x, sign * reference.position.y, 0.0);
			const Vec3D expectedV(sign * reference.canonicalVelocity.x, reference.canonicalVelocity.y, 0.0);
			EXPECT_LT((propagated.first - expectedR).magnitude(), 5.0e-11);
			EXPECT_LT((propagated.second - expectedV).magnitude(), 5.0e-12);

			const auto restored = Calculations::propagateKeplerState(propagated, -days);
			ASSERT_TRUE(restored.isValid());
			EXPECT_LT((restored.first - initial.first).magnitude(), 1.0e-10);
			EXPECT_LT((restored.second - initial.second).magnitude(), 1.0e-9);
		}
	}
}

TEST(Herget, TinyHyperbolicTimeStepsPreserveTheInitialState)
{
	const StateVectors initial(Vec3D(0.1, 0.0, 0.0), Vec3D(0.0, 5.0, 0.0));
	for (const double days : {-1.0e-40, 1.0e-40}) {
		SCOPED_TRACE(days);
		const auto result = Calculations::propagateKeplerState(initial, days);
		ASSERT_TRUE(result.isValid());
		EXPECT_LT((result.first - initial.first).magnitude(), 1e-15);
		EXPECT_LT((result.second - initial.second).magnitude(), 1e-15);
	}
}

TEST(Herget, CometaryPropagationRoundTripsAcrossExactParabola)
{
	for (const double e : {1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10}) {
		const StateVectors initial(Vec3D(2.0, 0.0, 0.0), Vec3D(0.0, std::sqrt((1.0 + e) / 2.0), 0.0));
		for (const double days : {-1000.0, -365.0, 365.0, 1000.0}) {
			SCOPED_TRACE(::testing::Message() << "e=" << e << " days=" << days);
			const auto propagated = Calculations::propagateKeplerState(initial, days);
			ASSERT_TRUE(propagated.isValid());
			const auto restored = Calculations::propagateKeplerState(propagated, -days);
			ASSERT_TRUE(restored.isValid());
			EXPECT_LT((restored.first - initial.first).magnitude(), 1.0e-10);
			EXPECT_LT((restored.second - initial.second).magnitude(), 1.0e-11);
			EXPECT_NEAR(Calculations::specificEnergy(propagated.first, propagated.second, 1.0),
			            Calculations::specificEnergy(initial.first, initial.second, 1.0), 1.0e-11);
			EXPECT_LT((propagated.first.cross(propagated.second) - initial.first.cross(initial.second)).magnitude(),
			          1.0e-10);
		}
	}
}

TEST(Herget, ExactParabolicPropagationMatchesBarkerEquation)
{
	const StateVectors initial(Vec3D(2.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0));
	for (const double days : {-1000.0, -365.0, 365.0, 1000.0}) {
		SCOPED_TRACE(days);
		const double tau = Constants::GAUSS_K * days;
		const double d = 2.0 * std::sinh(std::asinh(3.0 * tau / 8.0) / 3.0);
		const Vec3D expectedR(2.0 * (1.0 - d * d), 4.0 * d, 0.0);
		const Vec3D expectedV(-d / (1.0 + d * d), 1.0 / (1.0 + d * d), 0.0);
		const auto propagated = Calculations::propagateKeplerState(initial, days);
		ASSERT_TRUE(propagated.isValid());
		EXPECT_LT((propagated.first - expectedR).magnitude(), 1.0e-10);
		EXPECT_LT((propagated.second - expectedV).magnitude(), 1.0e-11);
	}
}

TEST(Herget, ConvergedEllipticPropagationMatchesEccentricAnomaly)
{
	constexpr double a = 2.0, e = 0.2;
	const double q = a * (1.0 - e);
	const StateVectors initial(Vec3D(q, 0.0, 0.0), Vec3D(0.0, std::sqrt((1.0 + e) / q), 0.0));
	for (const double anomaly : {-1.1, 1.1}) {
		SCOPED_TRACE(anomaly);
		const double days = a * std::sqrt(a) * (anomaly - e * std::sin(anomaly)) / Constants::GAUSS_K;
		const double beta = std::sqrt(1.0 - e * e);
		const Vec3D expectedR(a * (std::cos(anomaly) - e), a * beta * std::sin(anomaly), 0.0);
		const Vec3D expectedV =
		    Vec3D(-std::sin(anomaly), beta * std::cos(anomaly), 0.0) / (std::sqrt(a) * (1.0 - e * std::cos(anomaly)));
		const auto propagated = Calculations::propagateKeplerState(initial, days);
		ASSERT_TRUE(propagated.isValid());
		EXPECT_LT((propagated.first - expectedR).magnitude(), 1.0e-11);
		EXPECT_LT((propagated.second - expectedV).magnitude(), 1.0e-11);
	}
}
