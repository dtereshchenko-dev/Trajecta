#include <gtest/gtest.h>

#include "orbital_core/forces/IForceModel.h"
#include "orbital_core/forces/SunGravity.h"
#include "orbital_core/forces/EihRelativity.h"
#include "orbital_core/forces/SunJ2.h"
#include "orbital_core/forces/ForceField.h"
#include "orbital_core/perturbations/AsteroidPerturbations.h"
#include "orbital_core/Constants.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>

namespace ocf = orbital_core::forces;

static_assert(ocf::ForceModel<ocf::SunGravity>);
static_assert(ocf::ForceModel<ocf::EihRelativity>);
static_assert(ocf::ForceModel<ocf::SunJ2>);

TEST(AsteroidPerturbations, ExcludesOnlyTheIdentifiedAsteroidFromFallback)
{
	ocf::ForceContext input{};
	input.R = Vec3D(2.0, 0.5, 0.1);
	input.setEpochMjd(60000.0);
	input.bodyMask = 0x0700;
	input.computePartials = true;
	const auto evaluate = [](ocf::ForceContext context) {
		Vec3D acceleration;
		Eigen::Matrix3d dr = Eigen::Matrix3d::Zero(), dv = Eigen::Matrix3d::Zero();
		ocf::AsteroidPerturbations{}.evaluate(context, acceleration, dr, dv);
		return std::pair{acceleration, dr};
	};
	const auto full = evaluate(input);
	input.targetNaifId = 2000001;
	const auto withoutCeres = evaluate(input);
	input.targetNaifId = 0;
	input.bodyMask &= static_cast<std::uint16_t>(~(1u << 8));
	const auto manual = evaluate(input);
	EXPECT_DOUBLE_EQ((withoutCeres.first - manual.first).magnitude(), 0.0);
	EXPECT_TRUE(withoutCeres.second.isApprox(manual.second, 0.0));
	EXPECT_GT((full.first - withoutCeres.first).magnitude(), 0.0);
	EXPECT_GT(withoutCeres.first.magnitude(), 0.0);
}

TEST(SunJ2, MatchesReferenceEquatorialAndPolarSpecialCases)
{
	constexpr double radius = ocf::SunJ2::kRadiusAu;
	constexpr double j2 = ocf::SunJ2::kJ2;
	const double factor = 1.5 * Constants::GM_SUN * j2 * radius * radius;

	Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dadv = Eigen::Matrix3d::Zero();

	ocf::ForceContext equatorial;
	const Vec3D pole = ocf::SunJ2::poleEcliptic();
	const Vec3D transverse(-pole.y, pole.x, 0.0);
	equatorial.R = transverse / transverse.magnitude();
	equatorial.computePartials = false;
	Vec3D equatorialAcceleration;

	ocf::SunJ2{}.evaluate(equatorial, equatorialAcceleration, dadr, dadv);
	EXPECT_LT((equatorialAcceleration + equatorial.R * factor).magnitude(), factor * 1e-14);
	const double c = std::cos(Constants::EARTH_INCL), s = std::sin(Constants::EARTH_INCL);
	const Vec3D eq(pole.x, c * pole.y - s * pole.z, s * pole.y + c * pole.z);
	EXPECT_NEAR(std::atan2(eq.y, eq.x) / Constants::DEG_TO_RAD + 360.0, 286.13, 1e-12);
	EXPECT_NEAR(std::asin(eq.z) / Constants::DEG_TO_RAD, 63.87, 1e-12);

	ocf::ForceContext polar;
	polar.R = pole;
	polar.computePartials = false;
	Vec3D polarAcceleration;

	ocf::SunJ2{}.evaluate(polar, polarAcceleration, dadr, dadv);
	EXPECT_LT((polarAcceleration - pole * (2.0 * factor)).magnitude(), factor * 1e-14);
	EXPECT_TRUE(dadr.isZero());
	EXPECT_TRUE(dadv.isZero());
}

namespace {

	template <class Force>
	void checkForcePartials(ocf::ForceContext ctx, double positionStep)
	{
		const Force force;
		Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero(), dadv = Eigen::Matrix3d::Zero();
		Vec3D acceleration;
		ctx.computePartials = true;
		force.evaluate(ctx, acceleration, dadr, dadv);
		Eigen::Matrix<double, 3, 6> numerical;
		for (int col = 0; col < 6; ++col) {
			const double h = col < 3 ? positionStep : 1e-6;
			const auto sample = [&](double offset) {
				auto perturbed = ctx;
				perturbed.computePartials = false;
				double* coordinate = col < 3 ? (col == 0   ? &perturbed.R.x
				                                : col == 1 ? &perturbed.R.y
				                                           : &perturbed.R.z)
				                             : (col == 3   ? &perturbed.V.x
				                                : col == 4 ? &perturbed.V.y
				                                           : &perturbed.V.z);
				*coordinate += offset;
				Eigen::Matrix3d untouchedR = Eigen::Matrix3d::Identity();
				Eigen::Matrix3d untouchedV = Eigen::Matrix3d::Identity();
				Vec3D a;
				force.evaluate(perturbed, a, untouchedR, untouchedV);
				EXPECT_TRUE(untouchedR.isIdentity(0.0));
				EXPECT_TRUE(untouchedV.isIdentity(0.0));
				return Eigen::Vector3d(a.x, a.y, a.z);
			};
			numerical.col(col) = (8.0 * (sample(h) - sample(-h)) - (sample(2.0 * h) - sample(-2.0 * h))) / (12.0 * h);
		}
		EXPECT_LT((dadr - numerical.leftCols<3>()).norm(), std::max(dadr.norm(), 1e-30) * 2e-8);
		if (dadv.isZero(0.0))
			EXPECT_TRUE(numerical.rightCols<3>().isZero(1e-22));
		else
			EXPECT_LT((dadv - numerical.rightCols<3>()).norm(), dadv.norm() * 2e-8);

		ctx.computePartials = false;
		Vec3D stateOnly;
		Eigen::Matrix3d ignoredR = Eigen::Matrix3d::Zero(), ignoredV = Eigen::Matrix3d::Zero();
		force.evaluate(ctx, stateOnly, ignoredR, ignoredV);
		EXPECT_DOUBLE_EQ(stateOnly.x, acceleration.x);
		EXPECT_DOUBLE_EQ(stateOnly.y, acceleration.y);
		EXPECT_DOUBLE_EQ(stateOnly.z, acceleration.z);
		ctx.computePartials = true;
		ignoredR.setIdentity();
		ignoredV.setIdentity();
		force.evaluate(ctx, stateOnly, ignoredR, ignoredV);
		EXPECT_TRUE(ignoredR.isApprox(Eigen::Matrix3d::Identity() + dadr, 1e-15));
		EXPECT_TRUE(ignoredV.isApprox(Eigen::Matrix3d::Identity() + dadv, 1e-15));
	}

	ocf::ForceContext relativisticContext()
	{
		ocf::ForceContext ctx;
		ctx.havePlanets = ctx.havePlanetVelocities = true;
		ctx.bodyMask = 0x18ff;
		ctx.ephemerisTerms.sunVelocity = Vec3D(4e-6, -7e-6, 2e-6);
		ctx.ephemerisTerms.plutoPosition = Vec3D(25.0, -20.0, -3.0);
		ctx.ephemerisTerms.plutoVelocity = Vec3D(0.002, 0.001, -0.0002);
		ctx.plutoGm = ocf::detail::GM_PLUTO;
		ctx.sunGm = Constants::GM_SUN * 1.001;
		for (int i = 0; i < 8; ++i) {
			const double angle = 0.8 * i + 0.4;
			const double radius = 0.4 + 0.7 * i;
			ctx.planetPositions[i] =
			    Vec3D(radius * std::cos(angle), radius * std::sin(angle), 0.03 * radius * std::sin(2 * angle));
			ctx.planetVelocities[i] = Vec3D(-0.02 * std::sin(angle) / std::sqrt(radius),
			                                0.02 * std::cos(angle) / std::sqrt(radius), 0.0004 * std::cos(2 * angle));
			ctx.planetGms[i] = ocf::detail::kPerturbers[i].gm * 1.002;
		}
		ctx.moonPosition = ctx.planetPositions[2] + Vec3D(0.002, 0.001, -0.0004);
		ctx.moonVelocity = ctx.planetVelocities[2] + Vec3D(-0.0001, 0.0002, 0.00005);
		return ctx;
	}

}  // namespace

TEST(EihRelativity, StatePartialsMatchIndependentDifferences)
{
	for (const bool planets : {false, true}) {
		for (const Vec3D position :
		     {Vec3D(0.8, -0.3, 0.2), Vec3D(5.0, 2.0, -0.5), Vec3D(0.01, -0.002, 0.004), Vec3D(0.0, 0.0, 1.0)}) {
			auto ctx = relativisticContext();
			if (!planets)
				ctx.bodyMask = 0;
			ctx.R = position;
			ctx.V = Vec3D(0.004, 0.018, -0.003);
			SCOPED_TRACE(::testing::Message() << "planets=" << planets << " R=" << position.x);
			checkForcePartials<ocf::EihRelativity>(ctx, position.magnitude() * 1e-5);
		}
	}
	auto close = relativisticContext();
	close.R = close.moonPosition + Vec3D(0.0003, -0.0002, 0.0001);
	close.V = Vec3D(-0.004, 0.018, 0.003);
	checkForcePartials<ocf::EihRelativity>(close, 1e-8);
	close.bodyMask = (1u << 2) | orbital_core::planetaryBodyMask(EphemerisPlanet::Moon);
	checkForcePartials<ocf::EihRelativity>(close, 1e-8);
}

TEST(EihRelativity, SunOnlyMatchesSchwarzschildAccelerationAndVelocityDerivative)
{
	ocf::ForceContext ctx;
	ctx.R = Vec3D(0.7, -0.2, 0.1);
	ctx.V = Vec3D(0.01, 0.02, -0.003);
	Vec3D acceleration;
	Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero(), dadv = Eigen::Matrix3d::Zero();
	ocf::EihRelativity{}.evaluate(ctx, acceleration, dadr, dadv);
	const double radius = ctx.R.magnitude();
	const double factor =
	    Constants::GM_SUN / (radius * radius * radius * Constants::AU_PER_DAY * Constants::AU_PER_DAY);
	const Vec3D expected =
	    (ctx.R * (4 * Constants::GM_SUN / radius - ctx.V.dot(ctx.V)) + ctx.V * (4 * ctx.R.dot(ctx.V))) * factor;
	EXPECT_LT((acceleration - expected).magnitude(), expected.magnitude() * 1e-14);
	const Eigen::Vector3d r(ctx.R.x, ctx.R.y, ctx.R.z), v(ctx.V.x, ctx.V.y, ctx.V.z);
	const Eigen::Matrix3d expectedV =
	    factor * (-2 * r * v.transpose() + 4 * v * r.transpose() + 4 * r.dot(v) * Eigen::Matrix3d::Identity());
	EXPECT_TRUE(dadv.isApprox(expectedV, 1e-14));
}

TEST(SunJ2, StatePartialsMatchIndependentDifferences)
{
	for (const Vec3D position :
	     {Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 0.0, 1.0), Vec3D(0.01, -0.002, 0.004), Vec3D(5.0, 2.0, -0.5)}) {
		ocf::ForceContext ctx;
		ctx.R = position;
		ctx.V = Vec3D(0.004, 0.018, -0.003);
		ctx.sunGm = Constants::GM_SUN * 1.002;
		checkForcePartials<ocf::SunJ2>(ctx, position.magnitude() * 1e-5);
		Vec3D acceleration;
		Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero(), dadv = Eigen::Matrix3d::Zero();
		ocf::SunJ2{}.evaluate(ctx, acceleration, dadr, dadv);
		EXPECT_NEAR(dadr.trace(), 0.0, dadr.norm() * 1e-14);
		EXPECT_TRUE(dadr.isApprox(dadr.transpose(), 1e-14));
		EXPECT_TRUE(dadv.isZero(0.0));
	}
}

TEST(ForceField, EmptyFieldProducesZeroAcceleration)
{
	ocf::ForceField<> field;
	static_assert(field.empty());
	static_assert(field.size() == 0);

	Vec3D a(1.0, 1.0, 1.0);
	Eigen::Matrix3d dadr = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d dadv = Eigen::Matrix3d::Identity();

	ocf::ForceContext ctx{Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 0.02, 0.0)};
	ctx.setEpochMjd(2459266.348 - 2400000.5);
	field.evaluateZeroed(ctx, a, dadr, dadv);

	EXPECT_DOUBLE_EQ(a.x, 0.0);
	EXPECT_DOUBLE_EQ(a.y, 0.0);
	EXPECT_DOUBLE_EQ(a.z, 0.0);
	EXPECT_TRUE(dadr.isZero());
	EXPECT_TRUE(dadv.isZero());
}

TEST(SunGravity, AccelerationMatchesNewtonAtOneAU)
{
	ocf::SunGravity sun;
	ocf::ForceContext ctx;
	ctx.R = Vec3D(1.0, 0.0, 0.0);
	ctx.V = Vec3D(0.0, 0.0, 0.0);
	ctx.setEpochMjd(2459266.348 - 2400000.5);

	Vec3D a(0.0, 0.0, 0.0);
	Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dadv = Eigen::Matrix3d::Zero();

	sun.evaluate(ctx, a, dadr, dadv);

	EXPECT_NEAR(a.x, -Constants::GM_SUN, 1e-15);
	EXPECT_NEAR(a.y, 0.0, 1e-20);
	EXPECT_NEAR(a.z, 0.0, 1e-20);
}

TEST(SunGravity, JacobianMatchesAnalyticAtOneAU)
{
	ocf::SunGravity sun;
	ocf::ForceContext ctx;
	ctx.R = Vec3D(1.0, 0.0, 0.0);

	Vec3D a(0.0, 0.0, 0.0);
	Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dadv = Eigen::Matrix3d::Zero();

	sun.evaluate(ctx, a, dadr, dadv);

	EXPECT_NEAR(dadr(0, 0), 2.0 * Constants::GM_SUN, 1e-12);
	EXPECT_NEAR(dadr(1, 1), -Constants::GM_SUN, 1e-12);
	EXPECT_NEAR(dadr(2, 2), -Constants::GM_SUN, 1e-12);
	EXPECT_NEAR(dadr(0, 1), 0.0, 1e-15);
	EXPECT_NEAR(dadr(0, 2), 0.0, 1e-15);
	EXPECT_NEAR(dadr(1, 2), 0.0, 1e-15);
}
