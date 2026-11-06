#include <gtest/gtest.h>
#include "orbital_core/nongrav/NonGravModel.h"
#include "orbital_core/Constants.h"
#include "orbital_core/nongrav/NonGravForceWrapper.h"
#include <cmath>

namespace {
}

TEST(NonGravModelTest, MarsdenCOUsesNormalizedR0FiveLaw)
{
	NonGravModel model;
	model.setModel(NgModel::MarsdenCO);
	EXPECT_NEAR(model.gr(1.0), 1.0, 1e-14);
	EXPECT_NEAR(model.gr(2.0), 0.217214673332965, 1e-12);
	EXPECT_NEAR(model.gr(5.0), 0.00673564805, 1e-11);
}

TEST(NonGravModelTest, StatePartialsMatchIndependentCentralDifferences)
{
	NonGravModel model;
	model.setEnabled(true);
	model.setParams(1.2e-7, -2.3e-8, 3.1e-8);
	const Vec3D r(1.3, 0.4, -0.2);
	for (const auto kind : {NgModel::InverseSquare, NgModel::Yarkovsky, NgModel::MarsdenH2O, NgModel::MarsdenCO,
	                        NgModel::InverseSquareDT}) {
		model.setModel(kind);
		for (double speed : {0.012, 0.035}) {
			const Vec3D v(-0.003, speed, 0.002);
			for (double direction : {-1.0, 1.0}) {
				for (double dt : {-80.0, 0.0, 60.0}) {
					model.setDT(dt);
					const Vec3D velocity = v * direction;
					SCOPED_TRACE(::testing::Message() << int(kind) << " " << speed << " " << direction << " " << dt);
					Eigen::Matrix3d dr, dv;
					model.accelerationPartials(r, velocity, dr, dv);
					ASSERT_TRUE(dr.allFinite());
					ASSERT_TRUE(dv.allFinite());
					for (int col = 0; col < 6; ++col) {
						const double h = col < 3 ? 1e-5 : 1e-7;
						const Vec3D delta(col % 3 == 0 ? h : 0.0, col % 3 == 1 ? h : 0.0, col % 3 == 2 ? h : 0.0);
						const auto plus =
						    model.acceleration(col < 3 ? r + delta : r, col < 3 ? velocity : velocity + delta, 0.0);
						const auto minus =
						    model.acceleration(col < 3 ? r - delta : r, col < 3 ? velocity : velocity - delta, 0.0);
						const auto diff = (plus - minus) / (2.0 * h);
						const Eigen::Vector3d numeric(diff.x, diff.y, diff.z);
						const Eigen::Vector3d analytic =
						    col < 3 ? Eigen::Vector3d(dr.col(col)) : Eigen::Vector3d(dv.col(col - 3));
						EXPECT_LE((analytic - numeric).norm(), 2e-6 * std::max(numeric.norm(), 1e-15)) << col;
					}
				}
			}
		}
	}
}

TEST(NonGravModelTest, StatePartialsCoverTurningPointsAndZeroForce)
{
	NonGravModel model;
	model.setEnabled(true);
	model.setModel(NgModel::MarsdenH2O);
	model.setParams(1e-7, 2e-8, -3e-8);
	model.setDT(30.0);
	const Vec3D r(1.0, 0.0, 0.0);
	for (double speed : {0.012, 0.021, 0.03}) {
		const Vec3D v(0.0, speed, 0.0);
		Eigen::Matrix3d dr, dv;
		model.accelerationPartials(r, v, dr, dv);
		const double h = 1e-5;
		const auto numeric =
		    (model.acceleration(r, v + Vec3D(h, 0, 0), 0.0) * 8.0 -
		     model.acceleration(r, v - Vec3D(h, 0, 0), 0.0) * 8.0 - model.acceleration(r, v + Vec3D(2 * h, 0, 0), 0.0) +
		     model.acceleration(r, v - Vec3D(2 * h, 0, 0), 0.0)) /
		    (12 * h);
		EXPECT_LE((dv.col(0) - Eigen::Vector3d(numeric.x, numeric.y, numeric.z)).norm(),
		          2e-7 * std::max(numeric.magnitude(), 1e-15));
	}
	Eigen::Matrix3d dr, dv;
	model.setParams(0, 0, 0);
	model.accelerationPartials(r, Vec3D(0, 0.012, 0), dr, dv);
	EXPECT_TRUE(dr.isZero(0));
	EXPECT_TRUE(dv.isZero(0));
	model.setParams(1e-7, 0, 0);
	model.setEnabled(false);
	model.accelerationPartials(r, Vec3D(0, 0.012, 0), dr, dv);
	EXPECT_TRUE(dr.isZero(0));
	EXPECT_TRUE(dv.isZero(0));
}

TEST(NonGravModelTest, ForceWrapperAccumulatesStatePartialsIndependentlyOfFitMask)
{
	NonGravModel model;
	model.setEnabled(true);
	model.setModel(NgModel::InverseSquare);
	model.setParams(1e-7, 2e-8, 3e-8);
	model.setFitParameterMask(NgParameterMask::None);
	orbital_core::forces::ForceContext ctx;
	ctx.R = Vec3D(1.0, 0.2, -0.1);
	ctx.V = Vec3D(-0.001, 0.017, 0.002);
	ctx.nonGrav = &model;
	Eigen::Matrix3d dr, dv;
	model.accelerationPartials(ctx.R, ctx.V, dr, dv);
	ASSERT_GT(dr.norm(), 0);
	ASSERT_GT(dv.norm(), 0);
	Eigen::Matrix3d totalR = Eigen::Matrix3d::Identity(), totalV = Eigen::Matrix3d::Identity();
	Vec3D acceleration;
	orbital_core::forces::NonGravForceWrapper{}.evaluate(ctx, acceleration, totalR, totalV);
	EXPECT_EQ(totalR, Eigen::Matrix3d::Identity() + dr);
	EXPECT_EQ(totalV, Eigen::Matrix3d::Identity() + dv);
	ctx.computePartials = false;
	const auto savedR = totalR, savedV = totalV;
	orbital_core::forces::NonGravForceWrapper{}.evaluate(ctx, acceleration, totalR, totalV);
	EXPECT_EQ(totalR, savedR);
	EXPECT_EQ(totalV, savedV);
}

TEST(NonGravModelTest, DisabledReturnsZeroAcceleration)
{
	NonGravModel m;
	m.setModel(NgModel::MarsdenH2O);
	m.setEnabled(false);
	m.setParams(1.0, 0.0, 0.0);

	Vec3D r(2.0, 0.0, 0.0), v(0.0, 0.012, 0.0);
	Vec3D a = m.acceleration(r, v, 0.0);
	EXPECT_DOUBLE_EQ(a.magnitude(), 0.0);
}

TEST(NonGravModelTest, ZeroAValuesReturnZeroAcceleration)
{
	NonGravModel m;
	m.setModel(NgModel::MarsdenH2O);
	m.setEnabled(true);
	m.setParams(0.0, 0.0, 0.0);

	Vec3D r(2.0, 0.0, 0.0), v(0.0, 0.012, 0.0);
	Vec3D a = m.acceleration(r, v, 0.0);
	EXPECT_DOUBLE_EQ(a.magnitude(), 0.0);
}

TEST(NonGravModelTest, NumFitParamsPerModel)
{
	NonGravModel m;
	m.setEnabled(true);

	m.setModel(NgModel::None);
	EXPECT_EQ(m.numFitParams(), 0);
	m.setModel(NgModel::InverseSquare);
	EXPECT_EQ(m.numFitParams(), 3);
	m.setModel(NgModel::MarsdenH2O);
	EXPECT_EQ(m.numFitParams(), 4);
	m.setModel(NgModel::MarsdenCO);
	EXPECT_EQ(m.numFitParams(), 4);

	m.setEnabled(false);
	EXPECT_EQ(m.numFitParams(), 0);
}

TEST(NonGravModelTest, DTShiftsAccelerationMagnitudeConsistentWithKepler)
{
	const double mu = Constants::GM_SUN;
	const double a0 = 1.0;
	const double vCirc = std::sqrt(mu / a0);

	NonGravModel m;
	m.setModel(NgModel::MarsdenH2O);
	m.setEnabled(true);
	m.setParams(1e-7, 0.0, 0.0);

	Vec3D r(a0, 0.0, 0.0);
	Vec3D v(0.0, vCirc, 0.0);

	m.setDT(0.0);
	double mag0 = m.acceleration(r, v, 0.0).magnitude();

	m.setDT(30.0);
	double mag30 = m.acceleration(r, v, 0.0).magnitude();

	m.setDT(-60.0);
	double magM60 = m.acceleration(r, v, 0.0).magnitude();

	const double tolerance = 2e-8 * mag0;
	EXPECT_NEAR(mag0, mag30, tolerance);
	EXPECT_NEAR(mag0, magM60, tolerance);
	EXPECT_GT(mag0, 0.0);
}

TEST(NonGravModelTest, DTOnEllipticalOrbitChangesGConsistently)
{
	const double mu = Constants::GM_SUN;
	const double a = 2.0;
	const double e = 0.5;
	const double rPeri = a * (1.0 - e);
	const double vPeri = std::sqrt(mu * (1.0 + e) / (a * (1.0 - e)));

	NonGravModel m;
	m.setModel(NgModel::MarsdenH2O);
	m.setEnabled(true);
	m.setParams(1e-7, 0.0, 0.0);

	Vec3D r(rPeri, 0.0, 0.0);
	Vec3D v(0.0, vPeri, 0.0);

	m.setDT(0.0);
	double gPeri = m.acceleration(r, v, 0.0).magnitude();

	m.setDT(30.0);
	double gPast = m.acceleration(r, v, 0.0).magnitude();

	EXPECT_LT(gPast, gPeri);
	EXPECT_GT(gPast, 0.0);
}

TEST(NonGravModelTest, AnalyticDTSensitivityMatchesFD)
{
	NonGravModel ng;
	ng.setModel(NgModel::MarsdenH2O);
	ng.setEnabled(true);
	ng.setParams(1e-9, 5e-10, -3e-10);
	ng.setDT(2.0);

	Vec3D r(2.5, 0.3, -0.1);
	Vec3D v(-0.004, 0.014, 0.002);

	Vec3D analytic = ng.sensitivityPartialDT(r, v, 0.0);

	const double h = 1e-3;
	NonGravModel ngp = ng;
	ngp.setDT(2.0 + h);
	NonGravModel ngm = ng;
	ngm.setDT(2.0 - h);
	Vec3D ap = ngp.acceleration(r, v, 0.0);
	Vec3D am = ngm.acceleration(r, v, 0.0);
	Vec3D fd((ap.x - am.x) / (2.0 * h), (ap.y - am.y) / (2.0 * h), (ap.z - am.z) / (2.0 * h));

	EXPECT_NEAR(analytic.x, fd.x, 1e-4 * std::max(1e-14, std::abs(fd.x)))
	    << "∂a_NG/∂DT x mismatch: analytic=" << analytic.x << " fd=" << fd.x;
	EXPECT_NEAR(analytic.y, fd.y, 1e-4 * std::max(1e-14, std::abs(fd.y)))
	    << "∂a_NG/∂DT y mismatch: analytic=" << analytic.y << " fd=" << fd.y;
	EXPECT_NEAR(analytic.z, fd.z, 1e-4 * std::max(1e-14, std::abs(fd.z)))
	    << "∂a_NG/∂DT z mismatch: analytic=" << analytic.z << " fd=" << fd.z;
}

TEST(NonGravModelTest, InverseSquareDTSensitivityMatchesFiniteDifference)
{
	NonGravModel ng;
	ng.setModel(NgModel::InverseSquare);
	ng.setEnabled(true);
	ng.setParams(1e-9, 5e-10, -3e-10);

	Vec3D r(2.5, 0.3, -0.1);
	Vec3D v(-0.004, 0.014, 0.002);

	for (const double dt : {0.0, 2.0, -40.0}) {
		ng.setDT(dt);
		const Vec3D analytic = ng.sensitivityPartialDT(r, v, 0.0);
		constexpr double h = 0.001;
		ng.setDT(dt + h);
		const Vec3D plus = ng.acceleration(r, v, 0.0);
		ng.setDT(dt - h);
		const Vec3D minus = ng.acceleration(r, v, 0.0);
		const Vec3D numerical = (plus - minus) / (2.0 * h);
		EXPECT_LT((analytic - numerical).magnitude(), 1e-5 * numerical.magnitude());
	}
}

TEST(NonGravModelTest, LargeActivityDelayUsesDelayedRadialVelocity)
{
	const Vec3D r(-1.8522186522633379, 0.97983033652375306, 0.25658875037548068);
	const Vec3D v(-0.012109111092259412, -0.0070100643609283192, 0.00062397766597353582);
	NonGravModel model;
	model.setModel(NgModel::MarsdenH2O);
	model.setEnabled(true);
	model.setParams(1.01664960749814e-9, -5.0455905897310078e-11, 2.525142080175053e-10);
	for (double dt : {-100.0, -40.0, 39.985688386194205, 100.0}) {
		SCOPED_TRACE(dt);
		model.setDT(dt);
		const Vec3D analytic = model.sensitivityPartialDT(r, v, 0.0);
		constexpr double h = 0.001;
		model.setDT(dt + h);
		const Vec3D plus = model.acceleration(r, v, 0.0);
		model.setDT(dt - h);
		const Vec3D minus = model.acceleration(r, v, 0.0);
		const Vec3D numerical = (plus - minus) / (2.0 * h);
		EXPECT_LT((analytic - numerical).magnitude(), 1e-7 * numerical.magnitude());
	}
}

TEST(NonGravModelTest, BorisovHyperbolicH2OMatchesIndependentIntegration)
{
	const Vec3D r(-1.8708360327800977, -1.1481752691438272, -2.1095460609000192);
	const Vec3D v(-4.6254083702732560e-4, -1.9822243926802881e-2, -1.2197657542355590e-2);

	struct Reference {
		double direction, dt, g, dg;
	};

	const Reference cases[] = {
	    {1, -92.847801, 1.4245761704674773e-7, 1.3756365907381947e-8},
	    {1, 92.847801, 0.10178466555344309, 0.00081913085608877664},
	    {1, 300, 1.6036032993234853e-7, -1.5525742474508598e-8},
	    {-1, -92.847801, 0.10178466555344309, -0.00081913085608877664},
	    {-1, 92.847801, 1.4245761704674773e-7, -1.3756365907381947e-8},
	    {-1, 300, 2.615347211562683e-14, -1.498140027650049e-15},
	};
	NonGravModel model;
	model.setModel(NgModel::MarsdenH2O);
	model.setEnabled(true);
	model.setParams(1.0, 0.0, 0.0);
	const Vec3D radial = r / r.magnitude();
	for (const auto& reference : cases) {
		SCOPED_TRACE(::testing::Message() << "direction=" << reference.direction << " DT=" << reference.dt);
		model.setDT(reference.dt);
		const auto acceleration = model.acceleration(r, v * reference.direction, 0.0);
		const auto derivative = model.sensitivityPartialDT(r, v * reference.direction, 0.0);
		EXPECT_LT((acceleration - radial * reference.g).magnitude(), reference.g * 3e-11);
		EXPECT_LT((derivative - radial * reference.dg).magnitude(), std::abs(reference.dg) * 5e-8);
	}
}
