#include <gtest/gtest.h>

#include "hazard/LovCoordinateChart.h"

#include "orbital_core/Constants.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/LineOfVariationImpact.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <cmath>
#include <limits>

using namespace orbital_core;
using namespace orbital_core::hazard;

namespace {
	Orbit_p makeStateOrbit(const Vec3D& rAu, const Vec3D& vPhysAuPerDay, double epoch)
	{
		return make_orbit_from_states_vectors(StateVectors(rAu, vPhysAuPerDay * (1.0 / Constants::GAUSS_K)), epoch);
	}

	const Vec3D kR(1.2, 0.10, 0.03);
	const Vec3D kVphys(-0.0015, 0.0150, 0.0006);
	constexpr double kEpoch = 2459800.5;
}  // namespace

TEST(LovCoordinateChartTest, EquinoctialUnitsFollowCurrentNode)
{
	Eigen::VectorXd coordinates = Eigen::VectorXd::Zero(7);
	coordinates(0) = 1.1174239505458232;
	const auto units = lovCoordinateUnits(LovCoordinateSystem::Equinoctial, coordinates, true);
	ASSERT_EQ(units.size(), 7);
	EXPECT_DOUBLE_EQ(units(0), coordinates(0));
	EXPECT_DOUBLE_EQ(units(5), Constants::TWOPI);
	EXPECT_TRUE(units.segment<4>(1).isOnes(0.0));
	EXPECT_DOUBLE_EQ(units(6), 1.0);
	coordinates(0) = 1.3;
	EXPECT_DOUBLE_EQ(lovCoordinateUnits(LovCoordinateSystem::Equinoctial, coordinates, true)(0), 1.3);
}

TEST(LovCoordinateChartTest, KeplerianUnitsAndUnscaledMode)
{
	EXPECT_TRUE(LineOfVariationImpactConfig{}.coordinateScaling);
	Eigen::VectorXd coordinates = Eigen::VectorXd::Ones(6);
	const auto units = lovCoordinateUnits(LovCoordinateSystem::Keplerian, coordinates, true);
	ASSERT_EQ(units.size(), 6);
	EXPECT_DOUBLE_EQ(units(0), 1.0e-3);
	EXPECT_DOUBLE_EQ(units(1), 1.0);
	EXPECT_DOUBLE_EQ(units(2), Constants::PI);
	EXPECT_TRUE(units.tail<3>().isConstant(Constants::TWOPI, 0.0));
	for (const auto system : {LovCoordinateSystem::Equinoctial, LovCoordinateSystem::Keplerian}) {
		EXPECT_TRUE(lovCoordinateUnits(system, coordinates, false).isOnes(0.0));
	}
}

TEST(LovCoordinateChartTest, CoordinateUnitsRejectInvalidScale)
{
	using System = LovCoordinateSystem;
	EXPECT_EQ(lovCoordinateUnits(System::Equinoctial, Eigen::VectorXd::Ones(5), true).size(), 0);
	Eigen::VectorXd coordinates = Eigen::VectorXd::Ones(6);
	for (const double value :
	     {0.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
		coordinates(0) = value;
		EXPECT_EQ(lovCoordinateUnits(System::Equinoctial, coordinates, true).size(), 0);
	}
}

TEST(LovCoordinateChartTest, AutomaticCometaryChartCrossesParabolicBoundary)
{
	cometary::Vector6 cot;
	cot << 10.9, 1.0 - 1.0e-10, 1.6, 3.2, 5.7, -1.2;
	const auto state = cometary::toCartesian(cot);
	ASSERT_TRUE(state);
	auto orbit = makeStateOrbit(Vec3D((*state)(0), (*state)(1), (*state)(2)),
	                            Vec3D((*state)(3), (*state)(4), (*state)(5)), kEpoch);
	ASSERT_TRUE(orbit);
	Eigen::Matrix<double, 6, 1> sigma;
	sigma << 1.0e-7, 1.0e-7, 1.0e-7, Constants::GAUSS_K * 1.0e-7, Constants::GAUSS_K * 1.0e-7,
	    Constants::GAUSS_K * 1.0e-7;
	const FitParameterCovariance covariance = sigma.array().square().matrix().asDiagonal();
	const auto chart = makeAutomaticLovChart(*orbit, covariance, NonGravModel{});
	ASSERT_TRUE(chart.valid);
	EXPECT_EQ(chart.coordinateSystem, LovCoordinateSystem::Cometary);
	const auto jac = physicalJacobianFromChart(chart, chart.nominalCoordinates);
	ASSERT_TRUE(jac);
	const Eigen::MatrixXd recovered = *jac * chart.coordinateCovariance * jac->transpose();
	const Eigen::MatrixXd normalizedError =
	    sigma.cwiseInverse().asDiagonal() * (recovered - covariance) * sigma.cwiseInverse().asDiagonal();
	EXPECT_LT(normalizedError.cwiseAbs().maxCoeff(), 1.0e-10);
	for (const double e : {1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10}) {
		Eigen::VectorXd node = chart.nominalCoordinates;
		node(1) = e;
		const auto physical = physicalParametersFromChart(chart, node);
		ASSERT_TRUE(physical);
		EXPECT_TRUE(physical->allFinite());
		const auto next = cometary::fromCartesian(physical->head<6>());
		ASSERT_TRUE(next);
		EXPECT_NEAR((*next)(1), e, 2.0e-15);
	}
	const auto units = lovCoordinateUnits(chart.coordinateSystem, chart.nominalCoordinates, true);
	EXPECT_DOUBLE_EQ(units(0), chart.nominalCoordinates(0));
	EXPECT_DOUBLE_EQ(units(1), 1.0);
	EXPECT_DOUBLE_EQ(units(2), Constants::PI);
	EXPECT_TRUE(units.tail<3>().isConstant(Constants::TWOPI, 0.0));
}

TEST(LovCoordinateChartTest, AutomaticChartPreservesOrdinaryEquinoctialCoordinates)
{
	auto orbit = makeStateOrbit(kR, kVphys, kEpoch);
	ASSERT_TRUE(orbit);
	const FitParameterCovariance covariance = FitParameterCovariance::Identity(6, 6) * 1.0e-12;
	const auto explicitChart = makeEquinoctialLovChart(*orbit, covariance, NonGravModel{});
	const auto automaticChart = makeAutomaticLovChart(*orbit, covariance, NonGravModel{});
	ASSERT_TRUE(automaticChart.valid);
	EXPECT_EQ(automaticChart.coordinateSystem, LovCoordinateSystem::Equinoctial);
	EXPECT_TRUE((automaticChart.nominalCoordinates.array() == explicitChart.nominalCoordinates.array()).all());
	EXPECT_TRUE((automaticChart.coordinateCovariance.array() == explicitChart.coordinateCovariance.array()).all());
}

TEST(LovCoordinateChartTest, PhysicalEquinoctialRoundTrip)
{
	auto orbit = makeStateOrbit(kR, kVphys, kEpoch);
	ASSERT_TRUE(orbit);
	ASSERT_EQ(orbit->getType(), OrbitType::Elliptical);

	NonGravModel ng;
	const FitParameterCovariance cov = FitParameterCovariance::Identity(6, 6) * 1.0e-12;

	const LovCoordinateChart chart = makeEquinoctialLovChart(*orbit, cov, ng);
	ASSERT_TRUE(chart.valid);
	ASSERT_EQ(chart.nominalCoordinates.size(), 6);

	const auto phys = physicalParametersFromChart(chart, chart.nominalCoordinates);
	ASSERT_TRUE(phys.has_value());
	ASSERT_EQ(phys->size(), 6);

	EXPECT_NEAR((*phys)(0), kR.x, 1.0e-9);
	EXPECT_NEAR((*phys)(1), kR.y, 1.0e-9);
	EXPECT_NEAR((*phys)(2), kR.z, 1.0e-9);
	EXPECT_NEAR((*phys)(3), kVphys.x, 1.0e-9);
	EXPECT_NEAR((*phys)(4), kVphys.y, 1.0e-9);
	EXPECT_NEAR((*phys)(5), kVphys.z, 1.0e-9);
}

TEST(LovCoordinateChartTest, PhysicalKeplerianRoundTrip)
{
	auto orbit = makeStateOrbit(kR, kVphys, kEpoch);
	ASSERT_TRUE(orbit);
	ASSERT_EQ(orbit->getType(), OrbitType::Elliptical);

	NonGravModel ng;
	const FitParameterCovariance cov = FitParameterCovariance::Identity(6, 6) * 1.0e-12;

	const LovCoordinateChart chart = makeKeplerianLovChart(*orbit, cov, ng);
	ASSERT_TRUE(chart.valid);
	EXPECT_EQ(chart.coordinateSystem, LovCoordinateSystem::Keplerian);

	const auto phys = physicalParametersFromChart(chart, chart.nominalCoordinates);
	ASSERT_TRUE(phys.has_value());
	ASSERT_EQ(phys->size(), 6);
	EXPECT_NEAR((*phys)(0), kR.x, 1.0e-9);
	EXPECT_NEAR((*phys)(1), kR.y, 1.0e-9);
	EXPECT_NEAR((*phys)(2), kR.z, 1.0e-9);
	EXPECT_NEAR((*phys)(3), kVphys.x, 1.0e-9);
	EXPECT_NEAR((*phys)(4), kVphys.y, 1.0e-9);
	EXPECT_NEAR((*phys)(5), kVphys.z, 1.0e-9);
}

TEST(LovCoordinateChartTest, StepPerSigmaMatchesFiniteDifference)
{
	auto orbit = makeStateOrbit(kR, kVphys, kEpoch);
	ASSERT_TRUE(orbit);

	NonGravModel ng;
	const FitParameterCovariance cov = FitParameterCovariance::Identity(6, 6) * 1.0e-12;
	const LovCoordinateChart chart = makeEquinoctialLovChart(*orbit, cov, ng);
	ASSERT_TRUE(chart.valid);

	Eigen::VectorXd stepDir = Eigen::VectorXd::Zero(6);
	stepDir(0) = 1.0;

	const auto stepPhys = physicalStepPerSigmaFromChart(chart, chart.nominalCoordinates, stepDir);
	ASSERT_TRUE(stepPhys.has_value());
	ASSERT_EQ(stepPhys->size(), 6);

	const double h = 1.0e-5;
	const auto plus = physicalParametersFromChart(chart, chart.nominalCoordinates + h * stepDir);
	const auto minus = physicalParametersFromChart(chart, chart.nominalCoordinates - h * stepDir);
	ASSERT_TRUE(plus.has_value() && minus.has_value());
	const Eigen::VectorXd fd = (*plus - *minus) / (2.0 * h);

	EXPECT_GT(stepPhys->norm(), 0.0);
	for (int i = 0; i < 6; ++i)
		EXPECT_NEAR((*stepPhys)(i), fd(i), 1.0e-6 * (1.0 + std::abs(fd(i))));
}
