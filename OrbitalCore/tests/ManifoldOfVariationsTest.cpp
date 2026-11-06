#include <gtest/gtest.h>
#include "hazard/mov/MovCoordinates.h"
#include <Eigen/LU>

namespace {
	namespace mov = orbital_core::hazard::mov;

	TEST(ManifoldOfVariations, ChartRoundTripAndIndependentDerivatives)
	{
		mov::Chart chart{60000, [](double t) {
			                 const double dt = t - 60000;
			                 return NBodyPropagator::StateVel{{1 + .003 * dt, .01 * dt, .002 * dt}, {.003, .01, .002}};
		                 }};
		mov::Vector6 att;
		att << 2.4, .3, .004, -.003, .02, -.007;
		EXPECT_TRUE(chart.inverse(chart.state(att)).isApprox(att, 1e-12));
		const auto jacobian = chart.derivative(att);
		for (int col = 0; col < 6; ++col) {
			const double h = 1e-5;
			auto plus = att, minus = att;
			plus(col) += h;
			minus(col) -= h;
			const mov::Vector6 numeric = (chart.state(plus) - chart.state(minus)) / (2 * h);
			EXPECT_LT((numeric - jacobian.col(col)).norm() / jacobian.col(col).norm(), 2e-7) << col;
		}
	}

	TEST(ManifoldOfVariations, AdmissibleRegionMatchesPhysicalEnergy)
	{
		mov::Chart chart{60000, [](double) { return NBodyPropagator::StateVel{{1, 0, 0}, {0, .0172, 0}}; }};
		const Eigen::Vector4d angles(2.1, .1, .02, -.003);
		const auto region = mov::admissibleRegion(chart, angles);
		ASSERT_FALSE(region.ranges.empty());
		for (int i = 0; i < 100; ++i) {
			const double rho = .0001 * std::pow(1e6, i / 99.0);
			for (int j = -20; j <= 20; ++j) {
				const double rate = .005 * j;
				const Eigen::Vector3d r = region.position + rho * region.direction;
				const Eigen::Vector3d v = region.velocity + rate * region.direction + rho * region.angularVelocity;
				const double energy = v.squaredNorm() / 2 - Constants::GM_SUN / r.norm();
				EXPECT_EQ(region.contains(rho, rate), energy <= -Constants::GM_SUN / 200) << rho << ' ' << rate;
			}
		}
	}

	TEST(ManifoldOfVariations, VolumeElementUsesRegressionAndNotCovarianceInflation)
	{
		Eigen::MatrixXd design = Eigen::MatrixXd::Zero(8, 6);
		design.topLeftCorner<4, 4>() = Eigen::Matrix4d::Identity();
		design(0, 4) = 3;
		design(1, 5) = 4;
		design(4, 4) = 2;
		design(5, 5) = 1;
		EXPECT_NEAR(mov::manifoldJacobian(design), std::sqrt(170.0), 1e-12);
		EXPECT_NEAR(mov::manifoldJacobian(2 * design), std::sqrt(170.0), 1e-12);
	}

	TEST(ManifoldOfVariations, AcceleratingObserverRangeDerivative)
	{
		mov::Chart chart{
		    60000, [](double t) {
			    const double dt = t - 60000, omega = 6.28;
			    return NBodyPropagator::StateVel{
			        {1 + .00004 * std::cos(omega * dt), .017 * dt + .00004 * std::sin(omega * dt), 0},
			        {-.00004 * omega * std::sin(omega * dt), .017 + .00004 * omega * std::cos(omega * dt), 0}};
		    }};
		mov::Vector6 att;
		att << -3.1, .7, .004, -.002, .1, -.003;
		EXPECT_TRUE(chart.inverse(chart.state(att)).isApprox(att, 1e-12));
		const auto jac = chart.derivative(att);
		auto plus = att, minus = att;
		plus(4) += 1e-4;
		minus(4) -= 1e-4;
		const mov::Vector6 numeric = (chart.state(plus) - chart.state(minus)) / 2e-4;
		EXPECT_LT((numeric - jac.col(4)).norm(), 2e-9);
	}

	TEST(ManifoldOfVariations, RankDeficientVolumeRejected)
	{
		EXPECT_THROW(mov::manifoldJacobian(Eigen::MatrixXd::Zero(8, 6)), std::runtime_error);
		EXPECT_THROW(mov::manifoldJacobian(Eigen::MatrixXd::Identity(8, 5)), std::runtime_error);
	}

	TEST(ManifoldOfVariations, RegionBoundsContainIndependentEnergyScan)
	{
		mov::Chart chart{60000, [](double) { return NBodyPropagator::StateVel{{1, 0, 0}, {0, .0172, 0}}; }};
		for (double angle : {.2, 1.1, 2.1, 3.1})
			for (double speed : {-.03, -.003, .003, .03}) {
				const auto region = mov::admissibleRegion(chart, Eigen::Vector4d(angle, .1, speed, .001));
				for (int i = 0; i < 250; ++i) {
					const double rho = .001 * std::pow(1e5, i / 249.0);
					const double disc = region.discriminant(rho);
					if (disc < 0)
						continue;
					const double mid = -region.direction.dot(region.velocity);
					EXPECT_TRUE(region.contains(rho, mid)) << angle << ' ' << speed << ' ' << rho;
					EXPECT_LE(region.rateMin, mid - std::sqrt(disc) + 1e-12);
					EXPECT_GE(region.rateMax, mid + std::sqrt(disc) - 1e-12);
				}
			}
	}
}  // namespace
