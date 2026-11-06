#include <gtest/gtest.h>
#include "orbital_core/hazard/EquinoctialCoordinates.h"
#include "hazard/LovCoordinateChart.h"
#include "orbital_core/Constants.h"
#include <Eigen/Geometry>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numbers>

using namespace orbital_core;
using namespace orbital_core::hazard;
using namespace orbital_core::hazard::equinoctial;

namespace {
	Vector6 vvState()
	{
		Vector6 x;
		x << 0.62434444840226033, 0.83392068676579079, 0.0034709928450474057, -0.012853569243758767,
		    0.011750580852292960, 0.000093051552130687213;
		return x;
	}

	Vector6 elementDifference(Vector6 a, const Vector6& b)
	{
		a -= b;
		a(5) = std::remainder(a(5), 2.0 * std::numbers::pi);
		return a;
	}
}  // namespace

TEST(EquinoctialCoordinatesTest, RoundTripsIncludingCircularEquatorialAndLongitudeWrap)
{
	for (const double e : {0.0, 1.0e-10, 0.2, 0.95, 0.999}) {
		for (const double p : {0.0, 1.0e-10, 0.2, 4.0}) {
			for (const double longitude : {-1.0e-12, 0.0, 1.3, 2.0 * std::numbers::pi + 1.0e-12}) {
				Vector6 elements;
				elements << 1.4, 0.6 * e, 0.8 * e, p, -0.7 * p, longitude;
				SCOPED_TRACE(elements.transpose());
				const auto state = toCartesian(elements);
				ASSERT_TRUE(state);
				const auto restored = fromCartesian(*state);
				ASSERT_TRUE(restored);
				EXPECT_LT(elementDifference(*restored, elements).norm(), 3.0e-10);
				const auto again = toCartesian(*restored);
				ASSERT_TRUE(again);
				EXPECT_LT((*again - *state).norm(), 3.0e-10);
			}
		}
	}
}

TEST(EquinoctialCoordinatesTest, AnalyticJacobianMatchesIndependentFivePointDifference)
{
	for (const double e : {0.0, 1.0e-8, 0.2, 0.85}) {
		Vector6 elements;
		elements << 1.4, 0.6 * e, 0.8 * e, 0.001, -0.002, 1.3;
		Matrix6 j;
		ASSERT_TRUE(toCartesian(elements, &j));
		for (int column = 0; column < 6; ++column) {
			const double step = 1.0e-4;
			const Vector6 offset = step * Vector6::Unit(column);
			const Vector6 finite = (-*toCartesian(elements + 2.0 * offset) + 8.0 * *toCartesian(elements + offset) -
			                        8.0 * *toCartesian(elements - offset) + *toCartesian(elements - 2.0 * offset)) /
			                       (12.0 * step);
			EXPECT_LT((finite - j.col(column)).norm() / j.col(column).norm(), 1.0e-9);
		}
		const auto inverse = cartesianToElementsJacobian(elements);
		ASSERT_TRUE(inverse);
		EXPECT_LT((*inverse * j - Matrix6::Identity()).norm(), 2.0e-12);
	}
}

TEST(EquinoctialCoordinatesTest, CovariancePoleDerivativesStayStableUnderMillimetreStateChanges)
{
	const Matrix6 covariance = (Vector6() << 7e-12, 2e-12, 3e-14, 5e-15, 9e-15, 2e-17).finished().asDiagonal();
	for (int sample = 0; sample < 11; ++sample) {
		Vector6 state = vvState();
		state(0) += sample * 1.0e-15;
		const Eigen::Vector3d r = state.head<3>(), v = state.tail<3>();
		const Eigen::Vector3d angular = r.cross(v);
		const double norm = angular.norm(), divisor = norm + angular.z();
		Eigen::Matrix<double, 3, 6> dh;
		for (int column = 0; column < 3; ++column) {
			const Eigen::Vector3d unit = Eigen::Vector3d::Unit(column);
			dh.col(column) = unit.cross(v);
			dh.col(column + 3) = r.cross(unit);
		}
		const Eigen::RowVector<double, 6> dd = angular.transpose() * dh / norm + dh.row(2);
		Eigen::Matrix<double, 2, 6> direct;
		direct.row(0) = dh.row(0) / divisor - angular.x() * dd / (divisor * divisor);
		direct.row(1) = -dh.row(1) / divisor + angular.y() * dd / (divisor * divisor);
		const auto elements = fromCartesian(state);
		ASSERT_TRUE(elements);
		const auto j = cartesianToElementsJacobian(*elements);
		ASSERT_TRUE(j);
		const Eigen::Matrix2d expected = direct * covariance * direct.transpose();
		const Eigen::Matrix2d actual = j->middleRows<2>(3) * covariance * j->middleRows<2>(3).transpose();
		for (int i = 0; i < 2; ++i)
			EXPECT_NEAR(std::sqrt(actual(i, i) / expected(i, i)), 1.0, 2.0e-12);
	}
}

TEST(EquinoctialCoordinatesTest, NonGravIdentityAndCrossCovariancesArePreserved)
{
	NonGravModel ng;
	ng.setModel(NgModel::Yarkovsky);
	ng.setEnabled(true);
	ng.setParams(0.0, -3.0e-14, 0.0);
	const auto x = vvState();
	const auto orbit = make_orbit_from_states_vectors(
	    StateVectors(Vec3D(x(0), x(1), x(2)), Vec3D(x(3), x(4), x(5)) / Constants::GAUSS_K), 2459169.5);
	ASSERT_TRUE(orbit);
	Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(7, 7) * 1.0e-12;
	covariance(6, 6) = 1.0e-28;
	covariance(0, 6) = covariance(6, 0) = 1.0e-21;
	const auto chart = makeEquinoctialLovChart(*orbit, covariance, ng);
	ASSERT_TRUE(chart.valid);
	const auto state = physicalParametersFromChart(chart, chart.nominalCoordinates);
	const auto j = physicalJacobianFromChart(chart, chart.nominalCoordinates);
	ASSERT_TRUE(state && j);
	EXPECT_DOUBLE_EQ((*state)(6), ng.A2());
	EXPECT_DOUBLE_EQ((*j)(6, 6), 1.0);
	EXPECT_TRUE(j->topRightCorner(6, 1).isZero(0.0));
	EXPECT_TRUE(j->bottomLeftCorner(1, 6).isZero(0.0));
	const Eigen::MatrixXd restored = *j * chart.coordinateCovariance * j->transpose();
	EXPECT_LT((restored - covariance).norm() / covariance.norm(), 1.0e-10);
	EXPECT_DOUBLE_EQ(chart.coordinateCovariance(6, 6), covariance(6, 6));
	EXPECT_GT(chart.coordinateCovariance.topRightCorner(6, 1).norm(), 0.0);
	EXPECT_TRUE(chart.coordinateCovariance.isApprox(chart.coordinateCovariance.transpose(), 1.0e-15));
}

TEST(EquinoctialCoordinatesTest, InvalidAndSingularStatesFailExplicitly)
{
	EXPECT_FALSE(fromCartesian(Vector6::Zero()));
	Vector6 elements = Vector6::Zero();
	EXPECT_FALSE(toCartesian(elements));
	elements(0) = 1.0;
	elements(1) = 1.0;
	EXPECT_FALSE(toCartesian(elements));
	elements(1) = 0.0;
	elements(5) = std::numeric_limits<double>::quiet_NaN();
	EXPECT_FALSE(toCartesian(elements));
	Vector6 retrograde;
	retrograde << 1.0, 0.0, 0.0, 0.0, -Constants::GAUSS_K, 0.0;
	EXPECT_FALSE(fromCartesian(retrograde));
	Vector6 hyperbolic;
	hyperbolic << 1.0, 0.0, 0.0, 0.0, 0.04, 0.0;
	EXPECT_FALSE(fromCartesian(hyperbolic));
}
