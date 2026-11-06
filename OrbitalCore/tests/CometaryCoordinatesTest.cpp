#include <gtest/gtest.h>
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/Constants.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

using namespace orbital_core;
using namespace orbital_core::cometary;

namespace {
	constexpr double radiansPerDegree = std::numbers::pi / 180.0;

	Vector6 cotSeedK14UR1N()
	{
		Vector6 elements;
		elements << 1.09392887734970E+01, 0.998177173568370, 95.5383742001458, 190.0663927572708, 326.3748746309441,
		    264.7948222869081;
		elements.tail<4>() *= radiansPerDegree;
		return elements;
	}

	Vector6 carSeedK14UR1N()
	{
		Vector6 state;
		state << 1.5154308780457381E+01, 8.5491346989789696E-01, -1.8636168138824399E+01, -4.8372148116785155E-03,
		    -7.8508657258669805E-04, 7.4761153441349395E-04;
		return state;
	}

	Matrix6 symmetricCovariance(const std::array<double, 21>& upper)
	{
		Matrix6 covariance;
		int index = 0;
		for (int row = 0; row < 6; ++row)
			for (int column = row; column < 6; ++column)
				covariance(row, column) = covariance(column, row) = upper[index++];
		return covariance;
	}

	Matrix6 cotSeedCovarianceK14UR1N()
	{
		Matrix6 covariance = symmetricCovariance(
		    {8.945163550510338E-10,  1.043470031827537E-10,  6.791318449874664E-11, 1.298405911994385E-10,
		     -5.901710523430252E-09, 5.889750433018134E-09,  1.961690340725353E-11, 8.277832220365925E-12,
		     1.675999013057592E-11,  -8.029733069266978E-10, 8.110608123704280E-10, 3.600999541767908E-11,
		     7.769930462176083E-11,  -4.428186358322275E-10, 4.512100449237282E-10, 3.307115212313096E-10,
		     -8.437217484869453E-10, 8.783479219674374E-10,  4.091837154454448E-08, -4.085628087232593E-08,
		     4.088568865708540E-08});
		covariance.bottomRows<4>() *= radiansPerDegree;
		covariance.rightCols<4>() *= radiansPerDegree;
		return covariance;
	}

	Matrix6 carSeedCovarianceK14UR1N()
	{
		return symmetricCovariance({3.839164152290086E-10,  2.613854630144268E-11,  -4.822441710395081E-10,
		                            -1.565755247083013E-13, -2.325672716459748E-14, 1.000269729491152E-13,
		                            1.498658976411458E-11,  -3.335777550018451E-11, -1.057470606255184E-14,
		                            -6.738762644527730E-15, 5.797809481778580E-15,  6.442995713233677E-10,
		                            1.977851461750990E-13,  2.937728199217805E-14,  -1.424263053226216E-13,
		                            6.422188995983678E-17,  9.315204478680730E-18,  -4.249078633421696E-17,
		                            3.814013042901125E-18,  -4.161464882816909E-18, 4.704976622872050E-17});
	}

	Vector6 elementDifference(Vector6 a, const Vector6& b)
	{
		a -= b;
		for (int index = 3; index < 6; ++index)
			a(index) = std::remainder(a(index), 2.0 * std::numbers::pi);
		return a;
	}

	void expectStateNear(const Vector6& actual, const Vector6& expected, double relative)
	{
		EXPECT_LE((actual.head<3>() - expected.head<3>()).norm(), relative * expected.head<3>().norm());
		EXPECT_LE((actual.tail<3>() - expected.tail<3>()).norm(), relative * expected.tail<3>().norm());
	}

	void expectCovarianceNear(const Matrix6& actual, const Matrix6& expected, double relative)
	{
		for (int row = 0; row < 6; ++row)
			for (int column = 0; column < 6; ++column) {
				SCOPED_TRACE(::testing::Message() << "row=" << row << " column=" << column);
				EXPECT_NEAR(actual(row, column), expected(row, column),
				            relative * std::sqrt(expected(row, row) * expected(column, column)));
			}
	}
}  // namespace

TEST(CometaryCoordinatesTest, CartesianAndCovarianceRoundTripsAcrossParabolicBoundary)
{
	for (const double e : {0.2, 1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10, 1.8}) {
		for (const double inclination : {0.3, 1.67, 2.6}) {
			for (const double anomaly : {-1.4, 0.0, 1.2, 2.0 * std::numbers::pi - 1.0e-12}) {
				Vector6 elements;
				elements << 10.9, e, inclination, -0.2, 2.0 * std::numbers::pi + 0.4, anomaly;
				SCOPED_TRACE(elements.transpose());
				Matrix6 forward;
				const auto state = toCartesian(elements, &forward);
				ASSERT_TRUE(state);
				const auto restored = fromCartesian(*state);
				ASSERT_TRUE(restored);
				EXPECT_LT(elementDifference(*restored, elements).norm(), 2.0e-12);
				for (int angle = 3; angle < 6; ++angle) {
					EXPECT_GE((*restored)(angle), 0.0);
					EXPECT_LT((*restored)(angle), 2.0 * std::numbers::pi);
				}
				if (e < 1.0)
					EXPECT_LT((*restored)(1), 1.0);
				if (e > 1.0)
					EXPECT_GT((*restored)(1), 1.0);
				const auto again = toCartesian(*restored);
				ASSERT_TRUE(again);
				expectStateNear(*again, *state, 2.0e-13);
				const auto inverse = cartesianToElementsJacobian(*restored);
				ASSERT_TRUE(inverse);
				EXPECT_LT((*inverse * forward - Matrix6::Identity()).norm(), 2.0e-11);
				const Matrix6 cotCovariance = cotSeedCovarianceK14UR1N();
				const Matrix6 carCovariance = forward * cotCovariance * forward.transpose();
				const Matrix6 roundTrip = *inverse * carCovariance * inverse->transpose();
				expectCovarianceNear(roundTrip, cotCovariance, 2.0e-10);
				const Matrix6 sourceCarCovariance = carSeedCovarianceK14UR1N();
				const Matrix6 converted = *inverse * sourceCarCovariance * inverse->transpose();
				const Matrix6 restoredCar = forward * converted * forward.transpose();
				expectCovarianceNear(restoredCar, sourceCarCovariance, 2.0e-10);
			}
		}
	}
}

TEST(CometaryCoordinatesTest, ForwardJacobianMatchesIndependentFivePointDifferences)
{
	for (const double e : {0.2, 1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10, 1.8}) {
		for (const double anomaly : {0.0, -1.4, 1.2}) {
			Vector6 elements = cotSeedK14UR1N();
			elements(1) = e;
			elements(5) = anomaly;
			SCOPED_TRACE(elements.transpose());
			Matrix6 jacobian;
			ASSERT_TRUE(toCartesian(elements, &jacobian));
			for (int column = 0; column < 6; ++column) {
				SCOPED_TRACE(column);
				const double step = 1.0e-4 * (column == 0 ? elements(0) : 1.0);
				const Vector6 offset = step * Vector6::Unit(column);
				const auto plusTwo = toCartesian(elements + 2.0 * offset);
				const auto plus = toCartesian(elements + offset);
				const auto minus = toCartesian(elements - offset);
				const auto minusTwo = toCartesian(elements - 2.0 * offset);
				ASSERT_TRUE(plusTwo && plus && minus && minusTwo);
				const Vector6 finite = (-*plusTwo + 8.0 * *plus - 8.0 * *minus + *minusTwo) / (12.0 * step);
				for (int block = 0; block < 2; ++block) {
					const auto analytic = jacobian.col(column).segment<3>(3 * block);
					const double scale = std::max(analytic.norm(), block == 0 ? elements(0) : 0.01);
					EXPECT_LT((finite.segment<3>(3 * block) - analytic).norm() / scale, 2.0e-9);
				}
			}
		}
	}
}

TEST(CometaryCoordinatesTest, InverseJacobianMatchesWrappedFivePointDifferences)
{
	for (const double e : {0.2, 1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10, 1.8}) {
		Vector6 elements = cotSeedK14UR1N();
		elements(1) = e;
		elements(3) = 0.0;
		elements(4) = 0.0;
		elements(5) = 0.0;
		SCOPED_TRACE(elements.transpose());
		const auto state = toCartesian(elements);
		const auto inverse = cartesianToElementsJacobian(elements);
		ASSERT_TRUE(state && inverse);
		for (int column = 0; column < 6; ++column) {
			SCOPED_TRACE(column);
			const double step = 1.0e-4 * (column < 3 ? state->head<3>().norm() : state->tail<3>().norm());
			const Vector6 offset = step * Vector6::Unit(column);
			const auto plusTwo = fromCartesian(*state + 2.0 * offset);
			const auto plus = fromCartesian(*state + offset);
			const auto minus = fromCartesian(*state - offset);
			const auto minusTwo = fromCartesian(*state - 2.0 * offset);
			ASSERT_TRUE(plusTwo && plus && minus && minusTwo);
			const Vector6 finite =
			    (-elementDifference(*plusTwo, elements) + 8.0 * elementDifference(*plus, elements) -
			     8.0 * elementDifference(*minus, elements) + elementDifference(*minusTwo, elements)) /
			    (12.0 * step);
			EXPECT_LT((finite - inverse->col(column)).norm() / inverse->col(column).norm(), 2.0e-8);
		}
	}
}

TEST(CometaryCoordinatesTest, InverseFactorizationHandlesPhysicalUnitScales)
{
	for (const double q : {1.0e-4, 1.0, 1.0e4}) {
		Vector6 elements = cotSeedK14UR1N();
		elements(0) = q;
		elements(1) = 1.0;
		SCOPED_TRACE(q);
		Matrix6 forward;
		ASSERT_TRUE(toCartesian(elements, &forward));
		const auto inverse = cartesianToElementsJacobian(elements);
		ASSERT_TRUE(inverse);
		Matrix6 product = *inverse * forward;
		product.row(0) /= q;
		product.col(0) *= q;
		EXPECT_LT((product - Matrix6::Identity()).norm(), 2.0e-12);
	}
}

TEST(CometaryCoordinatesTest, MatchesNativeCotSeedK14UR1NAndCovariance)
{
	SCOPED_TRACE("Independent COT seed and Cartesian state reference");
	const Vector6 elements = cotSeedK14UR1N();
	Matrix6 forward;
	const auto state = toCartesian(elements, &forward);
	ASSERT_TRUE(state);
	expectStateNear(*state, carSeedK14UR1N(), 2.0e-13);
	const auto restored = fromCartesian(carSeedK14UR1N());
	ASSERT_TRUE(restored);
	EXPECT_LT(elementDifference(*restored, elements).norm(), 5.0e-13);
	const auto inverse = cartesianToElementsJacobian(elements);
	ASSERT_TRUE(inverse);
	const Matrix6 cotCovariance = cotSeedCovarianceK14UR1N();
	const Matrix6 carCovariance = carSeedCovarianceK14UR1N();
	const Matrix6 convertedCar = forward * cotCovariance * forward.transpose();
	const Matrix6 convertedCot = *inverse * carCovariance * inverse->transpose();
	expectCovarianceNear(convertedCar, carCovariance, 5.0e-10);
	expectCovarianceNear(convertedCot, cotCovariance, 5.0e-10);
}

TEST(CometaryCoordinatesTest, MatchesNativeCotSeedConversionsAtAndAroundExactParabola)
{
	SCOPED_TRACE("Independent CAR reference on both sides of the parabolic limit");
	const std::array<double, 3> eccentricities = {1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10};
	const std::array<Vector6, 3> nativeStates = {
	    (Vector6() << 1.51708918654355E+01, 8.55848986187418E-01, -1.86565613592866E+01, -4.83875921898275E-03,
	     -7.84847137073660E-04, 7.52826767464490E-04)
	        .finished(),
	    (Vector6() << 1.51708918663455E+01, 8.55848986238749E-01, -1.86565613604056E+01, -4.83875921906747E-03,
	     -7.84847137060535E-04, 7.52826767750472E-04)
	        .finished(),
	    (Vector6() << 1.51708918672554E+01, 8.55848986290081E-01, -1.86565613615246E+01, -4.83875921915218E-03,
	     -7.84847137047410E-04, 7.52826768036453E-04)
	        .finished()};
	for (int index = 0; index < 3; ++index) {
		SCOPED_TRACE(index);
		Vector6 elements = cotSeedK14UR1N();
		elements(1) = eccentricities[index];
		const auto state = toCartesian(elements);
		ASSERT_TRUE(state);
		expectStateNear(*state, nativeStates[index], 2.0e-13);
		const auto restored = fromCartesian(nativeStates[index]);
		ASSERT_TRUE(restored);
		EXPECT_LT(elementDifference(*restored, elements).norm(), 5.0e-12);
	}
}

TEST(CometaryCoordinatesTest, InvalidGeometryFailsWithoutWritingJacobian)
{
	auto expectInvalid = [](const Vector6& elements) {
		Matrix6 jacobian = Matrix6::Constant(17.0);
		EXPECT_FALSE(toCartesian(elements));
		EXPECT_FALSE(toCartesian(elements, &jacobian));
		EXPECT_TRUE((jacobian.array() == 17.0).all());
		EXPECT_FALSE(cartesianToElementsJacobian(elements));
	};
	for (int index = 0; index < 6; ++index) {
		for (const double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
		                             -std::numeric_limits<double>::infinity()}) {
			Vector6 elements = cotSeedK14UR1N();
			elements(index) = invalid;
			expectInvalid(elements);
			Vector6 state = carSeedK14UR1N();
			state(index) = invalid;
			EXPECT_FALSE(fromCartesian(state));
		}
	}
	for (const double q : {0.0, -1.0, std::numeric_limits<double>::max()}) {
		Vector6 elements = cotSeedK14UR1N();
		elements(0) = q;
		expectInvalid(elements);
	}
	Vector6 elements = cotSeedK14UR1N();
	elements(1) = -1.0e-10;
	expectInvalid(elements);
	for (const double e : {1.0, 1.8}) {
		elements(1) = e;
		elements(5) = std::numbers::pi;
		expectInvalid(elements);
	}
	EXPECT_FALSE(fromCartesian(Vector6::Zero()));
	Vector6 radial;
	radial << 1.0, 2.0, 3.0, 0.01, 0.02, 0.03;
	EXPECT_FALSE(fromCartesian(radial));
	EXPECT_FALSE(fromCartesian(Vector6::Constant(std::numeric_limits<double>::max())));
}

TEST(CometaryCoordinatesTest, CircularAndEquatorialAnglesFailExplicitlyForJacobians)
{
	for (const double e : {0.0, 1.0e-16, 0.5}) {
		for (const double inclination : {0.0, 1.0e-16, 1.2, std::numbers::pi}) {
			if (e == 0.5 && inclination == 1.2)
				continue;
			Vector6 elements = cotSeedK14UR1N();
			elements(1) = e;
			elements(2) = inclination;
			SCOPED_TRACE(elements.transpose());
			const auto state = toCartesian(elements);
			ASSERT_TRUE(state);
			const auto restored = fromCartesian(*state);
			if (e < 1.0e-14) {
				EXPECT_FALSE(restored);
			} else {
				ASSERT_TRUE(restored);
				EXPECT_DOUBLE_EQ((*restored)(3), 0.0);
				const auto again = toCartesian(*restored);
				ASSERT_TRUE(again);
				expectStateNear(*again, *state, 2.0e-13);
			}
			Matrix6 jacobian = Matrix6::Constant(17.0);
			EXPECT_FALSE(toCartesian(elements, &jacobian));
			EXPECT_TRUE((jacobian.array() == 17.0).all());
			EXPECT_FALSE(cartesianToElementsJacobian(elements));
		}
	}
}

TEST(CometaryCoordinatesTest, EquatorialNoncircularConicsUseZeroNodeConvention)
{
	for (const double e : {0.2, 1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10, 1.8}) {
		for (const double inclination : {0.0, std::numbers::pi}) {
			Vector6 elements = cotSeedK14UR1N();
			elements(1) = e;
			elements(2) = inclination;
			SCOPED_TRACE(elements.transpose());
			const auto state = toCartesian(elements);
			ASSERT_TRUE(state);
			const auto restored = fromCartesian(*state);
			ASSERT_TRUE(restored);
			EXPECT_DOUBLE_EQ((*restored)(3), 0.0);
			EXPECT_NEAR((*restored)(2), inclination, 2.0e-15);
			const double expectedPeri = elements(4) + (inclination == 0.0 ? elements(3) : -elements(3));
			EXPECT_NEAR(std::remainder((*restored)(4) - expectedPeri, 2.0 * std::numbers::pi), 0.0, 2.0e-14);
			const auto again = toCartesian(*restored);
			ASSERT_TRUE(again);
			expectStateNear(*again, *state, 2.0e-13);
			EXPECT_FALSE(cartesianToElementsJacobian(*restored));
		}
	}
}

TEST(CometaryCoordinatesTest, CoordinateSelectionIsNotConicClassification)
{
	static_assert(!usesCometaryElements(0.989999));
	static_assert(usesCometaryElements(0.99));
	static_assert(usesCometaryElements(1.0));
	static_assert(usesCometaryElements(1.2));
	EXPECT_FALSE(usesCometaryElements(std::nextafter(0.99, 0.0)));
	EXPECT_TRUE(usesCometaryElements(std::nextafter(0.99, 1.0)));
}

TEST(CometaryCoordinatesTest, PerihelionTimeIsContinuousAtExactParabola)
{
	for (const double q : {0.07, 1.0, 10.939288773497}) {
		for (const double nu : {-2.5, -1.4, -1.0e-8, 0.0, 1.0e-8, 1.4, 2.5}) {
			SCOPED_TRACE(::testing::Message() << "q=" << q << " nu=" << nu);
			const double d = std::tan(0.5 * nu);
			const double barker = std::sqrt(2.0 * q * q * q) / Constants::GAUSS_K * (d + d * d * d / 3.0);
			const double parabolic = timeSincePerihelion(q, 1.0, nu);
			EXPECT_NEAR(parabolic, barker, 5.0e-15 * std::max(1.0, std::abs(barker)));
			for (const double e : {1.0 - 1.0e-10, std::nextafter(1.0, 0.0), std::nextafter(1.0, 2.0), 1.0 + 1.0e-10}) {
				SCOPED_TRACE(e);
				const double time = timeSincePerihelion(q, e, nu);
				EXPECT_NEAR(time, parabolic, 2.0e-9 * std::max(1.0, std::abs(parabolic)));
				EXPECT_NEAR(timeSincePerihelion(q, e, -nu), -time, 2.0e-14 * std::max(1.0, std::abs(time)));
			}
		}
	}
}

TEST(CometaryCoordinatesTest, PerihelionTimeMatchesEllipticAndHyperbolicKeplerEquations)
{
	const double q = 1.7;
	const double ellipticE = 0.7;
	const double hyperbolicE = 1.7;
	for (const double anomaly : {-2.0, -0.4, 0.4, 2.0}) {
		SCOPED_TRACE(anomaly);
		const double ellipticNu = 2.0 * std::atan2(std::sqrt(1.0 + ellipticE) * std::sin(0.5 * anomaly),
		                                           std::sqrt(1.0 - ellipticE) * std::cos(0.5 * anomaly));
		const double a = q / (1.0 - ellipticE);
		const double ellipticTime =
		    (anomaly - ellipticE * std::sin(anomaly)) * std::sqrt(a * a * a) / Constants::GAUSS_K;
		EXPECT_NEAR(timeSincePerihelion(q, ellipticE, ellipticNu), ellipticTime, 2.0e-12 * std::abs(ellipticTime));
		const double hyperbolicNu =
		    2.0 * std::atan(std::sqrt((hyperbolicE + 1.0) / (hyperbolicE - 1.0)) * std::tanh(0.5 * anomaly));
		const double absA = q / (hyperbolicE - 1.0);
		const double hyperbolicTime =
		    (hyperbolicE * std::sinh(anomaly) - anomaly) * std::sqrt(absA * absA * absA) / Constants::GAUSS_K;
		EXPECT_NEAR(timeSincePerihelion(q, hyperbolicE, hyperbolicNu), hyperbolicTime,
		            2.0e-12 * std::abs(hyperbolicTime));
	}
}

TEST(CometaryCoordinatesTest, PerihelionTimeUsesNearestPassageAndNativeSeedEpochMinusTp)
{
	SCOPED_TRACE("Independent initial COM epoch and perihelion passage reference");
	const Vector6 elements = cotSeedK14UR1N();
	const double nativeTime = 58335.356851380 - 62893.7997488747;
	EXPECT_NEAR(timeSincePerihelion(elements(0), elements(1), elements(5)), nativeTime, 5.0e-8);
	for (const double e : {0.0, 0.7, 1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10, 1.7}) {
		const double expected = timeSincePerihelion(1.7, e, -1.4);
		EXPECT_LT(expected, 0.0);
		EXPECT_DOUBLE_EQ(timeSincePerihelion(1.7, e, 0.0), 0.0);
		for (const int turn : {-3, -1, 1, 3})
			EXPECT_NEAR(timeSincePerihelion(1.7, e, -1.4 + turn * 2.0 * std::numbers::pi), expected,
			            2.0e-13 * std::abs(expected));
	}
	const double a = 1.7 / (1.0 - 0.7);
	const double halfPeriod = std::numbers::pi * std::sqrt(a * a * a) / Constants::GAUSS_K;
	EXPECT_NEAR(timeSincePerihelion(1.7, 0.7, std::numbers::pi), halfPeriod, 2.0e-12 * halfPeriod);
	EXPECT_NEAR(timeSincePerihelion(1.7, 0.7, -std::numbers::pi), -halfPeriod, 2.0e-12 * halfPeriod);
}

TEST(CometaryCoordinatesTest, PerihelionTimeRejectsInvalidGeometry)
{
	for (const double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
	                             -std::numeric_limits<double>::infinity()}) {
		EXPECT_TRUE(std::isnan(timeSincePerihelion(invalid, 1.0, 0.4)));
		EXPECT_TRUE(std::isnan(timeSincePerihelion(1.0, invalid, 0.4)));
		EXPECT_TRUE(std::isnan(timeSincePerihelion(1.0, 1.0, invalid)));
	}
	EXPECT_TRUE(std::isnan(timeSincePerihelion(0.0, 1.0, 0.4)));
	EXPECT_TRUE(std::isnan(timeSincePerihelion(-1.0, 1.0, 0.4)));
	EXPECT_TRUE(std::isnan(timeSincePerihelion(1.0, -1.0, 0.4)));
	EXPECT_TRUE(std::isnan(timeSincePerihelion(1.0, 1.0, std::numbers::pi)));
	EXPECT_TRUE(std::isnan(timeSincePerihelion(1.0, 1.7, std::numbers::pi)));
	EXPECT_TRUE(std::isnan(timeSincePerihelion(std::numeric_limits<double>::max(), 1.0, 0.4)));
}
