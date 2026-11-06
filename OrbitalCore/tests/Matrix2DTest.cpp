#include <gtest/gtest.h>
#include "orbital_core/Matrix2D.h"
#include <cmath>

TEST(Matrix2D, DefaultConstructorZeros)
{
	Matrix2D<double, 3, 3> m;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			EXPECT_EQ(m(i, j), 0.0);
}

TEST(Matrix2D, InitializerListConstructor)
{
	Matrix2D<double, 2, 2> m({{1.0, 2.0}, {3.0, 4.0}});
	EXPECT_EQ(m(0, 0), 1.0);
	EXPECT_EQ(m(0, 1), 2.0);
	EXPECT_EQ(m(1, 0), 3.0);
	EXPECT_EQ(m(1, 1), 4.0);
}

TEST(Matrix2D, IndexOutOfRange)
{
	Matrix2D<double, 2, 2> m;
	EXPECT_THROW(m(2, 0), std::out_of_range);
	EXPECT_THROW(m(0, 2), std::out_of_range);
}

TEST(Matrix2D, ScalarMultiplication)
{
	Matrix2D<double, 2, 2> m({{1.0, 2.0}, {3.0, 4.0}});
	auto r = m * 3.0;
	EXPECT_EQ(r(0, 0), 3.0);
	EXPECT_EQ(r(0, 1), 6.0);
	EXPECT_EQ(r(1, 0), 9.0);
	EXPECT_EQ(r(1, 1), 12.0);
}

TEST(Matrix2D, Addition)
{
	Matrix2D<double, 2, 2> a({{1.0, 2.0}, {3.0, 4.0}});
	Matrix2D<double, 2, 2> b({{5.0, 6.0}, {7.0, 8.0}});
	auto c = a + b;
	EXPECT_EQ(c(0, 0), 6.0);
	EXPECT_EQ(c(0, 1), 8.0);
	EXPECT_EQ(c(1, 0), 10.0);
	EXPECT_EQ(c(1, 1), 12.0);
}

TEST(Matrix2D, MatrixMultiplication)
{
	Matrix2D<double, 2, 2> a({{1.0, 2.0}, {3.0, 4.0}});
	Matrix2D<double, 2, 2> b({{5.0, 6.0}, {7.0, 8.0}});
	auto c = a * b;
	EXPECT_NEAR(c(0, 0), 19.0, 1e-12);
	EXPECT_NEAR(c(0, 1), 22.0, 1e-12);
	EXPECT_NEAR(c(1, 0), 43.0, 1e-12);
	EXPECT_NEAR(c(1, 1), 50.0, 1e-12);
}

TEST(Matrix2D, MatrixVectorMultiplication)
{
	Matrix2D<double, 3, 3> m({{1.0, 0.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 1.0, 0.0}});
	Matrix2D<double, 3, 1> v({{1.0}, {2.0}, {3.0}});
	auto r = m * v;
	EXPECT_NEAR(r(0, 0), 1.0, 1e-15);
	EXPECT_NEAR(r(1, 0), -3.0, 1e-15);
	EXPECT_NEAR(r(2, 0), 2.0, 1e-15);
}

TEST(Matrix2D, Transpose)
{
	Matrix2D<double, 2, 3> m({{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
	auto t = m.transpose();
	EXPECT_EQ(t(0, 0), 1.0);
	EXPECT_EQ(t(1, 0), 2.0);
	EXPECT_EQ(t(2, 0), 3.0);
	EXPECT_EQ(t(0, 1), 4.0);
	EXPECT_EQ(t(1, 1), 5.0);
	EXPECT_EQ(t(2, 1), 6.0);
}

TEST(Matrix2D, IdentityInverse)
{
	Matrix2D<double, 3, 3> I({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
	auto inv = I.inverse();
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			EXPECT_NEAR(inv(i, j), (i == j) ? 1.0 : 0.0, 1e-12);
}

TEST(Matrix2D, InverseRoundtrip)
{
	Matrix2D<double, 3, 3> m({{2.0, 1.0, 0.0}, {1.0, 3.0, 1.0}, {0.0, 1.0, 2.0}});
	auto inv = m.inverse();
	auto product = m * inv;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			EXPECT_NEAR(product(i, j), (i == j) ? 1.0 : 0.0, 1e-10);
}

TEST(Matrix2D, RotationMatrixInverseEqualsTranspose)
{
	double angle = 0.7;
	Matrix2D<double, 3, 3> Rx(
	    {{1.0, 0.0, 0.0}, {0.0, std::cos(angle), std::sin(angle)}, {0.0, -std::sin(angle), std::cos(angle)}});
	auto inv = Rx.inverse();
	auto trans = Rx.transpose();
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			EXPECT_NEAR(inv(i, j), trans(i, j), 1e-12);
}

TEST(Matrix2D, SingularMatrixThrows)
{
	Matrix2D<double, 2, 2> m({{1.0, 2.0}, {2.0, 4.0}});
	EXPECT_THROW(m.inverse(), std::runtime_error);
}

TEST(Matrix2D, Norm)
{
	Matrix2D<double, 2, 2> m({{3.0, 0.0}, {4.0, 0.0}});
	EXPECT_NEAR(m.norm(), 5.0, 1e-12);
}

TEST(Matrix2D, NormalizedUnitNorm)
{
	Matrix2D<double, 2, 2> m({{3.0, 0.0}, {4.0, 0.0}});
	auto n = m.normalized();
	EXPECT_NEAR(n.norm(), 1.0, 1e-12);
}

TEST(Matrix2D, NormalizeZeroThrows)
{
	Matrix2D<double, 2, 2> m;
	EXPECT_THROW(m.normalized(), std::runtime_error);
}
