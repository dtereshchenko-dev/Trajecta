#include <gtest/gtest.h>
#include "orbital_core/Vec3D.h"
#include <cmath>

TEST(Vec3D, DefaultConstructor)
{
	Vec3D v;
	EXPECT_EQ(v.x, 0.0);
	EXPECT_EQ(v.y, 0.0);
	EXPECT_EQ(v.z, 0.0);
}

TEST(Vec3D, ParameterizedConstructor)
{
	Vec3D v(1.0, 2.0, 3.0);
	EXPECT_EQ(v.x, 1.0);
	EXPECT_EQ(v.y, 2.0);
	EXPECT_EQ(v.z, 3.0);
}

TEST(Vec3D, IsNull)
{
	EXPECT_TRUE(Vec3D().isNull());
	EXPECT_TRUE(Vec3D(0.0, 0.0, 0.0).isNull());
	EXPECT_FALSE(Vec3D(1.0, 0.0, 0.0).isNull());
}

TEST(Vec3D, Magnitude)
{
	EXPECT_NEAR(Vec3D(1.0, 0.0, 0.0).magnitude(), 1.0, 1e-15);
	EXPECT_NEAR(Vec3D(3.0, 4.0, 0.0).magnitude(), 5.0, 1e-15);
	EXPECT_NEAR(Vec3D(1.0, 1.0, 1.0).magnitude(), std::sqrt(3.0), 1e-15);
	EXPECT_NEAR(Vec3D().magnitude(), 0.0, 1e-15);
}

TEST(Vec3D, Normalized)
{
	Vec3D v(3.0, 4.0, 0.0);
	Vec3D n = v.normalized();
	EXPECT_NEAR(n.x, 0.6, 1e-15);
	EXPECT_NEAR(n.y, 0.8, 1e-15);
	EXPECT_NEAR(n.z, 0.0, 1e-15);
	EXPECT_NEAR(n.magnitude(), 1.0, 1e-15);
}

TEST(Vec3D, NormalizedZeroVector)
{
	Vec3D n = Vec3D().normalized();
	EXPECT_EQ(n.x, 0.0);
	EXPECT_EQ(n.y, 0.0);
	EXPECT_EQ(n.z, 0.0);
}

TEST(Vec3D, Addition)
{
	Vec3D a(1.0, 2.0, 3.0);
	Vec3D b(4.0, 5.0, 6.0);
	Vec3D c = a + b;
	EXPECT_EQ(c.x, 5.0);
	EXPECT_EQ(c.y, 7.0);
	EXPECT_EQ(c.z, 9.0);
}

TEST(Vec3D, Subtraction)
{
	Vec3D a(4.0, 5.0, 6.0);
	Vec3D b(1.0, 2.0, 3.0);
	Vec3D c = a - b;
	EXPECT_EQ(c.x, 3.0);
	EXPECT_EQ(c.y, 3.0);
	EXPECT_EQ(c.z, 3.0);
}

TEST(Vec3D, ScalarMultiplication)
{
	Vec3D v(1.0, 2.0, 3.0);
	Vec3D r = v * 2.0;
	EXPECT_EQ(r.x, 2.0);
	EXPECT_EQ(r.y, 4.0);
	EXPECT_EQ(r.z, 6.0);
}

TEST(Vec3D, ComponentMultiplication)
{
	Vec3D a(2.0, 3.0, 4.0);
	Vec3D b(5.0, 6.0, 7.0);
	Vec3D c = a * b;
	EXPECT_EQ(c.x, 10.0);
	EXPECT_EQ(c.y, 18.0);
	EXPECT_EQ(c.z, 28.0);
}

TEST(Vec3D, ScalarDivision)
{
	Vec3D v(4.0, 6.0, 8.0);
	Vec3D r = v / 2.0;
	EXPECT_EQ(r.x, 2.0);
	EXPECT_EQ(r.y, 3.0);
	EXPECT_EQ(r.z, 4.0);
}

TEST(Vec3D, DotProduct)
{
	Vec3D a(1.0, 2.0, 3.0);
	Vec3D b(4.0, 5.0, 6.0);
	EXPECT_NEAR(a.dot(b), 32.0, 1e-15);
}

TEST(Vec3D, DotProductOrthogonal)
{
	Vec3D a(1.0, 0.0, 0.0);
	Vec3D b(0.0, 1.0, 0.0);
	EXPECT_NEAR(a.dot(b), 0.0, 1e-15);
}

TEST(Vec3D, CrossProduct)
{
	Vec3D x(1.0, 0.0, 0.0);
	Vec3D y(0.0, 1.0, 0.0);
	Vec3D z = x.cross(y);
	EXPECT_NEAR(z.x, 0.0, 1e-15);
	EXPECT_NEAR(z.y, 0.0, 1e-15);
	EXPECT_NEAR(z.z, 1.0, 1e-15);
}

TEST(Vec3D, CrossProductAnticommutative)
{
	Vec3D a(1.0, 2.0, 3.0);
	Vec3D b(4.0, 5.0, 6.0);
	Vec3D ab = a.cross(b);
	Vec3D ba = b.cross(a);
	EXPECT_NEAR(ab.x, -ba.x, 1e-15);
	EXPECT_NEAR(ab.y, -ba.y, 1e-15);
	EXPECT_NEAR(ab.z, -ba.z, 1e-15);
}

TEST(Vec3D, CrossProductSelfIsZero)
{
	Vec3D v(1.0, 2.0, 3.0);
	Vec3D r = v.cross(v);
	EXPECT_NEAR(r.magnitude(), 0.0, 1e-14);
}

TEST(Vec3D, IndexOperator)
{
	Vec3D v(10.0, 20.0, 30.0);
	EXPECT_EQ(v[0], 10.0);
	EXPECT_EQ(v[1], 20.0);
	EXPECT_EQ(v[2], 30.0);
}

TEST(Vec3D, IndexOperatorWrite)
{
	Vec3D v;
	v[0] = 1.0;
	v[1] = 2.0;
	v[2] = 3.0;
	EXPECT_EQ(v.x, 1.0);
	EXPECT_EQ(v.y, 2.0);
	EXPECT_EQ(v.z, 3.0);
}

TEST(Vec3D, IndexOutOfRange)
{
	Vec3D v;
	EXPECT_THROW(v[3], std::out_of_range);
}

TEST(Vec3D, Nullify)
{
	Vec3D v(1.0, 2.0, 3.0);
	v.nullify();
	EXPECT_TRUE(v.isNull());
}
