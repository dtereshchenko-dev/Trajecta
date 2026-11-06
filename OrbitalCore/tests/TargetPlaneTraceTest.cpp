#include <gtest/gtest.h>

#include "hazard/TargetPlaneTrace.h"

#include <cmath>

using namespace orbital_core::hazard;

TEST(TargetPlaneTraceTest, PreservesMinorAxisOfHighlyStretchedEllipse)
{
	constexpr double major = 1.0e10;
	constexpr double minor = 1.0e-3;
	const double invSqrt2 = 1.0 / std::sqrt(2.0);

	TargetPlaneDifferential differential;
	differential.coordsKm = Eigen::Vector2d(1000.0, 2000.0);
	differential.dCoordsDParams.resize(2, 2);
	differential.dCoordsDParams << major * invSqrt2, -minor * invSqrt2, major * invSqrt2, minor * invSqrt2;
	differential.parameterCovariance = Eigen::Matrix2d::Identity();
	differential.parameterStepPerSigma = Eigen::Vector2d(1.0, 0.0);
	differential.bodyRadiusKm = 6378.0;

	const TargetPlaneTracePoint trace = makeTargetPlaneTracePoint(differential);

	ASSERT_TRUE(trace.valid);
	EXPECT_NEAR(trace.stretchKm, major, major * 1.0e-14);
	EXPECT_NEAR(trace.widthKm, minor, minor * 1.0e-12);
}

TEST(TargetPlaneTraceTest, UsesFocusedImpactCrossSection)
{
	TargetPlaneDifferential differential;
	differential.coordsKm = Eigen::Vector2d(3.0, 0.0);
	differential.dCoordsDParams = Eigen::Matrix2d::Identity();
	differential.parameterCovariance = Eigen::Matrix2d::Identity();
	differential.parameterStepPerSigma = Eigen::Vector2d(1.0, 0.0);
	differential.bodyRadiusKm = 1.0;
	differential.impactCrossSectionBodyRadii = 4.0;

	const TargetPlaneTracePoint trace = makeTargetPlaneTracePoint(differential);

	ASSERT_TRUE(trace.valid);
	EXPECT_DOUBLE_EQ(trace.impactCrossSectionBodyRadii, 4.0);
	EXPECT_TRUE(trace.impact);
}
