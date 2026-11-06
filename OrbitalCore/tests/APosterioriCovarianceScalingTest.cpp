#include <gtest/gtest.h>

#include "fit/APosterioriCovarianceScaling.h"

#include <limits>

TEST(APosterioriCovarianceScalingTest, FloorsResidualRmsAtOne)
{
	EXPECT_DOUBLE_EQ(orbital_core::fit::aPosterioriVarianceFactor(50.0, 100, 6), 100.0 / 94.0);
	EXPECT_DOUBLE_EQ(orbital_core::fit::aPosterioriVarianceFactor(50.0, 156, 6), 1.04);
}

TEST(APosterioriCovarianceScalingTest, UsesResidualRmsWhenAboveOne)
{
	EXPECT_DOUBLE_EQ(orbital_core::fit::aPosterioriVarianceFactor(150.0, 100, 6), 150.0 / 94.0);
}

TEST(APosterioriCovarianceScalingTest, InvalidInputsUseUnitScale)
{
	EXPECT_DOUBLE_EQ(orbital_core::fit::aPosterioriVarianceFactor(-1.0, 100, 6), 1.0);
	EXPECT_DOUBLE_EQ(orbital_core::fit::aPosterioriVarianceFactor(50.0, 6, 6), 1.0);
	EXPECT_DOUBLE_EQ(orbital_core::fit::aPosterioriVarianceFactor(std::numeric_limits<double>::quiet_NaN(), 100, 6),
	                 1.0);
}
