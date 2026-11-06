#include <gtest/gtest.h>

#include "fit/PostFitResidual.h"
#include "fit/HystereticRejector.h"
#include "fit/OppositionCoverage.h"

#include <array>
#include <cmath>
#include <limits>

using orbital_core::fit::postFitResidualStatistic;

TEST(PostFitResidualTest, IncludedObservationSubtractsPredictionCovariance)
{
	const auto result = postFitResidualStatistic(Eigen::Vector2d(2.0, 0.0), 0.75 * Eigen::Matrix2d::Identity(), true);
	ASSERT_TRUE(result.valid);
	EXPECT_NEAR(result.chiSquared, 16.0, 1e-13);
}

TEST(PostFitResidualTest, ExcludedObservationAddsPredictionCovariance)
{
	const auto result = postFitResidualStatistic(Eigen::Vector2d(4.0, 0.0), 3.0 * Eigen::Matrix2d::Identity(), false);
	ASSERT_TRUE(result.valid);
	EXPECT_NEAR(result.chiSquared, 4.0, 1e-13);
}

TEST(PostFitResidualTest, IsInvariantToWhiteningAxisRotation)
{
	Eigen::Matrix2d rotation;
	rotation << std::cos(0.7), -std::sin(0.7), std::sin(0.7), std::cos(0.7);
	Eigen::Matrix2d projected;
	projected << 0.45, 0.12, 0.12, 0.25;
	const Eigen::Vector2d residual(0.8, -1.1);
	for (const bool included : {false, true}) {
		const auto original = postFitResidualStatistic(residual, projected, included);
		const auto rotated =
		    postFitResidualStatistic(rotation * residual, rotation * projected * rotation.transpose(), included);
		ASSERT_TRUE(original.valid);
		ASSERT_TRUE(rotated.valid);
		EXPECT_NEAR(original.chiSquared, rotated.chiSquared, 2e-14);
	}
}

TEST(PostFitResidualTest, DoesNotManufactureStatisticFromSingularOrInvalidCovariance)
{
	EXPECT_FALSE(postFitResidualStatistic(Eigen::Vector2d::Ones(), Eigen::Matrix2d::Identity(), true).valid);
	EXPECT_FALSE(postFitResidualStatistic(Eigen::Vector2d::Ones(), 2.0 * Eigen::Matrix2d::Identity(), true).valid);
	EXPECT_FALSE(postFitResidualStatistic(Eigen::Vector2d::Constant(std::numeric_limits<double>::quiet_NaN()),
	                                      Eigen::Matrix2d::Zero(), true)
	                 .valid);
	const auto nearLimit =
	    postFitResidualStatistic(Eigen::Vector2d::Ones(), (1.0 - 1e-8) * Eigen::Matrix2d::Identity(), true);
	ASSERT_TRUE(nearLimit.valid);
	EXPECT_GT(nearLimit.chiSquared, 1e8);
}

TEST(HystereticLeverageTest, CorrectedStatisticsDriveRejectionAndRecovery)
{
	std::vector<ObsNode> observations(6);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459165.5 + i;
	std::vector<std::uint8_t> included{1, 0, 1, 1, 1, 1};
	std::array<double, 6> chiSquared{4.0, 16.0, 0.0, 0.0, 0.0, 0.0};
	const auto uncorrected = Calculations::applyHystereticReject(observations, included, chiSquared);
	EXPECT_EQ(uncorrected.status, Calculations::RejectStatus::NoChange);
	chiSquared[0] =
	    postFitResidualStatistic(Eigen::Vector2d(2.0, 0.0), 0.75 * Eigen::Matrix2d::Identity(), true).chiSquared;
	chiSquared[1] =
	    postFitResidualStatistic(Eigen::Vector2d(4.0, 0.0), 3.0 * Eigen::Matrix2d::Identity(), false).chiSquared;
	const auto corrected = Calculations::applyHystereticReject(observations, included, chiSquared);
	EXPECT_EQ(corrected.status, Calculations::RejectStatus::Changed);
	EXPECT_EQ(corrected.nFlippedOut, 1);
	EXPECT_EQ(corrected.nFlippedIn, 1);
	EXPECT_EQ(included[0], 0);
	EXPECT_EQ(included[1], 1);
}

TEST(HystereticLeverageTest, InvalidFieldLeavesMaskAndCountIntact)
{
	std::vector<ObsNode> observations(6);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459165.5 + i;
	const std::vector<std::uint8_t> original{1, 0, 1, 1, 1, 1};
	auto included = original;
	std::array<double, 6> invalid;
	invalid.fill(std::numeric_limits<double>::quiet_NaN());
	const auto result = Calculations::applyHystereticReject(observations, included, invalid);
	EXPECT_EQ(included, original);
	EXPECT_EQ(result.nKept, 5);
	const auto mismatched = Calculations::applyHystereticReject(observations, included, {});
	EXPECT_EQ(included, original);
	EXPECT_EQ(mismatched.nKept, 5);
}

TEST(HystereticLeverageTest, RecoveryGapDoesNotRaiseTheRejectionThreshold)
{
	std::vector<ObsNode> observations(16);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000.0 + i;
	for (const double sigma : {2.0, 3.0}) {
		SCOPED_TRACE(sigma);
		Calculations::HystereticRejectConfig cfg;
		cfg.rejectChiSquared = sigma * sigma;
		cfg.recoverChiSquared = cfg.rejectChiSquared * cfg.recoveryRatio;
		const double reject = cfg.rejectChiSquared;
		const double recover = cfg.recoverChiSquared;
		EXPECT_DOUBLE_EQ(recover, sigma == 2.0 ? 3.684 : 8.289);
		const double gap = (reject + recover) * 0.5;
		std::vector<std::uint8_t> included{1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1};
		const std::array<double, 16> chiSquared{reject,
		                                        std::nextafter(reject, 12.0),
		                                        gap,
		                                        recover,
		                                        std::nextafter(recover, 0.0),
		                                        std::nextafter(recover, reject),
		                                        gap,
		                                        reject,
		                                        0,
		                                        0,
		                                        0,
		                                        0,
		                                        0,
		                                        0,
		                                        0,
		                                        0};
		const auto result = Calculations::applyHystereticReject(observations, included, chiSquared, cfg);
		EXPECT_EQ(result.nFlippedOut, 1);
		EXPECT_EQ(result.nFlippedIn, 2);
		EXPECT_EQ(included, (std::vector<std::uint8_t>{1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}));
		EXPECT_DOUBLE_EQ(result.appliedRejectChiSquared, reject);
		EXPECT_DOUBLE_EQ(result.appliedRecoverChiSquared, recover);
		const auto stable = included;
		EXPECT_EQ(Calculations::applyHystereticReject(observations, included, chiSquared, cfg).status,
		          Calculations::RejectStatus::NoChange);
		EXPECT_EQ(included, stable);
	}
}

TEST(HystereticLeverageTest, EscalatesBothChiSquaredThresholdsAndHandlesUnsortedTimes)
{
	std::vector<ObsNode> observations(6);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000.0 + (5 - i);
	std::vector<std::uint8_t> included(6, 1);
	const std::array<double, 6> chiSquared{20, 20, 0, 0, 20, 20};
	const auto result = Calculations::applyHystereticReject(observations, included, chiSquared);
	EXPECT_EQ(included, (std::vector<std::uint8_t>(6, 1)));
	EXPECT_DOUBLE_EQ(result.appliedRejectChiSquared, 9.0 * 2.25);
	EXPECT_DOUBLE_EQ(result.appliedRecoverChiSquared, 8.289 * 2.25);
}

TEST(HystereticLeverageTest, RejectsInvalidThresholdsWithoutChangingMask)
{
	std::vector<ObsNode> observations(6);
	std::vector<std::uint8_t> included(6, 1);
	const std::array<double, 6> chiSquared{};
	for (double reject :
	     {0.0, 8.0, -1.0, 3600.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
		Calculations::HystereticRejectConfig cfg;
		cfg.rejectChiSquared = reject;
		EXPECT_THROW((void)Calculations::applyHystereticReject(observations, included, chiSquared, cfg),
		             std::invalid_argument);
		EXPECT_EQ(included, (std::vector<std::uint8_t>(6, 1)));
	}
}

TEST(HystereticLeverageTest, RadarIsNeverRejectedAutomaticallyButHonoursManualExclusion)
{
	std::vector<ObsNode> observations(12);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000.0 + i;
	for (const std::size_t i : {1u, 2u, 4u})
		observations[i].radar = orbital_core::obs::RadarObservation{};
	std::vector<std::uint8_t> included{1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1};
	const std::array<double, 12> chiSquared{0, 1.0e6, 1.0, 20, 1.0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<std::uint8_t> locked(12, 0);
	locked[4] = 1;
	const auto result = Calculations::applyHystereticReject(observations, included, chiSquared, {}, locked);
	EXPECT_EQ(included[1], 1);
	EXPECT_EQ(included[2], 1);
	EXPECT_EQ(included[3], 0);
	EXPECT_EQ(included[4], 0);
	EXPECT_EQ(result.nFlippedOut, 1);
	EXPECT_EQ(result.nFlippedIn, 1);
}

TEST(HystereticLeverageTest, RejectsWorstFirstAndCanBreakBatchCycles)
{
	std::vector<ObsNode> observations(12);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000.0 + i;
	std::vector<std::uint8_t> included(12, 1);
	std::array<double, 12> chiSquared{0, 100, 20, 11, 0, 0, 0, 0, 0, 0, 0, 0};
	const auto first = Calculations::applyHystereticReject(observations, included, chiSquared);
	EXPECT_EQ(first.nFlippedOut, 1);
	EXPECT_EQ(included[1], 0);
	EXPECT_EQ(included[2], 1);
	EXPECT_EQ(included[3], 1);
	Calculations::HystereticRejectConfig cfg;
	cfg.rejectLargestOnly = true;
	const auto single = Calculations::applyHystereticReject(observations, included, chiSquared, cfg);
	EXPECT_EQ(single.nFlippedOut, 1);
	EXPECT_EQ(included[2], 0);
	EXPECT_EQ(included[3], 1);
}

TEST(OppositionCoverageTest, LocalWindowHandlesUnsortedTimesDuplicatesAndMaskUpdates)
{
	std::vector<ObsNode> observations(4);
	observations[0].JD = 2459000 + 358;
	observations[1].JD = 2459000;
	observations[2].JD = 2459000 + 179;
	observations[3].JD = observations[0].JD;
	const std::array<std::uint8_t, 4> included{1, 1, 1, 0};
	orbital_core::fit::OppositionCoverage coverage(observations, included);
	EXPECT_FALSE(coverage.isLast(0));
	EXPECT_FALSE(coverage.isLast(1));
	coverage.exclude(2);
	EXPECT_TRUE(coverage.isLast(0));
	EXPECT_TRUE(coverage.isLast(1));
	EXPECT_FALSE(coverage.isLast(2));
	coverage.include(3);
	EXPECT_FALSE(coverage.isLast(0));
	EXPECT_FALSE(coverage.isLast(3));
	coverage.exclude(0);
	EXPECT_TRUE(coverage.isLast(3));
}

TEST(OppositionCoverageTest, BoundaryIsInclusiveAt180Days)
{
	std::vector<ObsNode> observations(2);
	observations[0].JD = 2459000;
	observations[1].JD = observations[0].JD + 180;
	const std::array<std::uint8_t, 2> included{1, 1};
	EXPECT_FALSE(orbital_core::fit::OppositionCoverage(observations, included).isLast(0));
	observations[1].JD = std::nextafter(observations[1].JD, std::numeric_limits<double>::infinity());
	EXPECT_TRUE(orbital_core::fit::OppositionCoverage(observations, included).isLast(0));
}

TEST(HystereticLeverageTest, ProtectedOppositionDoesNotRaiseTheBatchThreshold)
{
	std::vector<ObsNode> observations(12);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000 + i;
	observations[6].JD = 2459400;
	observations[6].note2 = 'E';
	std::vector<std::uint8_t> included(12, 1);
	std::array<double, 12> chiSquared{};
	chiSquared[6] = 1e10;
	chiSquared[3] = 20;
	const auto result = Calculations::applyHystereticReject(observations, included, chiSquared);
	EXPECT_EQ(result.nFlippedOut, 1);
	EXPECT_EQ(included[6], 1);
	EXPECT_EQ(included[3], 0);
	EXPECT_DOUBLE_EQ(result.appliedRejectChiSquared, 9.0);
}

TEST(HystereticLeverageTest, BatchCannotRemoveAnEntireMiddleOpposition)
{
	std::vector<ObsNode> observations(12);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000 + (i / 4) * 400 + i % 4;
	std::vector<std::uint8_t> included(12, 1);
	std::array<double, 12> chiSquared{};
	for (std::size_t i = 4; i < 8; ++i)
		chiSquared[i] = 100;
	const auto result = Calculations::applyHystereticReject(observations, included, chiSquared);
	EXPECT_EQ(result.nFlippedOut, 3);
	EXPECT_EQ(included, (std::vector<std::uint8_t>{1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1}));
	EXPECT_EQ(Calculations::applyHystereticReject(observations, included, chiSquared).status,
	          Calculations::RejectStatus::NoChange);
}

TEST(HystereticLeverageTest, ManualExclusionOfTheOnlyObservationRemainsExcluded)
{
	std::vector<ObsNode> observations(12);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000 + i;
	observations.back().JD += 400;
	std::vector<std::uint8_t> included(12, 1), locked(12, 0);
	included.back() = 0;
	locked.back() = 1;
	std::array<double, 12> chiSquared{};
	chiSquared[2] = 20;
	const auto result = Calculations::applyHystereticReject(observations, included, chiSquared, {}, locked);
	EXPECT_EQ(result.nFlippedOut, 1);
	EXPECT_EQ(result.nFlippedIn, 0);
	EXPECT_EQ(included.back(), 0);
	EXPECT_EQ(included[2], 0);
}
