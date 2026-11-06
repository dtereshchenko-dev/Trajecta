#include "strewn_core/FieldOutput.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace {
	using namespace strewn;

	MapOptions options()
	{
		MapOptions o;
		o.method = KdeMethod::Fixed;
		o.centerDeg = std::array<double, 2>{49.8, .8};
		o.halfSizeKm = {1, 1};
		o.fixedSigmaM = 50;
		return o;
	}

	TEST(FieldTest, FixedKernelMatchesDiscreteGaussianAndIncludesMultiplicity)
	{
		auto o = options();
		const std::vector<FieldPoint> points = {{49.8, .8, 2, 5, 3}};
		const auto g = estimateField(points, o);
		EXPECT_EQ(g.rows, 80U);
		EXPECT_EQ(g.columns, 80U);
		EXPECT_DOUBLE_EQ(g.count.totalWeight, 15);
		EXPECT_DOUBLE_EQ(g.mass.totalWeight, 30);
		EXPECT_EQ(g.count.probability, g.mass.probability);
		EXPECT_NEAR(std::accumulate(g.count.probability.begin(), g.count.probability.end(), 0.), 1., 1e-13);
		const auto center = (g.rows / 2) * g.columns + g.columns / 2;
		const auto& p = g.count.probability;
		EXPECT_NEAR(p[center + 1] / p[center], std::exp(-.125), 1e-12);
		EXPECT_NEAR(p[center + g.columns] / p[center], std::exp(-.125), 1e-12);
		for (std::size_t i = 0; i < 3; ++i) {
			double inside = 0;
			for (double value : p)
				if (value >= g.count.hpdThresholds[i])
					inside += value;
			const std::array<double, 3> levels{.5, .8, .99};
			EXPECT_GE(inside, levels[i] - 1e-12);
		}
	}

	TEST(FieldTest, FixedGroupedAndExpandedWeightsAgree)
	{
		const std::vector<FieldPoint> grouped = {{49.8, .8, .01, 4, 2}, {49.801, .802, .2, 1, 1}};
		const std::vector<FieldPoint> expanded = {{49.8, .8, .01, 1, 2},
		                                          {49.8, .8, .01, 1, 2},
		                                          {49.8, .8, .01, 1, 2},
		                                          {49.8, .8, .01, 1, 2},
		                                          {49.801, .802, .2, 1, 1}};
		const auto a = estimateField(grouped, options()), b = estimateField(expanded, options());
		EXPECT_EQ(a.count.probability, b.count.probability);
		EXPECT_EQ(a.mass.probability, b.mass.probability);
		EXPECT_NE(a.count.probability, a.mass.probability);
	}

	TEST(FieldTest, AdaptiveHasSeparateCountAndMassBandwidths)
	{
		auto o = options();
		o.method = KdeMethod::Adaptive;
		o.halfSizeKm = {4, 4};
		std::vector<FieldPoint> p;
		for (int i = 0; i < 40; ++i)
			p.push_back({49.8 + .001 * std::sin(i * 1.3), .8 + .003 * std::cos(i * .8), .01 + i * .005, 1. + i % 5, 1});
		const auto g = estimateField(p, o);
		EXPECT_FALSE(g.count.fixedKernelFallback);
		EXPECT_NE(g.count.covarianceM2, g.mass.covarianceM2);
		EXPECT_NEAR(std::accumulate(g.mass.probability.begin(), g.mass.probability.end(), 0.), 1., 1e-12);
		EXPECT_GT(g.count.hpdThresholds[0], g.count.hpdThresholds[1]);
		EXPECT_GT(g.count.hpdThresholds[1], g.count.hpdThresholds[2]);
	}

	TEST(FieldTest, SingularAdaptiveFallsBackExplicitly)
	{
		auto o = options();
		o.method = KdeMethod::Adaptive;
		const auto g = estimateField(std::vector<FieldPoint>{{49.8, .8, .1, 1, 1}}, o);
		EXPECT_TRUE(g.count.fixedKernelFallback);
		EXPECT_TRUE(g.mass.fixedKernelFallback);
		EXPECT_DOUBLE_EQ(g.count.covarianceM2[0], 2500);
	}

	TEST(FieldTest, OversizedGridAutomaticallyCoarsensWithoutChangingKernel)
	{
		auto o = options();
		o.halfSizeKm = {26, 26};
		const std::vector<FieldPoint> points = {{49.8, .8, .1, 1, 1}};
		const auto automatic = estimateField(points, o);
		EXPECT_DOUBLE_EQ(automatic.cellM, 50);
		EXPECT_LE(automatic.rows * automatic.columns, 4000000U);
		EXPECT_DOUBLE_EQ(o.cellM, 25);
		o.autoCellSize = false;
		EXPECT_THROW(estimateField(points, o), std::invalid_argument);
		o.cellM = 50;
		const auto explicitGrid = estimateField(points, o);
		EXPECT_EQ(automatic.count.probability, explicitGrid.count.probability);
		EXPECT_EQ(automatic.mass.probability, explicitGrid.mass.probability);
		EXPECT_EQ(automatic.count.covarianceM2, explicitGrid.count.covarianceM2);
	}

	TEST(FieldTest, RejectsClippingBadWeightsAndCancellation)
	{
		const auto o = options();
		EXPECT_THROW(estimateField({}, o), std::invalid_argument);
		EXPECT_THROW(estimateField(std::vector<FieldPoint>{{50, .8, .1, 1, 1}}, o), std::invalid_argument);
		EXPECT_THROW(estimateField(std::vector<FieldPoint>{{49.8, .8, .1, -1, 1}}, o), std::invalid_argument);
		auto small = o;
		small.halfSizeKm = {.1, .1};
		EXPECT_THROW(estimateField(std::vector<FieldPoint>{{49.8, .8, .1, 1, 1}}, small), std::invalid_argument);
		std::stop_source stop;
		stop.request_stop();
		EXPECT_THROW(estimateField(std::vector<FieldPoint>{{49.8, .8, .1, 1, 1}}, o, stop.get_token()),
		             std::runtime_error);
		auto huge = o;
		huge.fixedSigmaM = 1e300;
		EXPECT_THROW(estimateField(std::vector<FieldPoint>{{49.8, .8, .1, 1, 1}}, huge), std::invalid_argument);
	}

}  // namespace
