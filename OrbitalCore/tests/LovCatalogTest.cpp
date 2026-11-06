#include <gtest/gtest.h>

#include "hazard/lov/LovTraceDetail.h"
#include "hazard/LineOfVariation.h"
#include "orbital_core/hazard/LineOfVariationImpact.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

using namespace orbital_core::hazard::lov_detail;
using orbital_core::hazard::LovSamplingConfig;
using orbital_core::hazard::sampleLovSigmas;

TEST(LovCatalogTest, FixedSamplingDefaultsIncludeBothFiveSigmaEndpoints)
{
	const orbital_core::hazard::LineOfVariationImpactConfig defaults;
	EXPECT_FALSE(defaults.probabilitySampling);
	EXPECT_TRUE(defaults.coordinateScaling);
	EXPECT_EQ(defaults.encounterStepPoints, 100);
	EXPECT_DOUBLE_EQ(defaults.closeApproach.encounterRadiusAu, 0.1);
	EXPECT_DOUBLE_EQ(defaults.closeApproach.thresholdAu, 0.2);
	LovSamplingConfig cfg;
	cfg.sigmaMax = defaults.sigmaMax;
	cfg.sigmaStepMax = defaults.sigmaStepMax;
	cfg.maxPositiveSamples = defaults.maxPositiveSamples;
	const auto sigmas = sampleLovSigmas(cfg);
	ASSERT_EQ(sigmas.size(), 1001u);
	EXPECT_DOUBLE_EQ(sigmas.front(), -5.0);
	EXPECT_DOUBLE_EQ(sigmas[500], 0.0);
	EXPECT_DOUBLE_EQ(sigmas.back(), 5.0);
	for (std::size_t i = 0; i < sigmas.size(); ++i) {
		EXPECT_DOUBLE_EQ(sigmas[i], (static_cast<double>(i) - 500.0) * 0.01);
		EXPECT_DOUBLE_EQ(sigmas[i], -sigmas[sigmas.size() - 1 - i]);
	}
}

TEST(LovCatalogTest, FixedSamplingHonorsRangeAndNodeCap)
{
	LovSamplingConfig cfg;
	cfg.sigmaMax = 0.025;
	cfg.sigmaStepMax = 0.01;
	cfg.maxPositiveSamples = 100;
	EXPECT_EQ(sampleLovSigmas(cfg), (std::vector<double>{-0.02, -0.01, 0.0, 0.01, 0.02}));
	cfg.maxPositiveSamples = 1;
	EXPECT_EQ(sampleLovSigmas(cfg), (std::vector<double>{-0.01, 0.0, 0.01}));
	cfg.sigmaStepMax = std::numeric_limits<double>::quiet_NaN();
	EXPECT_EQ(sampleLovSigmas(cfg), (std::vector<double>{0.0}));
}

TEST(LovCatalogTest, ProbabilitySamplingRemainsExplicitAndNonuniform)
{
	LovSamplingConfig cfg;
	cfg.sigmaMax = 5.0;
	cfg.sigmaStepMax = 0.01;
	cfg.probabilityStep = 1.0e-7;
	cfg.diskRadiusBodyRadii = 1000.0;
	cfg.maxPositiveSamples = 10000;
	cfg.probabilitySampling = true;
	const auto sigmas = sampleLovSigmas(cfg);
	ASSERT_GT(sigmas.size(), 1001u);
	const auto center = sigmas.size() / 2;
	EXPECT_DOUBLE_EQ(sigmas[center], 0.0);
	EXPECT_DOUBLE_EQ(sigmas[center + 1],
	                 (0.5 * cfg.diskRadiusBodyRadii * cfg.probabilityStep) / 0.398942280401432677939946059934381868);
	EXPECT_LT(sigmas[center + 1], cfg.sigmaStepMax);
	EXPECT_LE(sigmas.back(), cfg.sigmaMax);
	for (std::size_t i = 1; i < sigmas.size(); ++i) {
		EXPECT_GT(sigmas[i], sigmas[i - 1]);
		EXPECT_DOUBLE_EQ(sigmas[i], -sigmas[sigmas.size() - 1 - i]);
	}
}

namespace {
	LovCatalogPoint point(int index, double sigma, bool valid = true)
	{
		LovCatalogPoint result;
		result.sampleIndex = index;
		result.sigma = sigma;
		result.valid = valid;
		return result;
	}
}  // namespace

TEST(LovCatalogTest, HalfIndexRoundsUpOnEitherSigmaBranch)
{
	for (double center : {-2.0, 0.0, 2.0}) {
		std::array catalog{point(10, center - 0.25), point(11, center + 0.25)};
		for (int ordering = 0; ordering < 2; ++ordering) {
			const auto* selected = catalogAnchorAtSigma(catalog, center);
			ASSERT_NE(selected, nullptr);
			EXPECT_EQ(selected->sampleIndex, 11);
			std::reverse(catalog.begin(), catalog.end());
		}
	}
}

TEST(LovCatalogTest, Native2050HalfIndexUsesUpperAnchor)
{
	std::array catalog{point(834, -1.07072262), point(835, -1.06968082)};
	constexpr double sigma = -1.07020172;
	ASSERT_EQ(std::abs(catalog[0].sigma - sigma), std::abs(catalog[1].sigma - sigma));
	for (int ordering = 0; ordering < 2; ++ordering) {
		const auto* selected = catalogAnchorAtSigma(catalog, sigma);
		ASSERT_NE(selected, nullptr);
		EXPECT_EQ(selected->sampleIndex, 835);
		std::reverse(catalog.begin(), catalog.end());
	}
}

TEST(LovCatalogTest, NonuniformCatalogIsIndependentOfStorageOrder)
{
	const std::array points{point(20, -4.0), point(21, -3.0), point(22, -2.75), point(23, 0.0)};
	std::array order{0, 1, 2, 3};
	do {
		std::vector<LovCatalogPoint> catalog;
		for (int index : order)
			catalog.push_back(points[index]);
		for (const auto& [sigma, index] : std::array{std::pair{-3.5, 21}, std::pair{-2.875, 22}, std::pair{-1.375, 23},
		                                             std::pair{-3.75, 20}, std::pair{-2.9, 21}, std::pair{-2.8, 22}}) {
			const auto* selected = catalogAnchorAtSigma(catalog, sigma);
			ASSERT_NE(selected, nullptr);
			EXPECT_EQ(selected->sampleIndex, index);
		}
	} while (std::next_permutation(order.begin(), order.end()));
}

TEST(LovCatalogTest, DistinctDistancesDoNotBecomeTies)
{
	const std::array catalog{point(30, -2.0), point(31, -1.0)};
	const auto* lower = catalogAnchorAtSigma(catalog, std::nextafter(-1.5, -2.0));
	const auto* upper = catalogAnchorAtSigma(catalog, std::nextafter(-1.5, -1.0));
	ASSERT_NE(lower, nullptr);
	ASSERT_NE(upper, nullptr);
	EXPECT_EQ(lower->sampleIndex, 30);
	EXPECT_EQ(upper->sampleIndex, 31);
}

TEST(LovCatalogTest, ExactNodesAndOutsideRangeKeepTheirAnchors)
{
	const std::array catalog{point(40, -1.0), point(41, 0.0), point(42, 1.0)};
	for (std::size_t i = 0; i < catalog.size(); ++i)
		EXPECT_EQ(catalogAnchorAtSigma(catalog, catalog[i].sigma), &catalog[i]);
	EXPECT_EQ(catalogAnchorAtSigma(catalog, -5.0), &catalog.front());
	EXPECT_EQ(catalogAnchorAtSigma(catalog, 5.0), &catalog.back());
}

TEST(LovCatalogTest, InvalidEntriesAndRequestsDoNotSelectAnchors)
{
	const double nan = std::numeric_limits<double>::quiet_NaN();
	const double infinity = std::numeric_limits<double>::infinity();
	const std::array catalog{point(50, 0.0, false), point(51, nan), point(52, infinity), point(53, 1.0)};
	EXPECT_EQ(catalogAnchorAtSigma({}, 0.0), nullptr);
	EXPECT_EQ(catalogAnchorAtSigma(catalog, nan), nullptr);
	EXPECT_EQ(catalogAnchorAtSigma(catalog, infinity), nullptr);
	EXPECT_EQ(catalogAnchorAtSigma(catalog, -infinity), nullptr);
	EXPECT_EQ(catalogAnchorAtSigma(std::span(catalog).first(3), 0.0), nullptr);
	EXPECT_EQ(catalogAnchorAtSigma(catalog, 0.0), &catalog.back());
}
