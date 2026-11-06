#include <gtest/gtest.h>

#include "hazard/LovReturnAnalysis.h"

#include <vector>

using namespace orbital_core::hazard;

namespace {
	std::vector<int> samplesOf(const LovReturnFilament& filament, const std::vector<LovTraceEntry>& entries)
	{
		std::vector<int> samples;
		for (const int index : filament.entryIndices)
			samples.push_back(entries[static_cast<std::size_t>(index)].sampleIndex);
		return samples;
	}
}  // namespace

TEST(LovReturnAnalysisTest, SplitsIndexRangesBeforeDuplicateReturns)
{
	const std::vector<LovTraceEntry> entries{
	    {0, 1, 100.0}, {1, 2, 101.0}, {2, 10, 101.5}, {3, 1, 102.0}, {4, 11, 102.5}};

	const std::vector<LovReturnFilament> returns = buildLovReturnFilaments(entries, LovReturnAnalysisConfig{});

	ASSERT_EQ(returns.size(), 3u);
	EXPECT_EQ(samplesOf(returns[0], entries), (std::vector<int>{1, 2}));
	EXPECT_EQ(samplesOf(returns[1], entries), (std::vector<int>{1}));
	EXPECT_EQ(samplesOf(returns[2], entries), (std::vector<int>{10, 11}));
}

TEST(LovReturnAnalysisTest, UsesShowerSpanAndGapCuts)
{
	const std::vector<LovTraceEntry> entries{{0, 1, 100.0}, {1, 2, 120.0}, {2, 3, 191.0}, {3, 4, 222.0}};

	const std::vector<LovReturnFilament> returns = buildLovReturnFilaments(entries, LovReturnAnalysisConfig{});

	ASSERT_EQ(returns.size(), 3u);
	EXPECT_EQ(samplesOf(returns[0], entries), (std::vector<int>{1, 2}));
	EXPECT_EQ(samplesOf(returns[1], entries), (std::vector<int>{3}));
	EXPECT_EQ(samplesOf(returns[2], entries), (std::vector<int>{4}));
}

TEST(LovReturnAnalysisTest, SplitsEqualEpochDuplicateIntoTwoReturnsDeterministically)
{
	const std::vector<LovTraceEntry> entries{
	    {0, 338, 70117.22599}, {1, 339, 70117.17576}, {2, 339, 70117.17576}, {3, 340, 70117.12580}};

	const std::vector<LovReturnFilament> returns = buildLovReturnFilaments(entries, LovReturnAnalysisConfig{});

	ASSERT_EQ(returns.size(), 2u);
	EXPECT_EQ(samplesOf(returns[0], entries), (std::vector<int>{339, 340}));
	EXPECT_EQ(returns[0].entryIndices.front(), 1);
	EXPECT_EQ(samplesOf(returns[1], entries), (std::vector<int>{338, 339}));
	EXPECT_EQ(returns[1].entryIndices.back(), 2);
}

TEST(LovReturnAnalysisTest, TwoRepeated2020VvEndpointsAddTwoReturns)
{
	const std::vector<LovTraceEntry> entries{{0, 335, 70117.22364}, {1, 336, 70117.17302}, {2, 337, 70117.12265},
	                                         {3, 337, 70117.12265}, {4, 338, 70117.07263}, {5, 338, 70117.07263},
	                                         {6, 339, 70117.02286}, {7, 340, 70116.97335}};
	const auto returns = buildLovReturnFilaments(entries, {});
	ASSERT_EQ(returns.size(), 3u);
	EXPECT_EQ(samplesOf(returns[0], entries), (std::vector<int>{338, 339, 340}));
	EXPECT_EQ(samplesOf(returns[1], entries), (std::vector<int>{337, 338}));
	EXPECT_EQ(samplesOf(returns[2], entries), (std::vector<int>{335, 336, 337}));

	const std::vector<LovTraceEntry> unique{entries[0], entries[1], entries[2], entries[4], entries[6], entries[7]};
	const auto withoutDuplicates = buildLovReturnFilaments(unique, {});
	ASSERT_EQ(withoutDuplicates.size(), 1u);
	EXPECT_EQ(samplesOf(withoutDuplicates[0], unique), (std::vector<int>{335, 336, 337, 338, 339, 340}));
}
