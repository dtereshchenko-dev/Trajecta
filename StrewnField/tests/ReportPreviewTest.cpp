#include "../output/ReportPreview.h"
#include <gtest/gtest.h>

namespace {
	using namespace strewn::output;

	TEST(ReportPreviewTest, RetainsAllSmallTreesAndBoundsLargeOnesWithAncestry)
	{
		std::vector<PreviewTrack> nodes{{0, 0, 1000}, {1, 0, 100}};
		for (std::uint64_t i = 2; i < 1002; ++i)
			nodes.push_back({i, 1, 1002. - static_cast<double>(i)});
		EXPECT_EQ(previewTracks(nodes, nodes.size()).size(), nodes.size());
		const auto chosen = previewTracks(nodes, 32);
		EXPECT_EQ(chosen.size(), 32U);
		EXPECT_TRUE(chosen.contains(0));
		EXPECT_TRUE(chosen.contains(1));
		EXPECT_TRUE(chosen.contains(2));
		EXPECT_TRUE(chosen.contains(1001));
		for (const auto& node : nodes)
			if (chosen.contains(node.id))
				EXPECT_TRUE(chosen.contains(node.parent));
		std::reverse(nodes.begin(), nodes.end());
		EXPECT_EQ(chosen, previewTracks(nodes, 32));
		EXPECT_EQ(previewLimit(1000), 32U);
		EXPECT_EQ(previewLimit(1), 500U);
	}

	TEST(ReportPreviewTest, DeepBranchNeverCreatesAnOrphanOrExceedsLimit)
	{
		std::vector<PreviewTrack> nodes{{0, 0, 100}};
		for (std::uint64_t i = 1; i < 100; ++i)
			nodes.push_back({i, i - 1, 100. - static_cast<double>(i)});
		const auto chosen = previewTracks(nodes, 8);
		EXPECT_LE(chosen.size(), 8U);
		for (const auto& node : nodes)
			if (chosen.contains(node.id))
				EXPECT_TRUE(chosen.contains(node.parent));
		EXPECT_TRUE(previewTracks({}, 32).empty());
	}

	TEST(ReportPreviewTest, HistogramUsesAllWeightsAndExactEndpoints)
	{
		ArchivedField a;
		a.points = {{0, 0, 1, 3, 2}, {0, 0, 1, 2, 1}, {0, 0, 10, 1, 4}, {0, 0, 100, 0, 1}};
		const auto h = impactHistogram(a, [&](auto i) { return a.points[i].massKg; }, true);
		ASSERT_EQ(h["y"].size(), 28U);
		EXPECT_DOUBLE_EQ(h["y"][0].get<double>(), 8. / 12.);
		EXPECT_DOUBLE_EQ(h["y"][27].get<double>(), 4. / 12.);
		EXPECT_DOUBLE_EQ(h["ranges"][0][0], 1.);
		EXPECT_NEAR(h["ranges"][27][1].get<double>(), 10., 1e-12);
		const auto single = impactHistogram(a, [](auto) { return 59.; }, false);
		EXPECT_DOUBLE_EQ(single["singleValue"], 59.);
		EXPECT_EQ(single["y"], Json::array({1.}));
		EXPECT_TRUE(impactHistogram(ArchivedField{}, [](auto) { return 1.; }, false)["x"].empty());
	}
}  // namespace
