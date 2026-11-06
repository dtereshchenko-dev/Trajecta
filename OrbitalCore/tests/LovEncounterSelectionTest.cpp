#include <gtest/gtest.h>
#include <algorithm>

#include "hazard/lov/LovTraceDetail.h"

using namespace orbital_core::hazard;
using namespace orbital_core::hazard::lov_detail;

namespace {
	CloseApproachEvent encounter(double jd, double distanceAu, int encounterIndex)
	{
		CloseApproachEvent event;
		event.jd = jd;
		event.distAu = distanceAu;
		event.encounterIndex = encounterIndex;
		return event;
	}

	std::vector<CloseApproachPassage> passages(const std::vector<CloseApproachEvent>& events)
	{
		std::vector<CloseApproachPassage> out;
		for (const auto& event : events) {
			auto it = std::find_if(out.begin(), out.end(),
			                       [&](const auto& passage) { return passage.encounterIndex == event.encounterIndex; });
			if (it == out.end()) {
				out.push_back({event.encounterIndex, event.jd, event.jd});
			} else {
				it->firstInsideJd = std::min(it->firstInsideJd, event.jd);
				it->lastInsideJd = std::max(it->lastInsideJd, event.jd);
			}
		}
		return out;
	}

	std::optional<TraceEvaluation> trace(double distanceKm)
	{
		TraceEvaluation out;
		out.trace.valid = true;
		out.trace.distanceKm = distanceKm;
		return out;
	}
}  // namespace

TEST(LovEncounterSelectionTest, SelectsLowestMinimumInLastForwardEpisode)
{
	const std::vector<CloseApproachEvent> events{encounter(90.0, 0.001, 0), encounter(100.0, 0.100, 1),
	                                             encounter(102.0, 0.050, 1)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(300.0), trace(200.0), trace(100.0)};

	const std::optional<std::size_t> selected =
	    selectReturnTargetPlaneMinimum(events, traces, passages(events), 85.0, 115.0, 50.0, 100.0);

	ASSERT_TRUE(selected.has_value());
	EXPECT_EQ(*selected, 2U);
}

TEST(LovEncounterSelectionTest, IgnoresEncountersOutsideWindow)
{
	const std::vector<CloseApproachEvent> events{encounter(80.0, 0.001, 0), encounter(100.0, 0.100, 1),
	                                             encounter(120.0, 0.001, 2)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(1.0), trace(100.0), trace(1.0)};

	const std::optional<std::size_t> selected =
	    selectReturnTargetPlaneMinimum(events, traces, passages(events), 90.0, 110.0, 50.0, 100.0);

	ASSERT_TRUE(selected.has_value());
	EXPECT_EQ(*selected, 1U);
}

TEST(LovEncounterSelectionTest, SelectsLowestMinimumInLastBackwardEpisode)
{
	const std::vector<CloseApproachEvent> events{encounter(80.0, 0.050, 2), encounter(90.0, 0.001, 1),
	                                             encounter(78.0, 0.100, 2)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(100.0), trace(1.0), trace(200.0)};

	const std::optional<std::size_t> selected =
	    selectReturnTargetPlaneMinimum(events, traces, passages(events), 70.0, 95.0, 100.0, 80.0);

	ASSERT_TRUE(selected.has_value());
	EXPECT_EQ(*selected, 0U);
}

TEST(LovEncounterSelectionTest, DoesNotFallBackFromInvalidTerminalEpisode)
{
	const std::vector<CloseApproachEvent> events{encounter(90.0, 0.001, 1), encounter(100.0, 0.100, 2)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(1.0), std::nullopt};

	EXPECT_FALSE(selectReturnTargetPlaneMinimum(events, traces, passages(events), 85.0, 105.0, 50.0, 95.0).has_value());
}

TEST(LovEncounterSelectionTest, Vv2020UnfinishedEntryResetsMinimum)
{
	const std::vector<CloseApproachEvent> events{encounter(67898.091326, 0.000523449, 1)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(100.0)};
	const std::vector<CloseApproachPassage> history{{1, 67800.0, 67997.306412}, {2, 68140.085019, 68199.852787}};
	EXPECT_FALSE(selectReturnTargetPlaneMinimum(events, traces, history, 67600.0, 68199.852787, 59169.7715, 67898.0));

	EXPECT_EQ(selectReturnTargetPlaneMinimum(events, traces, std::span(history).first(1), 67600.0, 68199.852787,
	                                         59169.7715, 67898.0),
	          0U);
}

TEST(LovEncounterSelectionTest, UnfinishedBackwardEntryDoesNotReuseLaterMinimum)
{
	const std::vector<CloseApproachEvent> events{encounter(90.0, .001, 1)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(100.0)};
	const std::vector<CloseApproachPassage> history{{2, 70.0, 80.0}, {1, 85.0, 95.0}};
	EXPECT_FALSE(selectReturnTargetPlaneMinimum(events, traces, history, 70.0, 98.0, 100.0, 85.0));
}

TEST(LovEncounterSelectionTest, LastMinimumOutsideWindowRejectsWholeEpisode)
{
	const std::vector<CloseApproachEvent> events{encounter(100.0, .001, 1), encounter(112.0, .1, 1)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(1.0), trace(10.0)};
	const std::vector<CloseApproachPassage> history{{1, 90.0, 110.0}};
	EXPECT_FALSE(selectReturnTargetPlaneMinimum(events, traces, history, 90.0, 110.0, 50.0, 100.0));
}

TEST(LovEncounterSelectionTest, WindowCheckPrecedesWholeEpisodeMinimumSelection)
{
	const std::vector<CloseApproachEvent> events{encounter(89.0, .001, 1), encounter(100.0, .1, 1)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(1.0), trace(10.0)};
	const std::vector<CloseApproachPassage> history{{1, 85.0, 105.0}};
	EXPECT_EQ(selectReturnTargetPlaneMinimum(events, traces, history, 90.0, 110.0, 50.0, 100.0), 0U);
}

TEST(LovEncounterSelectionTest, MissingEpisodeHistoryDoesNotGuessFromMinima)
{
	const std::vector<CloseApproachEvent> events{encounter(100.0, .001, 1)};
	const std::vector<std::optional<TraceEvaluation>> traces{trace(1.0)};
	EXPECT_FALSE(selectReturnTargetPlaneMinimum(events, traces, {}, 90.0, 110.0, 50.0, 100.0));
}

TEST(LovEncounterSelectionTest, BuildsForwardCheckpointChainBeforeTarget)
{
	const std::vector<CloseApproachEvent> events{encounter(130.0, 0.1, 3), encounter(110.0, 0.1, 1),
	                                             encounter(120.0, 0.1, 2), encounter(110.0, 0.2, 1),
	                                             encounter(90.0, 0.1, 0)};

	const std::vector<double> checkpoints = encounterCheckpointJds(events, 100.0, 130.0);

	EXPECT_EQ(checkpoints, (std::vector<double>{110.0, 120.0}));
}

TEST(LovEncounterSelectionTest, BuildsBackwardCheckpointChainBeforeTarget)
{
	const std::vector<CloseApproachEvent> events{encounter(70.0, 0.1, 3), encounter(90.0, 0.1, 1),
	                                             encounter(80.0, 0.1, 2), encounter(110.0, 0.1, 0)};

	const std::vector<double> checkpoints = encounterCheckpointJds(events, 100.0, 70.0);

	EXPECT_EQ(checkpoints, (std::vector<double>{90.0, 80.0}));
}
