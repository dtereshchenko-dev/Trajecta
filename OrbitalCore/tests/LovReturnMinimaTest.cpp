#include <gtest/gtest.h>

#include "hazard/lov/LovTraceDetail.h"

#include <array>
#include <cstddef>
#include <vector>

using namespace orbital_core::hazard::lov_detail;

namespace {
	constexpr double kSigmaMax = 5.0;

	TraceRecord point(int sample, double sigma, double distance, double derivative, double angle = 0.0)
	{
		TraceRecord out;
		out.sampleIndex = sample;
		out.event.localMoidAu = distance / kAuKm;
		out.event.distAu = out.event.localMoidAu;
		out.event.vInfKmS = 0.0;
		out.trace.valid = true;
		out.trace.sigma = sigma;
		out.trace.bodyRadiusKm = 1.0;
		out.trace.coordsKm = Eigen::Vector2d(distance, 0.0);
		out.trace.distanceKm = distance;
		out.trace.distanceBodyRadii = distance;
		out.trace.widthBodyRadii = 0.0;
		out.trace.stretchBodyRadii = 1.0;
		out.trace.stretchLovBodyRadiiPerSigma = 1.0;
		out.trace.dDistance2DSigmaBodyRadii2 = derivative;
		out.trace.alphaLovRad = angle;
		return out;
	}

	LovMinimumResult result(const TraceRecord& record, bool converged = true, bool noTargetPlane = false,
	                        bool reclassified = false)
	{
		LovMinimumResult out;
		out.record = record;
		out.converged = converged;
		out.noTargetPlane = noTargetPlane;
		out.reclassified = reclassified;
		return out;
	}

	std::vector<TraceRecord> interruptedPair()
	{
		return {point(1, -1.0, 100.0, -1.0), point(2, 1.0, 10.0, 1.0, Constants::PI)};
	}

	struct EndpointCall {
		int sample;
		double limit;
	};

	struct LimitCall {
		int sample;
		bool increasingHead;
	};

	struct ReturnReplay {
		LovMinimumResult bracketResult;
		std::vector<LovMinimumResult> endpointResults;
		std::vector<EndpointCall> endpointCalls;
		std::vector<LimitCall> limitCalls;
		int bracketCalls = 0;

		std::vector<TraceRecord> run(const std::vector<TraceRecord>& records)
		{
			return addRefinedMinima(
			    records, kSigmaMax,
			    [&](const TraceRecord& record, bool increasingHead) {
				    limitCalls.push_back({record.sampleIndex, increasingHead});
				    return 0.125;
			    },
			    [&](const TraceRecord& record, double limit) -> LovMinimumResult {
				    const std::size_t index = endpointCalls.size();
				    endpointCalls.push_back({record.sampleIndex, limit});
				    EXPECT_LT(index, endpointResults.size());
				    return index < endpointResults.size() ? endpointResults[index] : LovMinimumResult{};
			    },
			    [&](const TraceRecord&, const TraceRecord&) {
				    ++bracketCalls;
				    return bracketResult;
			    });
		}
	};
}  // namespace

TEST(LovReturnMinimaTest, AdmittedSimplePairRefinesBothEndsWithoutAnotherDistanceFilter)
{
	for (bool highLeft : {false, true}) {
		SCOPED_TRACE(highLeft);
		std::vector<TraceRecord> records{point(1, -1.0, highLeft ? 100.0 : 10.0, -1.0),
		                                 point(2, 1.0, highLeft ? 10.0 : 100.0, 1.0)};
		for (auto& record : records)
			record.trace.coordsKm = Eigen::Vector2d(0.0, record.trace.distanceKm);
		ASSERT_EQ(classifyLovPair(records[0], records[1]), LovPairType::SimpleMinimum);
		ASSERT_GT(minimumPossibleBodyRadii(records[highLeft ? 0 : 1], kSigmaMax), 60.0);
		ASSERT_LT(minimumPossibleBodyRadii(records[highLeft ? 1 : 0], kSigmaMax), 60.0);

		ReturnReplay replay;
		replay.endpointResults = {result(point(10, -0.5, 5.0, 0.0), highLeft),
		                          result(point(11, 0.5, 4.0, 0.0), !highLeft)};
		const auto output = replay.run(records);
		ASSERT_EQ(replay.endpointCalls.size(), 2U);
		EXPECT_EQ(replay.endpointCalls[0].sample, 1);
		EXPECT_DOUBLE_EQ(replay.endpointCalls[0].limit, 1.0);
		EXPECT_EQ(replay.endpointCalls[1].sample, 2);
		EXPECT_DOUBLE_EQ(replay.endpointCalls[1].limit, 0.0);
		EXPECT_TRUE(replay.limitCalls.empty());
		EXPECT_EQ(replay.bracketCalls, 1);
		ASSERT_EQ(output.size(), 1U);
		EXPECT_EQ(output[0].sampleIndex, highLeft ? 10 : 11);
	}
}

TEST(LovReturnMinimaTest, InterruptedExtraUsesLastRightFallbackEvenWithNoTpDiagnostics)
{
	const auto records = interruptedPair();
	ASSERT_EQ(classifyLovPair(records[0], records[1]), LovPairType::InterruptedMinimum);
	for (bool anchorNearRight : {false, true}) {
		for (bool rightNoTp : {false, true}) {
			SCOPED_TRACE(::testing::Message() << "nearRight=" << anchorNearRight << ", noTp=" << rightNoTp);
			ReturnReplay replay;
			const double anchorAngle = anchorNearRight ? Constants::PI : 0.0;
			const double previousAngle = anchorNearRight ? 0.0 : Constants::PI;
			replay.bracketResult = result(point(9, -0.7, 20.0, 0.0, previousAngle), false, true);
			replay.endpointResults = {result(point(10, -0.5, 15.0, 0.0, previousAngle)),
			                          result(point(11, 0.6, 8.0, 0.0, anchorAngle), true, rightNoTp),
			                          result(point(12, 0.0, 5.0, 0.0))};

			const auto output = replay.run(records);
			ASSERT_EQ(replay.endpointCalls.size(), 3U);
			EXPECT_EQ(replay.endpointCalls[0].sample, 1);
			EXPECT_EQ(replay.endpointCalls[1].sample, 2);
			EXPECT_DOUBLE_EQ(replay.endpointCalls[1].limit, 0.0);
			EXPECT_EQ(replay.endpointCalls[2].sample, anchorNearRight ? 1 : 2);
			EXPECT_DOUBLE_EQ(replay.endpointCalls[2].limit, anchorNearRight ? 0.8 : 0.2);
			ASSERT_EQ(output.size(), rightNoTp ? 2U : 3U);
			EXPECT_EQ(output.front().sampleIndex, 10);
			EXPECT_EQ(output.back().sampleIndex, 12);
			if (!rightNoTp)
				EXPECT_EQ(output[1].sampleIndex, 11);
		}
	}
}

TEST(LovReturnMinimaTest, EitherReclassifiedFallbackSuppressesInterruptedExtra)
{
	for (bool leftChanged : {false, true}) {
		for (bool rightChanged : {false, true}) {
			SCOPED_TRACE(::testing::Message() << leftChanged << ", " << rightChanged);
			ReturnReplay replay;
			replay.endpointResults = {result(point(10, -0.5, 15.0, 0.0), true, false, leftChanged),
			                          result(point(11, 0.6, 8.0, 0.0, Constants::PI), true, true, rightChanged)};
			if (!leftChanged && !rightChanged)
				replay.endpointResults.push_back(result(point(12, 0.0, 5.0, 0.0)));
			const auto output = replay.run(interruptedPair());
			const bool changed = leftChanged || rightChanged;
			EXPECT_EQ(replay.endpointCalls.size(), changed ? 2U : 3U);
			EXPECT_EQ(output.size(), changed ? 1U : 2U);
		}
	}
}

TEST(LovReturnMinimaTest, RightFallbackUsesZeroStepAndBothFallbacksRequireConvergedTp)
{
	struct Case {
		bool leftConverged;
		bool rightConverged;
		bool leftNoTp;
		bool rightNoTp;
		std::size_t accepted;
	};

	const std::array<Case, 7> cases{{{false, false, false, false, 0},
	                                 {true, false, false, false, 1},
	                                 {false, true, false, false, 1},
	                                 {true, true, false, false, 2},
	                                 {true, true, true, false, 1},
	                                 {true, true, false, true, 1},
	                                 {true, true, true, true, 0}}};
	auto records = interruptedPair();
	records[1].trace.alphaLovRad = 0.0;
	for (auto& record : records)
		record.trace.coordsKm = Eigen::Vector2d(0.0, record.trace.distanceKm);
	ASSERT_EQ(classifyLovPair(records[0], records[1]), LovPairType::SimpleMinimum);
	for (const auto& test : cases) {
		SCOPED_TRACE(::testing::Message() << test.leftConverged << ", " << test.rightConverged << ", " << test.leftNoTp
		                                  << ", " << test.rightNoTp);
		ReturnReplay replay;
		replay.endpointResults = {result(point(10, -0.5, 5.0, 0.0), test.leftConverged, test.leftNoTp),
		                          result(point(11, 0.5, 4.0, 0.0), test.rightConverged, test.rightNoTp)};
		const auto output = replay.run(records);
		ASSERT_EQ(replay.endpointCalls.size(), 2U);
		EXPECT_EQ(replay.endpointCalls[0].sample, 1);
		EXPECT_DOUBLE_EQ(replay.endpointCalls[0].limit, 1.0);
		EXPECT_EQ(replay.endpointCalls[1].sample, 2);
		EXPECT_DOUBLE_EQ(replay.endpointCalls[1].limit, 0.0);
		EXPECT_EQ(output.size(), test.accepted);
	}
}
