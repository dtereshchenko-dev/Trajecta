#include <gtest/gtest.h>

#include "hazard/lov/LovTraceDetail.h"
#include "hazard/HazardAudit.h"

#include <array>
#include <cmath>
#include <vector>

using namespace orbital_core::hazard;
using namespace orbital_core::hazard::lov_detail;

namespace {
	TraceRecord point(double sigma, double distance, double derivative, double stretch = 10.0)
	{
		TraceRecord out;
		out.sampleIndex = 153;
		out.event.jd = 2467898.0;
		out.trace.jd = out.event.jd;
		out.trace.sigma = sigma;
		out.trace.valid = true;
		out.trace.bodyRadiusKm = 1.0;
		out.trace.coordsKm = Eigen::Vector2d(distance, 0.0);
		out.trace.distanceKm = distance;
		out.trace.distanceBodyRadii = distance;
		out.trace.dDistance2DSigmaBodyRadii2 = derivative;
		out.trace.stretchLovBodyRadiiPerSigma = stretch;
		out.trace.stretchBodyRadii = stretch;
		out.trace.impactCrossSectionBodyRadii = 1.0;
		return out;
	}
}  // namespace

TEST(LovMinimumSearchTest, Vv2020FinalNoTpDoesNotLeakPreviousPoint)
{
	struct Step {
		double index;
		double mjd;
		double b;
		double derivative;
		double stretchLov;
	};

	const std::array<Step, 8> steps{
	    {{154.5, 0, 0, 0, 0},
	     {154.75, 67897.869587057095, 29.011060189786264, 218063.10899077804, 3760.1717690604323},
	     {154.625, 67897.984260347337, 24.303554287345225, 183440.71548228897, 3776.5584475678497},
	     {154.5625, 0, 0, 0, 0},
	     {154.59375, 67898.011793660917, 23.123633536278955, 174702.94686230819, 3780.4424008093938},
	     {154.578125, 0, 0, 0, 0},
	     {154.5859375, 67898.018600228024, 22.828473291104274, 172513.39676232144, 3781.3958113322919},
	     {154.58203125, 0, 0, 0, 0}}};
	auto start = point(-3.46729991, 38.361881310438669, 285697.41497098125, 3724.8984334037286);
	start.event.jd = 2400000.5 + 67897.62225;
	double referenceJd = start.event.jd;
	std::size_t calls = 0;
	const LovEncounterEvaluator evaluate = [&](int sample, double sigma, double fromJd,
	                                           double toJd) -> std::optional<TraceRecord> {
		EXPECT_EQ(sample, 153);
		EXPECT_DOUBLE_EQ(fromJd, referenceJd);
		EXPECT_DOUBLE_EQ(toJd, referenceJd);
		EXPECT_LT(calls, steps.size());
		if (calls >= steps.size())
			return std::nullopt;
		const auto& step = steps[calls++];
		EXPECT_NEAR(sigma, start.trace.sigma + 0.01 * (step.index - 155.0), 2e-14);
		if (step.mjd == 0.0)
			return std::nullopt;
		auto next = point(sigma, step.b, step.derivative, step.stretchLov);
		next.event.jd = 2400000.5 + step.mjd;
		referenceJd = next.event.jd;
		return next;
	};

	const auto result = refineEndpointMinimum(start, 0.005, evaluate);
	EXPECT_EQ(calls, 8U);
	EXPECT_TRUE(result.noTargetPlane);
	EXPECT_FALSE(result.converged);
	EXPECT_FALSE(result.usable());
	ASSERT_TRUE(result.record);
	EXPECT_NEAR(result.record->trace.sigma, -3.471440535, 2e-14);
	EXPECT_FALSE(result.record->riskCandidate);
}

TEST(LovMinimumSearchTest, AllNoTpAttemptsRejectTheStartingPoint)
{
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int, double sigma, double, double) -> std::optional<TraceRecord> {
		EXPECT_NEAR(sigma, std::ldexp(0.5, -calls), 1e-15);
		++calls;
		return std::nullopt;
	};
	const auto result = refineEndpointMinimum(point(0, 100, -1e8), 0.5, evaluate);
	EXPECT_EQ(calls, 8);
	EXPECT_FALSE(result.usable());
	EXPECT_TRUE(result.noTargetPlane);
}

TEST(LovMinimumSearchTest, ZeroLimitReevaluatesEndpointBeforeConverging)
{
	const auto start = point(-3.47, 10.0, -1000.0);
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int sample, double sigma, double fromJd,
	                                           double toJd) -> std::optional<TraceRecord> {
		++calls;
		EXPECT_EQ(sample, start.sampleIndex);
		EXPECT_DOUBLE_EQ(sigma, start.trace.sigma);
		EXPECT_DOUBLE_EQ(fromJd, start.event.jd);
		EXPECT_DOUBLE_EQ(toJd, start.event.jd);
		return start;
	};
	const auto result = refineEndpointMinimum(start, 0.0, evaluate);
	EXPECT_EQ(calls, 1);
	EXPECT_TRUE(result.usable(true));
	EXPECT_FALSE(result.reclassified);
}

TEST(LovMinimumSearchTest, SuccessfulRetryClearsNoTp)
{
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int, double sigma, double, double) -> std::optional<TraceRecord> {
		if (++calls < 8)
			return std::nullopt;
		return point(sigma, 0.5, -1e8);
	};
	const auto result = refineEndpointMinimum(point(0, 100, -1e8), 0.5, evaluate);
	EXPECT_EQ(calls, 8);
	EXPECT_FALSE(result.noTargetPlane);
	EXPECT_TRUE(result.converged);
	EXPECT_TRUE(result.usable(true));
}

TEST(LovMinimumSearchTest, NonconvergedTpIsUsableOnlyWhenCallerAllowsIt)
{
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int, double sigma, double, double) -> std::optional<TraceRecord> {
		++calls;
		return point(sigma, 100.0 - calls, -1e8);
	};
	const auto result = refineEndpointMinimum(point(0, 100, -1e8), 0.5, evaluate);
	EXPECT_EQ(calls, 8);
	EXPECT_FALSE(result.noTargetPlane);
	EXPECT_FALSE(result.converged);
	EXPECT_TRUE(result.usable());
	EXPECT_FALSE(result.usable(true));
}

TEST(LovMinimumSearchTest, AchillesChecksImpactBeforeSwitchingToFalsi)
{
	auto start = point(0, 1.1, -10);
	auto inside = point(0.05, 0.9, 10);
	inside.trace.alphaLovRad = Constants::PI;
	ASSERT_EQ(classifyLovPair(start, inside), LovPairType::InterruptedMinimum);
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int, double sigma, double, double) -> std::optional<TraceRecord> {
		++calls;
		EXPECT_NEAR(sigma, inside.trace.sigma, 1e-15);
		return inside;
	};
	const auto result = refineEndpointMinimum(start, 0.1, evaluate);
	EXPECT_EQ(calls, 1);
	EXPECT_TRUE(result.usable(true));
	EXPECT_FALSE(result.reclassified);
}

TEST(LovMinimumSearchTest, AchillesUsesRelativeDistanceCriterion)
{
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int, double sigma, double, double) -> std::optional<TraceRecord> {
		++calls;
		return point(sigma, 9.99, -3.996);
	};
	const auto result = refineEndpointMinimum(point(0, 10, -2000), 0.1, evaluate);
	EXPECT_EQ(calls, 1);
	EXPECT_TRUE(result.usable(true));
}

TEST(LovMinimumSearchTest, AchillesPropagatesNoTpFromNestedFalsi)
{
	auto start = point(0, 10, -1e5, 100);
	auto next = point(1, 8, 1000, 100);
	next.trace.alphaLovRad = Constants::PI;
	ASSERT_EQ(classifyLovPair(start, next), LovPairType::InterruptedMinimum);
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int, double, double, double) -> std::optional<TraceRecord> {
		return ++calls == 1 ? std::optional<TraceRecord>(next) : std::nullopt;
	};
	const auto result = refineEndpointMinimum(start, 1.0, evaluate);
	EXPECT_EQ(calls, 2);
	EXPECT_TRUE(result.noTargetPlane);
	EXPECT_FALSE(result.usable());
	EXPECT_TRUE(result.reclassified);
}

TEST(LovMinimumSearchTest, AchillesPreservesChangedTypeUntilLaterConvergence)
{
	auto start = point(0, 10, -1e5, 100);
	start.trace.coordsKm = Eigen::Vector2d(0, 10);
	auto next = point(1, 8, 1000, 100);
	next.trace.coordsKm = Eigen::Vector2d(0, 8);
	ASSERT_EQ(classifyLovPair(start, next), LovPairType::SimpleMinimum);
	int calls = 0;
	const LovEncounterEvaluator evaluate = [&](int, double sigma, double, double) -> std::optional<TraceRecord> {
		if (++calls == 1)
			return next;
		return point(sigma, 0.5, 1000, 100);
	};
	const auto result = refineEndpointMinimum(start, 1.0, evaluate);
	EXPECT_EQ(calls, 2);
	EXPECT_TRUE(result.usable(true));
	EXPECT_TRUE(result.reclassified);
}

TEST(LovMinimumSearchTest, FalsiUsesPreviousEndpointDerivativesAndSuperaccelerator)
{
	for (bool mirror : {false, true}) {
		const auto left = point(0, 10, mirror ? -6 : -2);
		const auto right = point(2, 12, mirror ? 2 : 6);
		const std::array<double, 3> expected =
		    mirror ? std::array<double, 3>{1.5, 1.125, 0.84375} : std::array<double, 3>{0.5, 0.875, 1.325};
		std::size_t calls = 0;
		const LovEncounterEvaluator evaluate = [&](int, double sigma, double, double) -> std::optional<TraceRecord> {
			EXPECT_LT(calls, expected.size());
			if (calls >= expected.size())
				return std::nullopt;
			EXPECT_NEAR(sigma, expected[calls], 1e-14);
			++calls;
			return point(sigma, calls == 3 ? 0.5 : 9.0 - calls, (mirror ? 1.0 : -1.0) / calls);
		};
		const auto result = refineBracketedMinimum(left, right, evaluate);
		EXPECT_EQ(calls, 3U);
		EXPECT_TRUE(result.usable(true));
	}
}

TEST(LovMinimumSearchTest, FalsiNoTpRetainsDiagnosticsButNotRiskCandidate)
{
	const LovEncounterEvaluator evaluate = [](int, double, double, double) -> std::optional<TraceRecord> {
		return std::nullopt;
	};
	const auto result = refineBracketedMinimum(point(0, 10, -2), point(2, 12, 6), evaluate);
	EXPECT_TRUE(result.noTargetPlane);
	EXPECT_FALSE(result.converged);
	EXPECT_FALSE(result.usable());
	ASSERT_TRUE(result.record);
	EXPECT_FALSE(result.record->riskCandidate);
}

TEST(LovMinimumSearchTest, FocusedAuditIsInheritedAndRestoresItsParent)
{
	const auto* previous = audit::current;
	{
		audit::Scope outer(153, -3.47, true);
		ASSERT_NE(audit::current, nullptr);
		EXPECT_EQ(audit::current->sample, 153);
		{
			audit::Scope inner(154, -3.48);
			ASSERT_NE(audit::current, nullptr);
			EXPECT_EQ(audit::current->sample, 154);
			EXPECT_DOUBLE_EQ(audit::current->sigma, -3.48);
		}
		EXPECT_EQ(audit::current->sample, 153);
	}
	EXPECT_EQ(audit::current, previous);
}
