#include <gtest/gtest.h>

#include "hazard/lov/LovTraceDetail.h"

#include "orbital_core/Constants.h"

using namespace orbital_core::hazard;
using namespace orbital_core::hazard::lov_detail;

namespace {
	TraceRecord traceRecord(const Eigen::Vector2d& coords, double alphaLov, double derivative)
	{
		TraceRecord record;
		record.event.vInfKmS = 10.0;
		record.trace.valid = true;
		record.trace.coordsKm = coords;
		record.trace.bodyRadiusKm = 6378.0;
		record.trace.distanceBodyRadii = coords.norm() / 6378.0;
		record.trace.stretchLovBodyRadiiPerSigma = 1.0;
		record.trace.alphaLovRad = alphaLov;
		record.trace.dDistance2DSigmaBodyRadii2 = derivative;
		return record;
	}
}  // namespace

TEST(LovTraceGeometryTest, TangentTestCanClassifyDerivativeBracketAsEntangled)
{
	const TraceRecord left = traceRecord({6378.0, 0.0}, 0.0, -1.0);
	const TraceRecord right = traceRecord({12756.0, 0.0}, 0.0, 1.0);

	EXPECT_EQ(classifyLovPair(left, right), LovPairType::EntangledAscending);
}

TEST(LovTraceGeometryTest, AlignedDerivativeBracketIsSimpleMinimum)
{
	const TraceRecord left = traceRecord({6378.0, 0.0}, 0.0, -1.0);
	const TraceRecord right = traceRecord({6378.0, 6378.0}, 0.0, 1.0);

	EXPECT_EQ(classifyLovPair(left, right), LovPairType::SimpleMinimum);
}

TEST(LovTraceGeometryTest, ReversedTangentIsInterruptedMinimum)
{
	const TraceRecord left = traceRecord({6378.0, 0.0}, 0.0, -1.0);
	const TraceRecord right = traceRecord({6378.0, 6378.0}, Constants::PI, 1.0);

	EXPECT_EQ(classifyLovPair(left, right), LovPairType::InterruptedMinimum);
}

TEST(LovTraceGeometryTest, CloseApproachWindowUsesConfiguredControlDistance)
{
	const double narrow = closeApproachWindowDays(0.1, 10.0);
	const double wide = closeApproachWindowDays(0.2, 10.0);
	EXPECT_NEAR(wide, 2.0 * narrow, 1.0e-12);
}

TEST(LovTraceGeometryTest, EndpointLimitUsesLocalAdaptiveCatalogStep)
{
	std::vector<LovCatalogPoint> catalog(3);
	catalog[0].sampleIndex = 40;
	catalog[0].sigma = 1.5000;
	catalog[1].sampleIndex = 41;
	catalog[1].sigma = 1.5020;
	catalog[2].sampleIndex = 42;
	catalog[2].sigma = 1.5045;

	EXPECT_NEAR(endpointSigmaStepLimit(catalog, 41, true, 0.05), 0.00125, 1.0e-15);
	EXPECT_NEAR(endpointSigmaStepLimit(catalog, 41, false, 0.05), 0.00100, 1.0e-15);
}

TEST(LovTraceGeometryTest, CatalogPredictorUsesOutgoingAdaptiveStep)
{
	const std::vector<double> sigmas{-0.07, -0.03, -0.01, 0.0, 0.01, 0.03, 0.07};
	constexpr std::size_t center = 3;

	EXPECT_DOUBLE_EQ(catalogPredictorStep(sigmas, center + 1, center), 0.02);
	EXPECT_DOUBLE_EQ(catalogPredictorStep(sigmas, center + 2, center), 0.04);
	EXPECT_DOUBLE_EQ(catalogPredictorStep(sigmas, center - 1, center), -0.02);
	EXPECT_DOUBLE_EQ(catalogPredictorStep(sigmas, center - 2, center), -0.04);
}

TEST(LovTraceGeometryTest, CatalogSigmaQUsesSignedObjectiveDistance)
{
	EXPECT_DOUBLE_EQ(catalogSigmaQ(34.0, 9.0), 5.0);
	EXPECT_DOUBLE_EQ(catalogSigmaQ(4.0, 9.0), -std::sqrt(5.0));
	EXPECT_DOUBLE_EQ(catalogSigmaQ(9.0, 9.0), 0.0);
}

TEST(LovTraceGeometryTest, DefaultWorkerCountLeavesOneThreadAvailable)
{
	EXPECT_EQ(defaultLovWorkerCount(1), 1);
	EXPECT_EQ(defaultLovWorkerCount(2), 1);
	EXPECT_EQ(defaultLovWorkerCount(16), 15);
}
