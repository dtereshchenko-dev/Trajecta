#include <gtest/gtest.h>

#include "hazard/LovImpactProbability.h"
#include "hazard/TargetPlaneTrace.h"

#include <algorithm>
#include <array>
#include <cmath>

using namespace orbital_core::hazard;

namespace {
	TargetPlaneTracePoint makeTube(double xKm, double yKm, double widthBodyRadii, double stretchBodyRadii,
	                               double alphaRad, double sigma)
	{
		TargetPlaneTracePoint t;
		t.coordsKm = Eigen::Vector2d(xKm, yKm);
		t.bodyRadiusKm = 1.0;
		t.widthBodyRadii = widthBodyRadii;
		t.stretchBodyRadii = stretchBodyRadii;
		t.alphaRad = alphaRad;
		t.sigma = sigma;
		t.valid = true;
		return t;
	}

	LovImpactProbabilityResult run(const TargetPlaneTracePoint& t, double r = 1.0)
	{
		LovImpactProbabilityInput in;
		in.tracePoint = t;
		in.impactCrossSectionBodyRadii = r;
		return lovImpactProbability(in);
	}

	constexpr double kNativeCsi = 1.8589562717887398;
	constexpr double kNativeZeta = -2.0481422973853115e-2;
	constexpr double kNativeSigma = -1.0697174233201485;
	constexpr double kNativeWidth = 4.3403076042778022e-4;
	constexpr double kNativeStretch = 2117.8613077568189;
	constexpr double kNativeCrossSection = 4.4904321008221268;
	constexpr double kNativeProbability = 8.6902472193722607e-4;

	constexpr double kProbabilityAbsTolerance = 1.0e-8;
	constexpr double kProbabilityRelTolerance = 1.0e-4;
}  // namespace

TEST(LovImpactProbabilityTest, DegenerateInputsInvalid)
{
	EXPECT_FALSE(run(makeTube(0, 0, 0.0, 1.0, 0.0, 0.0)).valid);
	EXPECT_FALSE(run(makeTube(0, 0, 1.0, 0.0, 0.0, 0.0)).valid);
	EXPECT_FALSE(run(makeTube(0, 0, 1.0, 1.0, 0.0, 0.0), 0.0).valid);

	TargetPlaneTracePoint invalid = makeTube(0, 0, 1.0, 1.0, 0.0, 0.0);
	invalid.valid = false;
	EXPECT_FALSE(run(invalid).valid);
}

TEST(LovImpactProbabilityTest, FarTubeNotCandidate)
{
	const auto res = run(makeTube(50.0, 0.0, 1.0, 1.0, 0.0, 0.0));
	EXPECT_TRUE(res.valid);
	EXPECT_FALSE(res.candidate);
	EXPECT_DOUBLE_EQ(res.probability, 0.0);
}

TEST(LovImpactProbabilityTest, Rejects2044PointOutsideRadialCorridor)
{
	constexpr double b = 23.039978021869469;
	constexpr double csi = 0.904003701287679;
	constexpr double width = 0.00045212202172250569;
	constexpr double radius = 4.575985398886724;
	TargetPlaneTracePoint tube = makeTube(-146950.58933924171, 674.28193078128345, width, 3806.2594385780631,
	                                      1.6054542482405858, -3.4717894768365274);
	tube.bodyRadiusKm = 6378.137;

	ASSERT_NEAR(tube.coordsKm.norm() / tube.bodyRadiusKm, b, 1e-13);
	ASSERT_NEAR(
	    (tube.coordsKm.x() * std::cos(tube.alphaRad) + tube.coordsKm.y() * std::sin(tube.alphaRad)) / tube.bodyRadiusKm,
	    csi, 1e-14);
	ASSERT_LT(std::abs(csi) - 10.0 * width, radius);
	ASSERT_GT(b - 10.0 * width, radius);

	const auto result = run(tube, radius);
	ASSERT_TRUE(result.valid);
	EXPECT_FALSE(result.candidate);
	EXPECT_DOUBLE_EQ(result.probability, 0.0);
	EXPECT_DOUBLE_EQ(result.sigmaImpact, 0.0);
}

TEST(LovImpactProbabilityTest, CandidateGateIsRotationInvariant)
{
	constexpr double bodyRadiusKm = 6378.137;
	for (const double b : {1.5, 3.0}) {
		for (const double angle : {0.0, 0.7, 1.5707963267948966}) {
			for (const double alpha : {0.0, angle}) {
				SCOPED_TRACE(b);
				SCOPED_TRACE(angle);
				SCOPED_TRACE(alpha);
				TargetPlaneTracePoint tube = makeTube(bodyRadiusKm * b * std::cos(angle),
				                                      bodyRadiusKm * b * std::sin(angle), 0.125, 10.0, alpha, 0.0);
				tube.bodyRadiusKm = bodyRadiusKm;

				const auto result = run(tube);
				ASSERT_TRUE(result.valid);
				EXPECT_EQ(result.candidate, b < 2.25);
			}
		}
	}
}

TEST(LovImpactProbabilityTest, RadialCandidateBoundaryIsStrict)
{
	constexpr double boundary = 2.25;
	const auto below = run(makeTube(std::nextafter(boundary, 0.0), 0.0, 0.125, 10.0, 0.0, 0.0));
	const auto at = run(makeTube(boundary, 0.0, 0.125, 10.0, 0.0, 0.0));
	const auto above = run(makeTube(std::nextafter(boundary, 3.0), 0.0, 0.125, 10.0, 0.0, 0.0));

	ASSERT_TRUE(below.valid);
	ASSERT_TRUE(at.valid);
	ASSERT_TRUE(above.valid);
	EXPECT_TRUE(below.candidate);
	EXPECT_FALSE(at.candidate);
	EXPECT_FALSE(above.candidate);
	EXPECT_DOUBLE_EQ(at.probability, 0.0);
	EXPECT_DOUBLE_EQ(above.probability, 0.0);
	EXPECT_DOUBLE_EQ(at.sigmaImpact, 10.0);
}

TEST(LovImpactProbabilityTest, CenteredTubeApproachesGaussianMass)
{
	const auto res = run(makeTube(0.0, 0.0, 1.0, 0.01, 0.0, 0.0), 1.0);
	const double gaussianMass = std::erf(1.0 / std::sqrt(2.0));

	EXPECT_TRUE(res.valid);
	EXPECT_TRUE(res.candidate);
	EXPECT_GE(res.probability, 0.0);
	EXPECT_LE(res.probability, 1.0);
	EXPECT_LE(res.probability, gaussianMass + 1.0e-9);
	EXPECT_NEAR(res.probability, gaussianMass, 0.03);
}

TEST(LovImpactProbabilityTest, OffsetReducesProbability)
{
	const auto centered = run(makeTube(0.0, 0.0, 1.0, 0.01, 0.0, 0.0));
	const auto offset = run(makeTube(0.6, 0.0, 1.0, 0.01, 0.0, 0.0));

	EXPECT_TRUE(centered.candidate);
	EXPECT_TRUE(offset.candidate);
	EXPECT_GT(centered.probability, offset.probability);
	EXPECT_GE(offset.probability, 0.0);
}

TEST(LovImpactProbabilityTest, FocusedCrossSectionIncreasesProbability)
{
	const TargetPlaneTracePoint tube = makeTube(0.0, 0.0, 0.1, 100.0, 0.0, 0.0);
	const auto physical = run(tube, 1.0);
	const auto focused = run(tube, 4.0);

	ASSERT_TRUE(physical.valid);
	ASSERT_TRUE(focused.valid);
	EXPECT_GT(focused.probability, physical.probability);
}

TEST(LovImpactProbabilityTest, MatchesFrozen2020Vv2050Probability)
{
	const auto result =
	    run(makeTube(kNativeCsi, kNativeZeta, kNativeWidth, kNativeStretch, 0.0, kNativeSigma), kNativeCrossSection);

	ASSERT_TRUE(result.valid);
	ASSERT_TRUE(result.candidate);
	EXPECT_DOUBLE_EQ(result.sigmaImpact, 0.0);
	EXPECT_NEAR(result.probability, kNativeProbability,
	            std::max(kProbabilityAbsTolerance, kProbabilityRelTolerance * kNativeProbability));
}

TEST(LovImpactProbabilityTest, ThinNativeTubeHasInverseStretchingDependence)
{
	const TargetPlaneTracePoint tube =
	    makeTube(kNativeCsi, kNativeZeta, kNativeWidth, kNativeStretch, 0.0, kNativeSigma);
	const auto reference = run(tube, kNativeCrossSection);
	ASSERT_TRUE(reference.valid);
	ASSERT_TRUE(reference.candidate);
	ASSERT_GT(reference.probability, 0.0);

	const std::array<double, 2> stretchings{2164.039225, 5731.263902};
	for (const double stretch : stretchings) {
		SCOPED_TRACE(stretch);
		TargetPlaneTracePoint changed = tube;
		changed.stretchBodyRadii = stretch;
		const auto result = run(changed, kNativeCrossSection);
		const double expected = reference.probability * kNativeStretch / stretch;

		ASSERT_TRUE(result.valid);
		ASSERT_TRUE(result.candidate);
		EXPECT_LT(result.probability, reference.probability);
		EXPECT_NEAR(result.probability, expected,
		            std::max(kProbabilityAbsTolerance, kProbabilityRelTolerance * expected));
	}
}

TEST(LovImpactProbabilityTest, MatchesReferenceSmallProbabilityQuadrature)
{
	struct Fixture {
		double csi, zeta, sigma, width, stretch, radius, probability;
	};

	const std::array<Fixture, 6> cases{
	    {{0.93053855629059179, 1.3397132344451865, -3.4776921139207539, 7.0746161027544145e-4, 207.09952758752726,
	      4.4813635166743904, 3.9077720386738237e-5},
	     {-1.7475497073567881, 0.72600347794801756, -3.4774947019603770, 8.6813686903645655e-4, 4405.3110234297919,
	      4.4844806192995632, 1.7689152939846604e-6},
	     {-3.2202266866819342, -0.64968581157805050, -3.4879994395824796, 1e-4, 1113298.5505912453, 4.5334250841583801,
	      5.2167961401635934e-9},
	     {-1.2899589903523496, 0.16146502018434988, -3.5338329186540172, 1e-4, 1103296.7848397670, 4.4506933292994688,
	      5.9826385703271673e-9},
	     {-1.6870700217156451, 0.32564861331691852, -3.4304248310579779, 1e-4, 732230.90949170536, 4.4878645579615730,
	      1.2615361250571082e-8},
	     {1.8555488421000155, -0.33955476550068786, -1.0702017200000000, 4.1630768058938778e-4, 5658.2818984398991,
	      4.4890104723316604, 3.2511872260630227e-4}}};
	for (const auto& f : cases) {
		SCOPED_TRACE(f.probability);
		const auto result = run(makeTube(f.csi, f.zeta, f.width, f.stretch, 0.0, f.sigma), f.radius);
		ASSERT_TRUE(result.valid);
		EXPECT_NEAR(result.probability, f.probability, f.probability * 1e-6);
	}
}

TEST(LovImpactProbabilityTest, ResolvesSmallProbabilitiesInBothGaussianTails)
{
	constexpr double radius = 4.5;
	constexpr double stretch = 1e6;
	constexpr double sigma = 5.0;
	const double narrowTubeMass =
	    2 * radius / stretch * std::exp(-sigma * sigma / 2) / std::sqrt(2 * 3.14159265358979323846);
	for (double sign : {-1.0, 1.0}) {
		const auto result = run(makeTube(0, 0, 1e-4, stretch, 0, sign * sigma), radius);
		ASSERT_TRUE(result.valid);
		EXPECT_GT(result.probability, 1e-11);
		EXPECT_NEAR(result.probability, narrowTubeMass, narrowTubeMass * 1e-6);
	}
}
