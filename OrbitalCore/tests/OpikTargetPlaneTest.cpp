#include <gtest/gtest.h>

#include "hazard/OpikTargetPlane.h"
#include "hazard/LovImpactProbability.h"

#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/Planet.h"
#include "orbital_core/hazard/PlanetData.h"

#include <cmath>

using namespace orbital_core;
using namespace orbital_core::hazard;

TEST(OpikTargetPlaneTest, FocusingVInfAndOrthonormalBasis)
{
	const Vec3D rGeoKm(50000.0, 10000.0, 2000.0);
	const Vec3D vGeoKmS(-9.0, 4.0, 1.0);
	const auto pt = computeOpikTargetPlanePoint(rGeoKm, vGeoKmS, Planet::Earth);

	ASSERT_TRUE(pt.valid);

	EXPECT_TRUE(std::isfinite(pt.gamma));
	EXPECT_GE(pt.gamma, 1.0);

	EXPECT_GT(pt.vInfKmS, 0.0);
	EXPECT_LT(pt.vInfKmS, vGeoKmS.magnitude());

	EXPECT_TRUE(pt.basis.valid);
	EXPECT_NEAR(pt.basis.xi.magnitude(), 1.0, 1.0e-9);
	EXPECT_NEAR(pt.basis.zeta.magnitude(), 1.0, 1.0e-9);
	EXPECT_NEAR(pt.basis.s.magnitude(), 1.0, 1.0e-9);
	EXPECT_NEAR(pt.basis.xi.dot(pt.basis.zeta), 0.0, 1.0e-9);
	EXPECT_NEAR(pt.basis.xi.dot(pt.basis.s), 0.0, 1.0e-9);
	EXPECT_NEAR(pt.basis.zeta.dot(pt.basis.s), 0.0, 1.0e-9);

	const double v2 = vGeoKmS.dot(vGeoKmS);
	const Vec3D rMtp = rGeoKm - vGeoKmS * (rGeoKm.dot(vGeoKmS) / v2);
	const double d = rMtp.magnitude();
	const double mu = getPlanetData(Planet::Earth).muKm3PerS2;
	const double bOverD = std::sqrt(d * v2 / (d * v2 - 2.0 * mu));
	EXPECT_NEAR(pt.coordsKm.norm(), d * bOverD, d * bOverD * 1.0e-12);
}

TEST(OpikTargetPlaneTest, DegenerateGeometryInvalid)
{
	const auto pt = computeOpikTargetPlanePoint(Vec3D(50000.0, 0.0, 0.0), Vec3D(0.0, 0.0, 0.0), Planet::Earth);
	EXPECT_FALSE(pt.valid);
}

TEST(OpikTargetPlaneTest, EarthImpactDiskUsesOnePhysicalBoundaryForReferenceFixtures)
{
	constexpr double auKm = 149597870.7;

	struct Fixture {
		Vec3D rAu;
		Vec3D vAuDay;
		double diskBodyRadii;
		double sigma;
		Eigen::Vector2d coordsNativeRadii;
		double stretchNativeRadii;
		double widthNativeRadii;
		double alpha;
		double probability;
	};

	const Fixture fixtures[]{{{2.5222450502984935e-6, 5.321308513493716e-6, 1.416022231945161e-6},
	                          {0.004874284021826811, 0.002044015599826263, -0.016363426956257237},
	                          4.510777014832009,
	                          -3.477692113920754,
	                          {-0.3262157185253612, 1.5982230320821762},
	                          207.09952758752726,
	                          0.0007074616102754415,
	                          0.8084251736687272,
	                          3.907772038673824e-5},
	                         {{2.792377548358555e-6, 7.396834224467981e-6, 1.8460002291090134e-6},
	                          {0.009145381770562944, -0.0005289429856361285, -0.01171443795474426},
	                          4.513914576608060,
	                          -3.477494701960377,
	                          {-1.0858337090758456, 1.5498310185016961},
	                          4405.311023429792,
	                          0.0008681368690364565,
	                          -0.5659120734312353,
	                          1.7689152939846604e-6}};
	for (const auto& fixture : fixtures) {
		const auto point =
		    computeOpikTargetPlanePoint(fixture.rAu * auKm, fixture.vAuDay * (auKm / 86400.0), Planet::Earth);
		ASSERT_TRUE(point.valid);
		const double expectedDisk =
		    6420.0 / planetRadiusKm(Planet::Earth) *
		    std::sqrt(1.0 + 2.0 * planetMuKm3PerS2(Planet::Earth) / (6420.0 * point.vInfKmS * point.vInfKmS));
		EXPECT_NEAR(point.gamma, expectedDisk, 2.0e-12);
		EXPECT_LT(point.gamma, fixture.diskBodyRadii);
		LovImpactProbabilityInput input;
		auto& trace = input.tracePoint;
		trace.valid = true;
		trace.bodyRadiusKm = planetRadiusKm(Planet::Earth);
		trace.coordsKm = fixture.coordsNativeRadii * 6420.0;
		trace.sigma = fixture.sigma;
		trace.stretchBodyRadii = fixture.stretchNativeRadii * 6420.0 / trace.bodyRadiusKm;
		trace.widthBodyRadii = fixture.widthNativeRadii * 6420.0 / trace.bodyRadiusKm;
		trace.alphaRad = fixture.alpha;
		input.impactCrossSectionBodyRadii = point.gamma;
		const auto probability = lovImpactProbability(input);
		ASSERT_TRUE(probability.valid);
		ASSERT_TRUE(probability.candidate);
		EXPECT_GT(probability.probability, 0.0);
		EXPECT_LT(probability.probability, fixture.probability);
	}
	EXPECT_DOUBLE_EQ(planetRadiusKm(Planet::Earth), 6378.137);
	const double expectedPhysicalFocusing =
	    std::sqrt(1.0 + 2.0 * planetMuKm3PerS2(Planet::Earth) / (6378.137 * 3.0 * 3.0));
	EXPECT_DOUBLE_EQ(gravitationalFocusing(3.0, Planet::Earth), expectedPhysicalFocusing);
}

TEST(OpikTargetPlaneTest, EarthImpactDiskAgreesWithPhysicalGrazingPericenter)
{
	constexpr double boundaryKm = 6420.0;
	const double mu = planetMuKm3PerS2(Planet::Earth);
	for (double vInf : {2.5, 10.0, 40.0}) {
		for (double radius : {boundaryKm - 1.0, boundaryKm, boundaryKm + 1.0}) {
			const double speed = std::sqrt(vInf * vInf + 2.0 * mu / radius);
			const auto point = computeOpikTargetPlanePoint({radius, 0.0, 0.0}, {0.0, speed, 0.0}, Planet::Earth);
			ASSERT_TRUE(point.valid);
			const double diskKm = point.gamma * planetRadiusKm(Planet::Earth);
			const double expectedKm = std::sqrt(boundaryKm * boundaryKm + 2.0 * mu * boundaryKm / (vInf * vInf));
			EXPECT_NEAR(diskKm, expectedKm, expectedKm * 1e-12);
			if (radius == boundaryKm)
				EXPECT_NEAR(point.coordsKm.norm(), diskKm, diskKm * 1e-12);
			else if (radius < boundaryKm)
				EXPECT_LT(point.coordsKm.norm(), diskKm);
			else
				EXPECT_GT(point.coordsKm.norm(), diskKm);
		}
	}
}

TEST(OpikTargetPlaneTest, OtherPlanetsRetainPhysicalImpactDisk)
{
	const auto point = computeOpikTargetPlanePoint({50000.0, 10000.0, 2000.0}, {-9.0, 4.0, 1.0}, Planet::Mars);
	ASSERT_TRUE(point.valid);
	EXPECT_DOUBLE_EQ(point.gamma, gravitationalFocusing(point.vInfKmS, Planet::Mars));
}

TEST(OpikTargetPlaneTest, ReproducesReferenceAnalyticJacobianFixture)
{
	constexpr double auKm = 149597870.7;
	constexpr double secondsPerDay = 86400.0;
	const Vec3D rAu(2.5410791863263239e-6, 5.8147624856710678e-6, 1.5105108124139222e-6);
	const Vec3D vAuDay(7.0364139085299215e-3, 8.1759230879430720e-4, -1.4984460772526690e-2);
	const auto differential =
	    computeOpikTargetPlaneDifferential(rAu * auKm, vAuDay * (auKm / secondsPerDay), Planet::Earth);
	ASSERT_TRUE(differential.valid);

	const Eigen::Vector2d expectedCoordsAu(-2.6168555933143281e-5, 6.7810998771159042e-5);
	const Eigen::Vector2d coordsAu = differential.point.coordsKm / auKm;
	for (int row = 0; row < 2; ++row)
		EXPECT_NEAR(coordsAu(row), expectedCoordsAu(row), 1.0e-10) << "coordinate row " << row;

	Eigen::Matrix<double, 2, 6> jacobianAu;
	jacobianAu.block<2, 3>(0, 0) = differential.dCoordsDGeoState.block<2, 3>(0, 0);
	jacobianAu.block<2, 3>(0, 3) = differential.dCoordsDGeoState.block<2, 3>(0, 3) / secondsPerDay;

	Eigen::Matrix<double, 2, 6> expected;
	expected << 9.7782176292445371e1, 2.1797479555897394e2, 5.7809913535138350e1, 8.0179786277276344e-2,
	    1.1404003794850470e-2, -1.7921018289757132e-1, -2.4413266019415053e2, -5.6088057266367787e2,
	    -1.4524313694815169e2, -2.1483418191176068e-1, -2.4156944406678763e-2, 4.5423638992696064e-1;
	for (int row = 0; row < 2; ++row) {
		for (int col = 0; col < 6; ++col) {
			EXPECT_NEAR(jacobianAu(row, col), expected(row, col), 1.0e-8)
			    << "Jacobian row " << row << ", column " << col;
		}
	}
}
