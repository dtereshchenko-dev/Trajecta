#include <gtest/gtest.h>
#include "obs/AstrometricResidual.h"
#include "OrbitOptimizer.h"
#include "orbital_core/Types.h"
#include <cmath>

TEST(OrbitOptimizerConfig, Defaults)
{
	const OrbitOptimizerConfig cfg;
	EXPECT_FALSE(cfg.useWeights);
	EXPECT_FALSE(cfg.evaluateOnly);
}

TEST(AstrometricAngularJacobian, MatchesCentralDifferencesAcrossRaWrap)
{
	const double c = std::cos(Constants::EARTH_INCL);
	const double s = std::sin(Constants::EARTH_INCL);
	const double metric = std::cos(0.3);
	for (const Vec3D eq : {Vec3D(0.04, 0.02, 0.01), Vec3D(1.0, 0.0, 0.2), Vec3D(-0.5, -0.2, -0.1)}) {
		const auto analytic = orbital_core::obs::astrometricAngularJacobianEcl(eq, metric);
		ASSERT_TRUE(analytic.allFinite());
		const Vec3D ecl(eq.x, c * eq.y + s * eq.z, -s * eq.y + c * eq.z);
		const auto angles = [=](const Vec3D& p) {
			const Vec3D q(p.x, c * p.y - s * p.z, s * p.y + c * p.z);
			return Eigen::Vector2d(std::atan2(q.y, q.x), std::asin(q.z / q.magnitude()));
		};
		constexpr double h = 1e-7;
		for (int k = 0; k < 3; ++k) {
			Vec3D plus = ecl, minus = ecl;
			(k == 0 ? plus.x : k == 1 ? plus.y : plus.z) += h;
			(k == 0 ? minus.x : k == 1 ? minus.y : minus.z) -= h;
			const auto a = angles(plus), b = angles(minus);
			const double dra = std::remainder(a[0] - b[0], 2.0 * Constants::PI);
			EXPECT_NEAR(analytic(0, k), -metric * dra / (2 * h), 2e-8);
			EXPECT_NEAR(analytic(1, k), -(a[1] - b[1]) / (2 * h), 2e-8);
		}
	}
}

TEST(AstrometricAngularJacobian, RejectsSingularDirections)
{
	EXPECT_FALSE(orbital_core::obs::astrometricAngularJacobianEcl(Vec3D(0, 0, 0), 1).allFinite());
	EXPECT_FALSE(orbital_core::obs::astrometricAngularJacobianEcl(Vec3D(0, 0, 1), 0).allFinite());
	EXPECT_FALSE(orbital_core::obs::astrometricAngularJacobianEcl(Vec3D(0, 0, -1), 0).allFinite());
}
