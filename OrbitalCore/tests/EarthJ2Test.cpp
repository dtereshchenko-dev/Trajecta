#include <gtest/gtest.h>
#include "NBodyPropagator.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Types.h"
#include "orbital_core/Constants.h"
#include "orbital_core/forces/EarthJ2.h"
#include "orbital_core/TimeScale.h"
#include "earth_orientation/EarthOrientation.h"
#include "earth_orientation/EarthRotation.h"
#include <cmath>
#include <filesystem>

class EarthJ2Test : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (getExecutableDir().empty()) {
			auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
			std::filesystem::path candidates[] = {
			    srcDir.parent_path() / "data",
			    srcDir / "data",
			    srcDir / "release",
			    srcDir / "debug",
			};
			for (const auto& dir : candidates) {
				if (std::filesystem::exists(dir / "linux_p1550p2650.440")) {
					getExecutableDir() = dir;
					break;
				}
			}
		}
		m_de440 = std::make_unique<DE440Source>();
		if (m_de440->isAvailable())
			setDefaultEphemerisSource(m_de440.get());
	}

	void TearDown() override { setDefaultEphemerisSource(nullptr); }

	static constexpr double kApophis2029JD = 2462240.4;

	std::unique_ptr<DE440Source> m_de440;
	NBodyPropagator m_prop;
};

TEST(NutationReferenceModel, DynamicalPoleMatchesRotsys2044)
{
	const Vec3D pole =
	    Helpers::ecvToEcl(EarthRotation::Iau1976_1980::dynamicalPoleEquatorial(2400000.5 + 67898.436647806448));
	EXPECT_NEAR(pole.x, 0.00436009458995341350, 2e-12);
	EXPECT_NEAR(pole.y, 0.397792956747345217, 2e-12);
	EXPECT_NEAR(pole.z, 0.917464851172714591, 2e-12);
}

TEST_F(EarthJ2Test, RotatingPoleAndZonalForceMatchPotential)
{
	if (!m_de440->isAvailable() || !orbital_core::earthOrientation().isAvailable())
		GTEST_SKIP() << "DE440 or C04 data unavailable";

	orbital_core::forces::ForceContext ctx;
	ctx.R = {2.5e-5, -1.7e-5, 2.0e-5};
	ctx.earthGm = 8.8876924467071022e-10;
	ctx.ephemerisSource = m_de440.get();
	ctx.havePlanets = true;
	ctx.planetPositions[2] = {};
	ctx.bodyMask = 1u << 2;
	ctx.computePartials = true;
	Eigen::Vector3d firstAcceleration;

	for (const double mjd : {60000.0, 61000.0}) {
		ctx.setEpochMjd(mjd);
		const double tt = orbital_core::time::fromMjdTdb(mjd, orbital_core::time::Scale::TT);
		const auto eop = orbital_core::earthOrientation().lookupMjdTt(tt);
		ASSERT_TRUE(eop.valid);
		const Vec3D p = Helpers::ecvToEcl(EarthRotation::WithEop::poleEquatorialMjd(tt, eop));
		const Eigen::Vector3d axis = Eigen::Vector3d(p.x, p.y, p.z).normalized();
		const Eigen::Vector3d r(ctx.R.x, ctx.R.y, ctx.R.z);
		const auto potential = [&](const Eigen::Vector3d& pos) {
			const double distance = pos.norm(), s = pos.dot(axis) / distance;
			const double rho = (6378.1366 / 149597870.7) / distance;
			const double p2 = (3 * s * s - 1) / 2;
			const double p3 = (5 * s * s * s - 3 * s) / 2;
			const double p4 = (35 * s * s * s * s - 30 * s * s + 3) / 8;
			return -ctx.earthGm / distance *
			       (1.08262539e-3 * rho * rho * p2 - 2.532410e-6 * rho * rho * rho * p3 -
			        1.619898e-6 * rho * rho * rho * rho * p4);
		};
		Vec3D acceleration;
		Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero(), dadv = Eigen::Matrix3d::Zero();
		orbital_core::forces::EarthJ2{}.evaluate(ctx, acceleration, dadr, dadv);
		const Eigen::Vector3d actual(acceleration.x, acceleration.y, acceleration.z);
		Eigen::Vector3d expected;
		const double step = r.norm() * 1e-4;
		for (int j = 0; j < 3; ++j) {
			Eigen::Vector3d d = Eigen::Vector3d::Zero();
			d[j] = step;
			expected[j] = (potential(r - 2 * d) - 8 * potential(r - d) + 8 * potential(r + d) - potential(r + 2 * d)) /
			              (12 * step);
		}
		EXPECT_LT((actual - expected).norm() / expected.norm(), 1e-8);
		EXPECT_DOUBLE_EQ(dadv.norm(), 0.0);
		if (mjd == 60000.0)
			firstAcceleration = actual;
		else
			EXPECT_GT((actual - firstAcceleration).norm() / actual.norm(), 1e-6);
	}
}

TEST_F(EarthJ2Test, BothCodePathsReturnGeocenter)
{
	if (!m_de440->isAvailable())
		GTEST_SKIP() << "DE440 kernel not present — skipping";

	const double jd = kApophis2029JD;

	Vec3D planets[8];
	Vec3D moonHelio;
	m_de440->getAllPlanetPositions(jd, planets, &moonHelio);
	const Vec3D rEarthBatch = planets[2];

	const Vec3D rEarthSingle = m_de440->getPlanetPosition(EphemerisPlanet::Earth, jd);

	const Vec3D diff = rEarthBatch - rEarthSingle;
	const double diffKm = diff.magnitude() * 149597870.7;
	EXPECT_LT(diffKm, 1.0) << "Hot/cold path Earth-position mismatch: " << diffKm
	                       << " km — one path is returning EMB instead of geocenter";
}

TEST_F(EarthJ2Test, EarthBodyMaskControlsTerm)
{
	if (!m_de440->isAvailable())
		GTEST_SKIP() << "DE440 kernel not present - skipping";

	namespace ocf = orbital_core::forces;
	const double jd = kApophis2029JD;
	const Vec3D rEarth = m_de440->getPlanetPosition(EphemerisPlanet::Earth, jd);

	ocf::ForceContext ctx;
	ctx.R = rEarth + Vec3D(6.0 * Constants::EARTH_EQ_RADIUS_AU, 0.0, 0.0);
	ctx.V = Vec3D(0.0, 0.0172, 0.0);
	ctx.setEpochMjd(jd - 2400000.5);
	ctx.ephemerisSource = m_de440.get();
	ctx.planetPositions[2] = rEarth;
	ctx.havePlanets = true;
	ctx.computePartials = false;

	Vec3D disabledAcceleration;
	Eigen::Matrix3d dAdR = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dAdV = Eigen::Matrix3d::Zero();
	ocf::EarthJ2{}.evaluate(ctx, disabledAcceleration, dAdR, dAdV);
	EXPECT_DOUBLE_EQ(disabledAcceleration.magnitude(), 0.0);

	ctx.bodyMask = static_cast<std::uint16_t>(1u << 2);
	Vec3D enabledAcceleration;
	ocf::EarthJ2{}.evaluate(ctx, enabledAcceleration, dAdR, dAdV);
	EXPECT_GT(enabledAcceleration.magnitude(), 0.0);
}

TEST_F(EarthJ2Test, PartialMatchesFiniteDifference)
{
	if (!m_de440->isAvailable())
		GTEST_SKIP() << "DE440 kernel not present — skipping";

	const double jd = kApophis2029JD;
	const Vec3D rEarth = m_de440->getPlanetPosition(EphemerisPlanet::Earth, jd);

	const double h = 4.0 * Constants::EARTH_EQ_RADIUS_AU;
	const Vec3D R(rEarth.x + h * 0.7, rEarth.y + h * 0.5, rEarth.z + h * 0.5);
	const Vec3D V(0.0, 0.0172, 0.0);

	const Eigen::Matrix3d analytic = m_prop.accelerationPartialR(R, V, jd);

	const double step = h * 1e-6;
	Eigen::Matrix3d fd;
	for (int j = 0; j < 3; ++j) {
		Vec3D Rp = R, Rm = R;
		switch (j) {
			case 0:
				Rp.x += step;
				Rm.x -= step;
				break;
			case 1:
				Rp.y += step;
				Rm.y -= step;
				break;
			case 2:
				Rp.z += step;
				Rm.z -= step;
				break;
		}
		const Vec3D ap = m_prop.accelerationAt(Rp, V, jd);
		const Vec3D am = m_prop.accelerationAt(Rm, V, jd);
		fd(0, j) = (ap.x - am.x) / (2.0 * step);
		fd(1, j) = (ap.y - am.y) / (2.0 * step);
		fd(2, j) = (ap.z - am.z) / (2.0 * step);
	}

	const double tol = 1.0e-5;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			const double scale = std::max(1.0e-30, std::abs(analytic(i, j)) + std::abs(fd(i, j)));
			EXPECT_NEAR(analytic(i, j), fd(i, j), tol * scale) << "mismatch at (" << i << ", " << j << ")";
		}
	}
}
