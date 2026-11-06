#include <gtest/gtest.h>
#include "orbital_core/DE440Source.h"
#include "orbital_core/Types.h"
#include "orbital_core/forces/EihRelativity.h"
#include "orbital_core/forces/EarthJ2.h"
#include "orbital_core/forces/SunJ2.h"
#include "orbital_core/forces/SunGravity.h"
#include "orbital_core/perturbations/PlanetaryPerturbations.h"
#include <atomic>
#include <cmath>
#include <filesystem>
#include <thread>
#include <vector>

namespace {

	Vec3D legacySolarJ2(const Vec3D& r, double gm)
	{
		const double r2 = r.dot(r), z2 = r.z * r.z;
		const double factor = 1.5 * gm * 2e-7 * 4.6527174e-3 * 4.6527174e-3 / std::pow(r.magnitude(), 7);
		return Vec3D(r.x * (5 * z2 - r2), r.y * (5 * z2 - r2), r.z * (5 * z2 - 3 * r2)) * factor;
	}

	Vec3D reconstructedSolarVelocity(const orbital_core::forces::ForceContext& c)
	{
		double mass = c.sunGm + c.moonGm;
		Vec3D momentum = c.moonVelocity * c.moonGm;
		for (int i = 0; i < 8; ++i) {
			mass += c.planetGms[i];
			momentum += c.planetVelocities[i] * c.planetGms[i];
		}
		return momentum * (-1.0 / mass);
	}

}  // namespace

class DE440SourceTest : public ::testing::Test {
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
		source = std::make_unique<DE440Source>();
	}

	std::unique_ptr<DE440Source> source;
};

TEST_F(DE440SourceTest, RecreatedSourcesDoNotReuseFreedMappedContexts)
{
	if (!source->isAvailable())
		GTEST_SKIP();
	const auto expected = source->queryJD(EphemerisPlanet::Earth, 2460000.5);
	source.reset();
	for (int i = 0; i < 16; ++i) {
		DE440Source recreated;
		ASSERT_TRUE(recreated.isAvailable());
		const auto actual = recreated.queryJD(EphemerisPlanet::Earth, 2460000.5);
		EXPECT_NEAR((actual.R - expected.R).magnitude(), 0, 1e-15);
		Vec3D positions[8], velocities[8];
		recreated.getAllPlanetStates(2460000.5, positions, velocities);
		EXPECT_NEAR((positions[2] - expected.R).magnitude(), 0, 1e-15);
	}
}

TEST_F(DE440SourceTest, IsAvailableWithFile)
{
	if (!source->isAvailable())
		GTEST_SKIP() << "DE440 file not found";
	EXPECT_TRUE(source->isAvailable());
	EXPECT_EQ(source->name(), "DE440");
}

TEST_F(DE440SourceTest, MjdQueryPreservesFramesAndSubJdTimes)
{
	if (!source->isAvailable())
		GTEST_SKIP();
	constexpr double mjd = 59139.25;
	for (int i = 0; i <= 9; ++i) {
		const auto planet = static_cast<EphemerisPlanet>(i);
		const auto jdState = source->queryJD(planet, mjd + 2400000.5);
		const auto mjdState = source->queryMjd(planet, mjd);
		SCOPED_TRACE(i);
		EXPECT_LT((jdState.R - mjdState.R).magnitude(), 2e-14);
		EXPECT_LT((jdState.V - mjdState.V).magnitude(), 2e-16);
	}
	const double later = mjd + 1e-10;
	ASSERT_EQ(mjd + 2400000.5, later + 2400000.5);
	const auto a = source->queryMjd(EphemerisPlanet::Earth, mjd);
	const auto b = source->queryMjd(EphemerisPlanet::Earth, later);
	EXPECT_GT((b.R - a.R).magnitude(), 1e-12);
	EXPECT_LT((b.R - a.R - a.V * (later - mjd)).magnitude(), 5e-16);
}

TEST_F(DE440SourceTest, EarthPositionJ2000Epoch)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	Vec3D pos = source->getPlanetPosition(EphemerisPlanet::Earth, Constants::J2000_JD);

	double dist = pos.magnitude();
	EXPECT_NEAR(dist, 0.983, 0.01) << "Earth heliocentric distance at J2000";
}

TEST_F(DE440SourceTest, SunPositionIsZero)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	Vec3D pos = source->getPlanetPosition(EphemerisPlanet::Sun, Constants::J2000_JD);
	EXPECT_DOUBLE_EQ(pos.x, 0.0);
	EXPECT_DOUBLE_EQ(pos.y, 0.0);
	EXPECT_DOUBLE_EQ(pos.z, 0.0);
}

TEST_F(DE440SourceTest, JupiterDistanceReasonable)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	Vec3D pos = source->getPlanetPosition(EphemerisPlanet::Jupiter, 2460000.5);
	double dist = pos.magnitude();

	EXPECT_GT(dist, 4.9);
	EXPECT_LT(dist, 5.5);
}

TEST_F(DE440SourceTest, MarsDistanceReasonable)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	Vec3D pos = source->getPlanetPosition(EphemerisPlanet::Mars, 2460000.5);
	double dist = pos.magnitude();

	EXPECT_GT(dist, 1.3);
	EXPECT_LT(dist, 1.7);
}

TEST_F(DE440SourceTest, FetchWindowReturnsCorrectCount)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	auto records = source->fetchWindow(EphemerisPlanet::Earth, 2460000.5, 5);

	EXPECT_EQ(records.size(), 11u);

	for (size_t i = 1; i < records.size(); ++i)
		EXPECT_NEAR(records[i].jd - records[i - 1].jd, 1.0, 1e-10);
}

TEST_F(DE440SourceTest, EarthVelocityReasonable)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	auto records = source->fetchWindow(EphemerisPlanet::Earth, 2460000.5, 1);
	ASSERT_GE(records.size(), 2u);

	double vMag = records[1].V.magnitude();
	EXPECT_NEAR(vMag, 0.0172, 0.001);
}

TEST_F(DE440SourceTest, BatchedStatesMatchSingleBodyQueries)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	constexpr double jd = 2460000.5;
	Vec3D positions[8];
	Vec3D velocities[8];
	Vec3D moonPosition;
	Vec3D moonVelocity;
	source->getAllPlanetStates(jd, positions, velocities, &moonPosition, &moonVelocity);

	constexpr EphemerisPlanet bodies[8] = {EphemerisPlanet::Mercury, EphemerisPlanet::Venus,   EphemerisPlanet::Earth,
	                                       EphemerisPlanet::Mars,    EphemerisPlanet::Jupiter, EphemerisPlanet::Saturn,
	                                       EphemerisPlanet::Uranus,  EphemerisPlanet::Neptune};
	for (int index = 0; index < 8; ++index) {
		const EphemerisRecord single = source->queryJD(bodies[index], jd);
		EXPECT_LT((positions[index] - single.R).magnitude(), 1.0e-12);
		EXPECT_LT((velocities[index] - single.V).magnitude(), 1.0e-12);
	}

	const EphemerisRecord moon = source->queryJD(EphemerisPlanet::Moon, jd);
	EXPECT_LT((moonPosition - moon.R).magnitude(), 1.0e-12);
	EXPECT_LT((moonVelocity - moon.V).magnitude(), 1.0e-12);
}

TEST_F(DE440SourceTest, MjdPreservesTimesThatRoundToTheSameJd)
{
	if (!source->isAvailable())
		GTEST_SKIP();
	constexpr double firstMjd = 67898.5;
	const double secondMjd = firstMjd + 1.0e-10;
	volatile double firstJd = firstMjd + 2400000.5;
	volatile double secondJd = secondMjd + 2400000.5;
	ASSERT_EQ(firstJd, secondJd);
	Vec3D firstR[8], firstV[8], secondR[8], secondV[8], repeatedR[8], repeatedV[8];
	source->getAllPlanetStatesMjd(firstMjd, firstR, firstV);
	source->getAllPlanetStatesMjd(secondMjd, secondR, secondV);
	source->getAllPlanetStatesMjd(firstMjd, repeatedR, repeatedV);
	const Vec3D displacement = secondR[2] - firstR[2];
	EXPECT_GT(displacement.magnitude(), 1.0e-12);
	EXPECT_LT((displacement - firstV[2] * (secondMjd - firstMjd)).magnitude(), 1.0e-14);
	EXPECT_EQ((repeatedR[2] - firstR[2]).magnitude(), 0.0);
}

TEST_F(DE440SourceTest, MjdAndJdQueriesAgreeAtRepresentableEpochs)
{
	if (!source->isAvailable())
		GTEST_SKIP();
	for (double mjd : {0.0, 51544.5, 59169.77153276652, 67898.5, 90000.0}) {
		volatile double roundedJd = mjd + 2400000.5;
		const double representableMjd = roundedJd - 2400000.5;
		Vec3D r[8], v[8], mr[8], mv[8];
		source->getAllPlanetStates(roundedJd, r, v);
		source->getAllPlanetStatesMjd(representableMjd, mr, mv);
		for (int i = 0; i < 8; ++i) {
			EXPECT_EQ((r[i] - mr[i]).magnitude(), 0.0);
			EXPECT_EQ((v[i] - mv[i]).magnitude(), 0.0);
		}
	}
}

TEST_F(DE440SourceTest, LegacyIrel2TraceWithReconstructedSun)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	namespace ocf = orbital_core::forces;
	ocf::ForceContext ctx;
	ctx.setEpochMjd(2459167.6151007407 - 2400000.5);
	ctx.R = Vec3D(0.657898625443739671, 0.8019419912883772772, 0.00322060501030155136);
	ctx.V = Vec3D(-0.0124041212969580542, 0.01232340477254079694, 9.54401605106113476e-5);
	ctx.bodyMask = 0x00ffu | static_cast<std::uint16_t>(1u << 11);
	ctx.sunGm = source->gravitationalParameter(EphemerisPlanet::Sun);
	for (int index = 0; index < 8; ++index) {
		constexpr EphemerisPlanet bodies[8] = {
		    EphemerisPlanet::Mercury, EphemerisPlanet::Venus,  EphemerisPlanet::Earth,  EphemerisPlanet::Mars,
		    EphemerisPlanet::Jupiter, EphemerisPlanet::Saturn, EphemerisPlanet::Uranus, EphemerisPlanet::Neptune};
		ctx.planetGms[index] = source->gravitationalParameter(bodies[index]);
	}
	ctx.moonGm = source->gravitationalParameter(EphemerisPlanet::Moon);
	source->getAllPlanetStates(ctx.jd, ctx.planetPositions, ctx.planetVelocities, &ctx.moonPosition, &ctx.moonVelocity);
	ctx.havePlanets = true;
	ctx.havePlanetVelocities = true;
	ctx.ephemerisTerms.sunVelocity = reconstructedSolarVelocity(ctx);

	Vec3D acceleration;
	Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dadv = Eigen::Matrix3d::Zero();
	ocf::EihRelativity{}.evaluate(ctx, acceleration, dadr, dadv);
	acceleration += legacySolarJ2(ctx.R, ctx.sunGm);

	constexpr double interpolationTolerance = 2.0e-18;
	EXPECT_NEAR(acceleration.x, 4.12397702686420930e-12, interpolationTolerance);
	EXPECT_NEAR(acceleration.y, 6.68468642039907827e-12, interpolationTolerance);
	EXPECT_NEAR(acceleration.z, 3.07974727039830496e-14, interpolationTolerance);
}

TEST_F(DE440SourceTest, EarthEncounterTraceWithExplicitOriginConversion)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	namespace ocf = orbital_core::forces;
	ocf::ForceContext ctx;
	ctx.setEpochMjd(67899.95816614527);
	ctx.R = Vec3D(0.950243844223652379, 0.305009598984274433, -7.71352295634108384e-6);
	ctx.V = Vec3D(-4.85502783007580544e-3, 1.76649751565652394e-2, 9.97002368463838521e-5);
	ctx.bodyMask = 0x00ffu | static_cast<std::uint16_t>(1u << 11);
	ctx.sunGm = source->gravitationalParameter(EphemerisPlanet::Sun);
	for (int index = 0; index < 8; ++index) {
		constexpr EphemerisPlanet bodies[8] = {
		    EphemerisPlanet::Mercury, EphemerisPlanet::Venus,  EphemerisPlanet::Earth,  EphemerisPlanet::Mars,
		    EphemerisPlanet::Jupiter, EphemerisPlanet::Saturn, EphemerisPlanet::Uranus, EphemerisPlanet::Neptune};
		ctx.planetGms[index] = source->gravitationalParameter(bodies[index]);
	}
	ctx.moonGm = source->gravitationalParameter(EphemerisPlanet::Moon);
	ctx.earthGm = ctx.planetGms[2];
	ctx.ephemerisSource = source.get();
	ctx.computePartials = false;

	Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dadv = Eigen::Matrix3d::Zero();
	Vec3D planetary;
	ocf::PlanetaryPerturbations{}.evaluate(ctx, planetary, dadr, dadv);
	Vec3D reconstructedOrigin = ctx.moonPosition * (ctx.moonGm / std::pow(ctx.moonPosition.magnitude(), 3));
	for (int i = 0; i < 8; ++i)
		reconstructedOrigin +=
		    ctx.planetPositions[i] * (ctx.planetGms[i] / std::pow(ctx.planetPositions[i].magnitude(), 3));
	const Vec3D legacyPlanetary = planetary + ctx.ephemerisTerms.sunAcceleration - reconstructedOrigin;
	const Vec3D expectedMajorBodyAcceleration(-1.2460974135361438e-6, 3.862131541101663e-6, -9.173935632177603e-9);
	EXPECT_NEAR(legacyPlanetary.x, expectedMajorBodyAcceleration.x, 2.0e-16);
	EXPECT_NEAR(legacyPlanetary.y, expectedMajorBodyAcceleration.y, 2.0e-16);
	EXPECT_NEAR(legacyPlanetary.z, expectedMajorBodyAcceleration.z, 2.0e-16);
	const auto sunBefore = source->barycentricSunMjd(ctx.mjd - .002);
	const auto sunAfter = source->barycentricSunMjd(ctx.mjd + .002);
	ASSERT_TRUE(sunBefore && sunAfter);
	EXPECT_LT(((sunAfter->V - sunBefore->V) / .004 - ctx.ephemerisTerms.sunAcceleration).magnitude(),
	          ctx.ephemerisTerms.sunAcceleration.magnitude() * 1e-6);

	Vec3D solar;
	ctx.centralNewtonianGm = ctx.sunGm;
	ocf::SunGravity{}.evaluate(ctx, solar, dadr, dadv);
	const Vec3D expectedSolar(-2.82886824979624587e-4, -9.08011112826131024e-5, 2.29630955442759152e-9);
	EXPECT_LT((solar - expectedSolar).magnitude(), 2.0e-18);

	const Vec3D earthIndirect(-8.45338286149393528e-10, -2.852899559519991166e-10, 3.83183092251907736e-14);
	const Vec3D earthDirect(-1.23475873082755121e-6, 3.829384059215757794e-6, -9.51666168486391147e-9);
	EXPECT_NEAR(ctx.earthGm, 8.88769244670710224e-10, 1.0e-22);
	const Vec3D q = ctx.planetPositions[2] - ctx.R;
	const double qMag = q.magnitude();
	const double earthMag = ctx.planetPositions[2].magnitude();
	const Vec3D actualEarthDirect = q * (ctx.earthGm / (qMag * qMag * qMag));
	const Vec3D actualEarthIndirect = ctx.planetPositions[2] * (-ctx.earthGm / (earthMag * earthMag * earthMag));
	constexpr double pointMassTolerance = 2.0e-16;
	const Vec3D expectedEarthPosition(0.945682807492555066, 0.3191548412211076180, -4.28668224782643477e-5);
	EXPECT_NEAR(ctx.planetPositions[2].x, expectedEarthPosition.x, 5.0e-13);
	EXPECT_NEAR(ctx.planetPositions[2].y, expectedEarthPosition.y, 5.0e-13);
	EXPECT_NEAR(ctx.planetPositions[2].z, expectedEarthPosition.z, 5.0e-13);
	EXPECT_NEAR(actualEarthDirect.x, earthDirect.x, pointMassTolerance);
	EXPECT_NEAR(actualEarthDirect.y, earthDirect.y, pointMassTolerance);
	EXPECT_NEAR(actualEarthDirect.z, earthDirect.z, pointMassTolerance);
	EXPECT_LT((actualEarthDirect - earthDirect).magnitude(), pointMassTolerance);
	EXPECT_LT((actualEarthIndirect - earthIndirect).magnitude(), pointMassTolerance);

	Vec3D relativistic;
	ctx.ephemerisTerms.sunVelocity = reconstructedSolarVelocity(ctx);
	ocf::EihRelativity{}.evaluate(ctx, relativistic, dadr, dadv);
	relativistic += legacySolarJ2(ctx.R, ctx.sunGm);
	const Vec3D expectedRelativistic(7.96710732699091634e-12, 2.823323567840220308e-12, 3.39526183271393626e-15);
	EXPECT_LT((relativistic - expectedRelativistic).magnitude(), 2.0e-18);

	Vec3D earthJ2;
	ocf::EarthJ2{}.evaluate(ctx, earthJ2, dadr, dadv);
	const Vec3D expectedEarthJ2(-4.71733490587118050e-15, 3.12210868382546784e-14, 3.69696738421058433e-14);
	EXPECT_LT((earthJ2 - expectedEarthJ2).magnitude(), 3.0e-16);
}

TEST_F(DE440SourceTest, PositionReturnsEclipticCoordinates)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	Vec3D pos = source->getPlanetPosition(EphemerisPlanet::Earth, Constants::J2000_JD);
	double dist = pos.magnitude();
	double eclLat = std::asin(pos.z / dist);

	EXPECT_NEAR(eclLat, 0.0, 0.001) << "Earth ecliptic latitude should be ~0";
}

TEST_F(DE440SourceTest, InvalidPathReturnsUnavailable)
{
	DE440Source bad("/nonexistent/path.440");
	EXPECT_FALSE(bad.isAvailable());
}

TEST_F(DE440SourceTest, AllPlanetsQueryable)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	EphemerisPlanet planets[] = {EphemerisPlanet::Mercury, EphemerisPlanet::Venus,   EphemerisPlanet::Earth,
	                             EphemerisPlanet::Mars,    EphemerisPlanet::Jupiter, EphemerisPlanet::Saturn,
	                             EphemerisPlanet::Uranus,  EphemerisPlanet::Neptune};

	for (auto p : planets) {
		Vec3D pos = source->getPlanetPosition(p, 2460000.5);
		EXPECT_GT(pos.magnitude(), 0.1) << "Planet position should be nonzero";
	}
}

TEST_F(DE440SourceTest, ConcurrentQueriesProduceIdenticalResults)
{
	if (!source->isAvailable())
		GTEST_SKIP();

	std::vector<double> jds;
	jds.reserve(200);
	for (int i = 0; i < 200; ++i)
		jds.push_back(2460000.5 + i * 0.37);

	std::vector<Vec3D> baseline(jds.size());
	for (size_t i = 0; i < jds.size(); ++i)
		baseline[i] = source->getPlanetPosition(EphemerisPlanet::Jupiter, jds[i]);

	constexpr int nThreads = 8;
	std::vector<std::thread> workers;
	std::atomic<bool> ok{true};

	for (int t = 0; t < nThreads; ++t) {
		workers.emplace_back([&]() {
			for (size_t i = 0; i < jds.size(); ++i) {
				Vec3D p = source->getPlanetPosition(EphemerisPlanet::Jupiter, jds[i]);
				if (p.x != baseline[i].x || p.y != baseline[i].y || p.z != baseline[i].z) {
					ok = false;
					return;
				}
			}
		});
	}
	for (auto& w : workers)
		w.join();
	EXPECT_TRUE(ok) << "Per-thread ephem contexts must produce bit-identical results";
}
