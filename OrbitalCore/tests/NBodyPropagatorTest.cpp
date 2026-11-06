#include <gtest/gtest.h>
#include "NBodyPropagator.h"
#include "integrators/RA15Integrator.h"
#include "integrators/IAS15Integrator.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/hazard/ObsMonteCarlo.h"
#include "orbital_core/perturbations/PlanetaryPerturbations.h"
#include "orbital_core/perturbations/AsteroidPerturbations.h"
#include "hazard/CloseApproachKinematics.h"
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>

static_assert((PLANET_BODY_MASK & orbital_core::planetaryBodyMask(EphemerisPlanet::Pluto)) != 0);
static_assert((PLANET_BODY_MASK | MOON_BODY_MASK) == 0x18FF);

TEST(IntegratorCancellationTest, StopsInsideCorrectorStagesInBothDirectionsAndAllModes)
{
	class CancellingDynamics final : public IIntegrationDynamics {
	   public:
		Vec3D accelerationAtMjd(const Vec3D&, const Vec3D&, double, const NonGravModel&) const override
		{
			if (++evaluations == 3)
				stop.request_stop();
			return Vec3D(0, 0, 0);
		}

		void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D& v, double t, const NonGravModel& ng, Vec3D& a,
		                                Eigen::Matrix3d& dr, Eigen::Matrix3d& dv) const override
		{
			a = accelerationAtMjd(r, v, t, ng);
			dr.setZero();
			dv.setZero();
		}

		mutable std::stop_source stop;
		mutable int evaluations = 0;
	};

	RA15Integrator ra15;
	IAS15Integrator ias15;
	for (const IIntegrator* integrator : {static_cast<IIntegrator*>(&ra15), static_cast<IIntegrator*>(&ias15)}) {
		for (const double direction : {-1., 1.}) {
			for (int mode = 0; mode < 3; ++mode) {
				for (const bool preCancelled : {false, true}) {
					SCOPED_TRACE(integrator->name());
					SCOPED_TRACE(direction);
					SCOPED_TRACE(mode);
					SCOPED_TRACE(preCancelled);
					CancellingDynamics dynamics;
					if (preCancelled)
						dynamics.stop.request_stop();
					const IntegrationOptions options{.cancel = dynamics.stop.get_token()};
					const std::array<double, 1> targets{59000 + direction * 100};
					const Vec3D r(1, 0, 0), v(0, .017, 0);
					NonGravModel ng;
					IntegrationReport report;
					if (mode == 0) {
						std::array<IntegratedState, 1> output;
						report = integrator->integrateState(dynamics, ng, {r, v}, 59000, targets, output, options);
					} else if (mode == 1) {
						std::array<IntegratedStateWithStm, 1> output;
						report =
						    integrator->integrateStateWithStm(dynamics, ng, {r, v}, 59000, targets, output, options);
					} else {
						std::array<IntegratedStateFull, 1> output;
						report = integrator->integrateStateWithStmAndSens(dynamics, ng, {r, v}, 59000, targets, output,
						                                                  options);
					}
					EXPECT_TRUE(report.cancelled);
					EXPECT_FALSE(report.stoppedByObserver);
					EXPECT_EQ(report.completedOutputs, 0u);
					EXPECT_DOUBLE_EQ(report.reachedMjd, 59000);
					EXPECT_EQ(dynamics.evaluations, preCancelled ? 0 : 3);
					EXPECT_DOUBLE_EQ((report.finalState.R - r).magnitude(), 0);
				}
			}
		}
	}
}

class NBodyPropagatorTest : public ::testing::Test {
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
		de440 = std::make_unique<DE440Source>();
		if (de440->isAvailable())
			setDefaultEphemerisSource(de440.get());

		sb441 = std::make_unique<SpkReader>();
		if (sb441->isAvailable())
			setDefaultAsteroidSource(sb441.get());
	}

	void TearDown() override
	{
		setDefaultAsteroidSource(nullptr);
		setDefaultEphemerisSource(nullptr);
	}

	std::unique_ptr<DE440Source> de440;
	std::unique_ptr<SpkReader> sb441;
	NBodyPropagator propagator;
};

TEST_F(NBodyPropagatorTest, CopiesPreserveConfigurationAndPropagateIndependently)
{
	NBodyPropagator original(IntegratorKind::RA15);
	original.setStepSize(2.0);
	original.setPlanetMask(0);
	original.setUseSb441(false);
	original.setTargetNaifId(2000001);
	original.setInactiveInteriorMassPolicy(orbital_core::InactiveInteriorMassPolicy::SmoothCentralMass);
	original.nonGrav().setModel(NgModel::Yarkovsky);
	original.nonGrav().setEnabled(true);
	original.nonGrav().setParams(0.0, -1e-14, 0.0);
	auto configured = original.forEncounterMonitoring(*de440, 64, 0.25);

	NBodyPropagator copied = configured;
	NBodyPropagator assigned(IntegratorKind::IAS15);
	assigned = configured;
	assigned = assigned;
	for (const auto* p : {&copied, &assigned}) {
		EXPECT_EQ(p->integratorKind(), IntegratorKind::RA15);
		EXPECT_DOUBLE_EQ(p->stepSize(), 2.0);
		EXPECT_EQ(p->planetMask(), 0);
		EXPECT_FALSE(p->useSb441());
		EXPECT_EQ(p->targetNaifId(), 2000001);
		EXPECT_EQ(p->inactiveInteriorMassPolicy(), orbital_core::InactiveInteriorMassPolicy::SmoothCentralMass);
		EXPECT_TRUE(p->nonGravEnabled());
		EXPECT_EQ(p->nonGrav().model(), NgModel::Yarkovsky);
		EXPECT_DOUBLE_EQ(p->A2(), -1e-14);
		EXPECT_EQ(p->encounterStepPointCount(), 64);
		EXPECT_DOUBLE_EQ(p->earthEncounterRadiusAu(), 0.25);
	}
	configured.setIntegratorKind(IntegratorKind::IAS15);
	configured.setStepSize(5.0);
	configured.nonGrav().setParams(0.0, 0.0, 0.0);
	EXPECT_EQ(copied.integratorKind(), IntegratorKind::RA15);
	EXPECT_DOUBLE_EQ(assigned.stepSize(), 2.0);
	EXPECT_DOUBLE_EQ(copied.A2(), -1e-14);

	const Vec3D r(1.2, 0.2, 0.1), v(-0.001, 0.015, 0.0);
	constexpr double epoch = 2460000.5;
	const auto expected = copied.propagatePosition(r, v, epoch, epoch + 30.0);
	auto first = std::async(std::launch::async, [&] { return copied.propagatePosition(r, v, epoch, epoch + 30.0); });
	auto second = std::async(std::launch::async, [&] { return assigned.propagatePosition(r, v, epoch, epoch + 30.0); });
	EXPECT_LT((first.get() - expected).magnitude(), 1e-13);
	EXPECT_LT((second.get() - expected).magnitude(), 1e-13);
}

TEST_F(NBodyPropagatorTest, PlutoTargetMatchesExplicitPlutoDisabledIncludingPartials)
{
	if (!de440->isAvailable())
		GTEST_SKIP() << "DE ephemeris unavailable";
	auto automatic = propagator;
	automatic.setPlanetMask(0x18FF);
	automatic.setTargetNaifId(999);
	auto manual = propagator;
	manual.setPlanetMask(static_cast<std::uint16_t>(0x18FF & ~orbital_core::planetaryBodyMask(EphemerisPlanet::Pluto)));
	const Vec3D r(1.2, 0.2, 0.1), v(-0.001, 0.015, 0.0);
	constexpr double jd = 2460000.5;
	EXPECT_DOUBLE_EQ((automatic.accelerationAt(r, v, jd) - manual.accelerationAt(r, v, jd)).magnitude(), 0.0);
	EXPECT_TRUE(automatic.accelerationPartialR(r, v, jd).isApprox(manual.accelerationPartialR(r, v, jd), 0.0));
	EXPECT_TRUE(automatic.accelerationPartialV(r, v, jd).isApprox(manual.accelerationPartialV(r, v, jd), 0.0));
	EXPECT_EQ(NBodyPropagator(automatic).targetNaifId(), 999);
}

TEST_F(NBodyPropagatorTest, MoonTargetKeepsEarthAndMatchesManualMoonExclusion)
{
	if (!de440->isAvailable())
		GTEST_SKIP() << "DE ephemeris unavailable";
	constexpr double mjd = 60000.0, jd = mjd + 2400000.5;
	const auto earth = de440->queryMjd(EphemerisPlanet::Earth, mjd);
	const auto moon = de440->queryMjd(EphemerisPlanet::Moon, mjd);
	EXPECT_GT((moon.R - earth.R).magnitude(), 0.002);
	EXPECT_LT((moon.R - earth.R).magnitude(), 0.003);
	EXPECT_NEAR(de440->gravitationalParameter(EphemerisPlanet::Earth), Constants::GM_EARTH_ALONE, 1e-20);
	auto automatic = propagator;
	automatic.setPlanetMask(0x18FF);
	automatic.setTargetNaifId(301);
	auto manual = propagator;
	manual.setPlanetMask(static_cast<std::uint16_t>(0x18FF & ~MOON_BODY_MASK));
	ASSERT_TRUE(automatic.planetEnabled(EphemerisPlanet::Earth));
	ASSERT_FALSE(automatic.planetEnabled(EphemerisPlanet::Moon));
	EXPECT_DOUBLE_EQ(
	    (automatic.accelerationAt(moon.R, moon.V, jd) - manual.accelerationAt(moon.R, moon.V, jd)).magnitude(), 0.0);
	EXPECT_TRUE(automatic.accelerationPartialR(moon.R, moon.V, jd)
	                .isApprox(manual.accelerationPartialR(moon.R, moon.V, jd), 0.0));
	for (const double days : {1.0, 7.0, 30.0}) {
		const Vec3D actual = automatic.propagatePosition(moon.R, moon.V, jd, jd + days);
		const Vec3D expected = manual.propagatePosition(moon.R, moon.V, jd, jd + days);
		EXPECT_LT((actual - expected).magnitude(), 1e-13);
	}
}

TEST_F(NBodyPropagatorTest, Sb441CeresExclusionRemovesOnlyItsDirectAcceleration)
{
	if (!de440->isAvailable() || !sb441->isAvailable())
		GTEST_SKIP() << "DE or SB441 unavailable";
	constexpr double mjd = 60000.0;
	Vec3D ceres, velocity;
	ASSERT_TRUE(sb441->queryStateMjd(2000001, mjd, ceres, velocity));
	const Vec3D r = ceres + Vec3D(.05, .02, .01);
	const Vec3D v(.001, .01, .002);
	auto full = propagator;
	full.setUseSb441(true);
	auto excluded = full;
	excluded.setTargetNaifId(2000001);
	const Vec3D difference =
	    full.accelerationAt(r, v, mjd + 2400000.5) - excluded.accelerationAt(r, v, mjd + 2400000.5);
	const Vec3D q = ceres - r;
	const double kernelGm = de440->asteroidGravitationalParameter(2000001);
	const double sunGm = de440->gravitationalParameter(EphemerisPlanet::Sun);
	const double gm = kernelGm > 0 ? kernelGm : orbital_core::forces::detail::GM_CERES * (sunGm / Constants::GM_SUN);
	const Vec3D expected = q * (gm / std::pow(q.magnitude(), 3));
	EXPECT_LT((difference - expected).magnitude(), expected.magnitude() * 1e-7 + 1e-19);
	auto unrelated = full;
	unrelated.setTargetNaifId(2000999);
	EXPECT_DOUBLE_EQ(
	    (unrelated.accelerationAt(r, v, mjd + 2400000.5) - full.accelerationAt(r, v, mjd + 2400000.5)).magnitude(),
	    0.0);
}

TEST_F(NBodyPropagatorTest, Ra15NonGravControlPreservesStateStmAndSensitivityAccuracy)
{
	class SunNonGravDynamics final : public IIntegrationDynamics {
	   public:
		Vec3D accelerationAtMjd(const Vec3D& r, const Vec3D& v, double t, const NonGravModel& ng) const override
		{
			return r * (-Constants::GM_SUN / std::pow(r.magnitude(), 3)) + ng.acceleration(r, v, t);
		}

		void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D& v, double t, const NonGravModel& ng, Vec3D& a,
		                                Eigen::Matrix3d& dr, Eigen::Matrix3d& dv) const override
		{
			a = accelerationAtMjd(r, v, t, ng);
			ng.accelerationPartials(r, v, dr, dv);
			const Eigen::Vector3d x(r.x, r.y, r.z);
			dr += Constants::GM_SUN / std::pow(r.magnitude(), 3) *
			      (3 * x * x.transpose() / x.squaredNorm() - Eigen::Matrix3d::Identity());
		}
	} dynamics;

	const Vec3D r(1.2, .1, -.02), v(-.002, .015, .0003);
	constexpr double epoch = 59169.771532767;
	for (const auto kind : {NgModel::Yarkovsky, NgModel::InverseSquare, NgModel::MarsdenH2O}) {
		for (const bool zeroAmplitude : {false, true}) {
			SCOPED_TRACE(static_cast<int>(kind));
			SCOPED_TRACE(zeroAmplitude);
			NonGravModel ng;
			ng.setModel(kind);
			ng.setEnabled(true);
			ng.setParams(zeroAmplitude ? 0 : 1e-9, zeroAmplitude ? 0 : -2e-10, zeroAmplitude ? 0 : 3e-11);
			ng.setDT(60.0);
			RA15Integrator ra15;
			for (double direction : {-1.0, 1.0}) {
				SCOPED_TRACE(direction);
				const std::array<double, 2> targets{epoch + direction * 120, epoch + direction * 1500};
				std::array<IntegratedStateFull, 2> full, reference;
				std::array<IntegratedStateWithStm, 2> stm;
				const auto fullReport =
				    ra15.integrateStateWithStmAndSens(dynamics, ng, {r, v}, epoch, targets, full, {});
				const auto stmReport = ra15.integrateStateWithStm(dynamics, ng, {r, v}, epoch, targets, stm, {});
				IntegrationOptions fixed;
				fixed.initialStepDays = .5;
				fixed.nextEncounterStep = [direction](double, double, const IntegratedState&) -> std::optional<double> {
					return direction * .5;
				};
				const auto referenceReport =
				    ra15.integrateStateWithStmAndSens(dynamics, ng, {r, v}, epoch, targets, reference, fixed);
				EXPECT_EQ(referenceReport.unconvergedEncounterSteps, 0u);
				EXPECT_EQ(fullReport.completedOutputs, targets.size());
				EXPECT_LT(fullReport.acceptedSteps, 4 * stmReport.acceptedSteps);
				for (std::size_t i = 0; i < targets.size(); ++i) {
					SCOPED_TRACE(i);
					EXPECT_LT((full[i].R - stm[i].R).magnitude(), 1e-10);
					EXPECT_TRUE(full[i].Phi.isApprox(stm[i].Phi, 1e-9));
					EXPECT_LT((full[i].R - reference[i].R).magnitude(), 1e-10);
					EXPECT_LT((full[i].V - reference[i].V).magnitude(), 1e-11);
					EXPECT_TRUE(full[i].Phi.isApprox(reference[i].Phi, 1e-9));
					for (int col = 0; col < 4; ++col) {
						const double norm = reference[i].Sens.col(col).norm();
						EXPECT_LE((full[i].Sens.col(col) - reference[i].Sens.col(col)).norm(), 1e-9 * norm + 1e-15);
					}
				}
				if (zeroAmplitude) {
					constexpr double delta = 1e-9;
					for (int col = 0; col < 3; ++col) {
						SCOPED_TRACE(col);
						auto plus = ng, minus = ng;
						plus.setParams(col == 0 ? delta : 0, col == 1 ? delta : 0, col == 2 ? delta : 0);
						minus.setParams(-plus.A1(), -plus.A2(), -plus.A3());
						const std::array<double, 1> shortArc{targets.front()};
						std::array<IntegratedStateWithStm, 1> upper, lower;
						ra15.integrateStateWithStm(dynamics, plus, {r, v}, epoch, shortArc, upper, {});
						ra15.integrateStateWithStm(dynamics, minus, {r, v}, epoch, shortArc, lower, {});
						for (std::size_t i = 0; i < shortArc.size(); ++i) {
							const Vec3D dr = (upper[i].R - lower[i].R) / (2 * delta);
							const Vec3D dv = (upper[i].V - lower[i].V) / (2 * delta);
							const Eigen::Vector3d diffR(dr.x, dr.y, dr.z), diffV(dv.x, dv.y, dv.z);
							EXPECT_LE((full[i].Sens.col(col).head<3>() - diffR).norm(), 1e-6 * diffR.norm() + 1e-15);
							EXPECT_LE((full[i].Sens.col(col).tail<3>() - diffV).norm(), 1e-6 * diffV.norm() + 1e-15);
						}
					}
				}
				EXPECT_GT(full.front().Sens.leftCols<3>().norm(), 1.0);
				if ((kind == NgModel::MarsdenH2O || kind == NgModel::InverseSquare) && !zeroAmplitude) {
					constexpr double deltaDays = 0.01;
					auto plus = ng, minus = ng;
					plus.setDT(ng.DT() + deltaDays);
					minus.setDT(ng.DT() - deltaDays);
					const std::array<double, 1> shortArc{targets.front()};
					std::array<IntegratedStateWithStm, 1> upper, lower;
					ra15.integrateStateWithStm(dynamics, plus, {r, v}, epoch, shortArc, upper, fixed);
					ra15.integrateStateWithStm(dynamics, minus, {r, v}, epoch, shortArc, lower, fixed);
					const Vec3D dr = (upper.front().R - lower.front().R) / (2 * deltaDays);
					const Vec3D dv = (upper.front().V - lower.front().V) / (2 * deltaDays);
					const Eigen::Vector3d diffR(dr.x, dr.y, dr.z), diffV(dv.x, dv.y, dv.z);
					EXPECT_GT(full.front().Sens.col(3).norm(), 0.0);
					EXPECT_LE((full.front().Sens.col(3).head<3>() - diffR).norm(), 1e-3 * diffR.norm() + 1e-13);
					EXPECT_LE((full.front().Sens.col(3).tail<3>() - diffV).norm(), 1e-3 * diffV.norm() + 1e-14);
				} else
					EXPECT_TRUE(full.front().Sens.col(3).isZero());
			}
		}
	}
}

TEST_F(NBodyPropagatorTest, SurfaceStopMatches2020VvGrid)
{
	if (!de440->isAvailable() || !sb441->isAvailable())
		GTEST_SKIP();
	std::ifstream input(std::filesystem::path(__FILE__).parent_path() / "data/2020vv_surface_grid.txt");
	ASSERT_TRUE(input);
	std::string header;
	std::getline(input, header);
	std::ifstream hits(std::filesystem::path(__FILE__).parent_path() / "data/surface_stop_expected.txt");
	ASSERT_TRUE(hits);
	while (hits.peek() == '#')
		std::getline(hits, header);
	const double radius = orbital_core::hazard::planetRadiusKm(orbital_core::hazard::Planet::Earth);
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	propagator.setUseSb441(true);
	int index, legacyHit, count = 0;
	double epoch, x[6];
	while (input >> index >> legacyHit >> epoch) {
		for (double& value : x)
			input >> value;
		ASSERT_TRUE(input);
		int expectedIndex, expectedHit;
		ASSERT_TRUE(hits >> expectedIndex >> expectedHit);
		ASSERT_EQ(expectedIndex, index);
		SCOPED_TRACE(index);
		const auto terminal = propagator.propagateUntilSurface(
		    Vec3D(x[0], x[1], x[2]), Vec3D(x[3], x[4], x[5]), epoch, 2468000, EphemerisPlanet::Earth, radius, *de440,
		    propagator.nonGrav(), NBodyPropagator::SurfaceShape::Sphere);
		EXPECT_EQ(terminal.hitSurface, expectedHit != 0);
		if (terminal.hitSurface) {
			const auto earth = de440->queryJD(EphemerisPlanet::Earth, terminal.jd);
			const auto r = terminal.R - earth.R;
			EXPECT_NEAR(r.magnitude() * auKm, radius, .02);
			EXPECT_LT(r.dot(terminal.V - earth.V), 0);
		} else {
			EXPECT_DOUBLE_EQ(terminal.jd, 2468000);
		}
		++count;
	}
	EXPECT_EQ(count, 309);
}

TEST(SurfaceStopTest, DetectsInteriorCrossingInEitherTimeDirection)
{
	class FixedEarth final : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, Vec3D(1, 0, 0), Vec3D()}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "fixed-earth-test"; }
	} earth;

	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	const double radius = orbital_core::hazard::planetRadiusKm(orbital_core::hazard::Planet::Earth);
	constexpr double epoch = 2460000.5;
	for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
		for (double direction : {-1., 1.}) {
			for (double altitude : {-1., 1., 100000.}) {
				SCOPED_TRACE(altitude);
				NBodyPropagator propagator(kind);
				propagator.setAllPlanetsEnabled(false);
				propagator.setUseSb441(false);
				propagator.setStepSize(.002);
				const auto terminal = propagator.propagateUntilSurface(
				    Vec3D(1, -.0001, (radius + altitude) / auKm), Vec3D(0, direction * .1, 0), epoch,
				    epoch + direction * .002, EphemerisPlanet::Earth, radius, earth,
				    NBodyPropagator::SurfaceShape::Sphere);
				EXPECT_EQ(terminal.hitSurface, altitude < 0);
				if (terminal.hitSurface) {
					EXPECT_GT((terminal.jd - epoch) * direction, 0);
					EXPECT_LT((terminal.jd - epoch) * direction, .001);
					EXPECT_NEAR((terminal.R - Vec3D(1, 0, 0)).magnitude() * auKm, radius, .02);
				}
			}
		}
	}
}

TEST(ImpactMonteCarloTest, SurfaceCrossingOverridesConicAndKeepsEventWindow)
{
	namespace hz = orbital_core::hazard;

	class FixedMoon final : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, Vec3D(1, 0, 0), Vec3D()}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "fixed-moon-mc-test"; }

		std::optional<LunarLibration> lunarLibrationMjdTdb(double) override { return LunarLibration{{0, 0, 0}, 440}; }
	} moon;

	constexpr double epoch = 2460000.5;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	const double radius = hz::planetRadiusKm(hz::Planet::Moon);
	NBodyPropagator propagator(IntegratorKind::RA15);
	propagator.setAllPlanetsEnabled(false);
	propagator.setUseSb441(false);
	propagator.setStepSize(.002);
	const orbital_core::StateCov6x6 cov = orbital_core::StateCov6x6::Identity() * 1e-30;
	for (double altitude : {-1., 1.}) {
		SCOPED_TRACE(altitude);
		Orbit orbit;
		orbit.setEpoch(epoch);
		orbit.setR(Vec3D(1, -.0001, (radius + altitude) / auKm));
		orbit.setV(Vec3D(0, .02 / Constants::GAUSS_K, 0));
		hz::CloseApproachEvent event;
		event.planet = hz::Planet::Moon;
		event.jd = epoch + .0025;
		event.rPlanet = Vec3D(1, 0, 0);
		event.vInfKmS = .02 * auKm / 86400;
		const auto sample =
		    propagator.propagateAllWithVelocity(orbit.getR(), orbit.getV() * Constants::GAUSS_K, epoch, {event.jd})
		        .front();
		event.rTarget = sample.R;
		event.vTarget = sample.V;
		const Vec3D r = (sample.R - event.rPlanet) * auKm;
		const Vec3D v = sample.V * (auKm / 86400);
		const double mu = hz::planetMuKm3PerS2(event.planet);
		const Vec3D h = r.cross(v);
		const Vec3D eccentricity = v.cross(h) * (1. / mu) - r * (1. / r.magnitude());
		const double q = h.dot(h) / (mu * (1. + eccentricity.magnitude()));
		EXPECT_LT(q, radius);
		auto later = event;
		later.jd += 3;
		const std::array events{later, event};
		hz::ObsMcConfig cfg;
		cfg.numClones = 4;
		const auto state = hz::runStateCovarianceMonteCarlo(orbit, cov, propagator, events, cfg, &moon);
		NonGravModel ng;
		ng.setModel(NgModel::Yarkovsky);
		ng.setEnabled(true);
		orbital_core::FitParameterCovariance paramCov = orbital_core::FitParameterCovariance::Identity(7, 7) * 1e-30;
		const auto parameters =
		    hz::runParameterCovarianceMonteCarlo(orbit, paramCov, ng, propagator, events, cfg, &moon);
		for (const auto* result : {&state, &parameters}) {
			EXPECT_EQ(result->numConverged, cfg.numClones);
			EXPECT_EQ(result->numFailedFit, 0);
			EXPECT_EQ(result->hitsPerEvent[0], 0);
			EXPECT_EQ(result->hitsPerEvent[1], altitude < 0 ? cfg.numClones : 0);
			EXPECT_EQ(result->corridor.size(), result->numAnyImpactClones);
			for (const auto& hit : result->corridor) {
				EXPECT_EQ(hit.eventIndex, 1);
				EXPECT_GT(hit.tcaJd, event.jd);
				EXPECT_LT(hit.tcaJd, epoch + .006);
				EXPECT_TRUE(hit.bodyFixedValid);
				EXPECT_NEAR(hit.rBodyFixedKm.magnitude(), radius, .02);
			}
		}
	}
}

TEST(ImpactMonteCarloTest, MissingEphemeridesDoNotCountAsSuccessfulMisses)
{
	namespace hz = orbital_core::hazard;

	class Missing final : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double) override { throw std::runtime_error("missing"); }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return false; }

		std::string name() const override { return "missing-mc-test"; }
	} missing;

	Orbit orbit;
	orbit.setEpoch(2460000.5);
	orbit.setR(Vec3D(1, 0, 0));
	orbit.setV(Vec3D(0, 1, 0));
	hz::CloseApproachEvent event;
	event.jd = orbit.getEpoch() + 1;
	NBodyPropagator prop(IntegratorKind::RA15);
	prop.setAllPlanetsEnabled(false);
	prop.setUseSb441(false);
	hz::ObsMcConfig cfg;
	cfg.numClones = 4;
	const auto result = hz::runStateCovarianceMonteCarlo(orbit, orbital_core::StateCov6x6::Identity() * 1e-30, prop,
	                                                     {&event, 1}, cfg, &missing);
	EXPECT_EQ(result.numConverged, 0);
	EXPECT_EQ(result.numFailedFit, cfg.numClones);
	EXPECT_FALSE(result.isReliable());
	EXPECT_TRUE(result.corridor.empty());
}

TEST(NBodyInteriorStopTest, Ra15RefinesSphereWithStateStmAndSens)
{
	class FixedEarth final : public IEphemerisSource {
	   public:
		double lastMjd{};

		EphemerisRecord queryMjd(EphemerisPlanet, double mjd) override
		{
			lastMjd = mjd;
			return {mjd + 2400000.5, Vec3D(1, 0, 0), Vec3D()};
		}

		EphemerisRecord queryJD(EphemerisPlanet planet, double jd) override { return queryMjd(planet, jd - 2400000.5); }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "fixed-earth-interior-stop"; }
	} earth;

	constexpr double epoch = 2460000.5;
	constexpr double radiusKm = 6378.137;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	const std::array<double, 2> checkpoints{epoch + 0.000625, epoch + 0.00175};
	for (const bool withSensitivities : {false, true}) {
		SCOPED_TRACE(withSensitivities);
		NBodyPropagator propagator(IntegratorKind::RA15);
		propagator.setAllPlanetsEnabled(false);
		propagator.setUseSb441(false);
		propagator.setStepSize(0.0004);
		NonGravModel ng;
		if (withSensitivities) {
			ng.setModel(NgModel::InverseSquare);
			ng.setEnabled(true);
		}
		IntegrationReport report;
		const auto stopped =
		    propagator.propagateWithStmToSurface(Vec3D(1.0002, 0.00001, 0), Vec3D(-0.1, 0, 0), epoch, checkpoints,
		                                         EphemerisPlanet::Earth, radiusKm, earth, ng, &report);
		ASSERT_TRUE(stopped.hitSurface);
		EXPECT_TRUE(report.stoppedByObserver);
		EXPECT_DOUBLE_EQ(stopped.mjd, earth.lastMjd);
		EXPECT_DOUBLE_EQ(stopped.mjd, report.reachedMjd);
		EXPECT_DOUBLE_EQ(stopped.jd, stopped.mjd + 2400000.5);
		EXPECT_GT(stopped.jd, checkpoints.front());
		EXPECT_LE(stopped.jd, checkpoints.back());
		const double timeQuantum = std::nextafter(stopped.mjd, std::numeric_limits<double>::infinity()) - stopped.mjd;
		EXPECT_NEAR((stopped.state.R - Vec3D(1, 0, 0)).magnitude() * auKm, radiusKm,
		            stopped.state.V.magnitude() * timeQuantum * auKm + 2e-7);
		EXPECT_DOUBLE_EQ(stopped.state.R.x, report.finalState.R.x);
		EXPECT_DOUBLE_EQ(stopped.state.V.x, report.finalState.V.x);
		EXPECT_TRUE(stopped.state.Phi.allFinite());
		EXPECT_GT(stopped.state.Phi(0, 3), 0.0);
		if (withSensitivities)
			EXPECT_GT(stopped.state.SensA.norm(), 0.0);
		else
			EXPECT_TRUE(stopped.state.SensA.isZero());
	}
}

TEST(NBodyObservationTimeTest, SubJdTargetsRemainDistinctForStateStmAndSens)
{
	constexpr double epoch = 59139.25;
	const std::vector<double> targets{epoch + 1e-10, epoch, epoch - 1e-10};
	const Vec3D r(1, 0, 0), v(0, .0172, 0);
	for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
		SCOPED_TRACE(static_cast<int>(kind));
		NBodyPropagator p(kind);
		p.setAllPlanetsEnabled(false);
		p.setUseSb441(false);
		IntegrationReport stmReport, fullReport;
		const auto states = p.propagateAllWithVelocityMjd(r, v, epoch, targets, p.nonGrav());
		const auto stm = p.propagateAllWithStmMjd(r, v, epoch, targets, p.nonGrav(), &stmReport);
		const auto full = p.propagateAllWithStmAndSensMjd(r, v, epoch, targets, p.nonGrav(), &fullReport);
		for (std::size_t i = 0; i < targets.size(); ++i) {
			ASSERT_EQ(2400000.5 + targets[i], 2400000.5 + epoch);
			const double dt = targets[i] - epoch;
			EXPECT_NEAR(states[i].R.y, v.y * dt, 1e-24);
			EXPECT_NEAR(stm[i].R.y, v.y * dt, 1e-24);
			EXPECT_NEAR(full[i].R.y, v.y * dt, 1e-24);
			EXPECT_NEAR(stm[i].Phi(1, 4), dt, 1e-22);
			EXPECT_NEAR(full[i].Phi(1, 4), dt, 1e-22);
			EXPECT_LT((full[i].Phi - stm[i].Phi).norm(), 1e-20);
		}
		EXPECT_EQ(stmReport.completedOutputs, targets.size());
		EXPECT_EQ(fullReport.completedOutputs, targets.size());
		EXPECT_DOUBLE_EQ(stmReport.reachedMjd, epoch - 1e-10);
	}
}

TEST(NBodyObservationTimeTest, MjdAndJdAgreeAtRepresentableTargets)
{
	constexpr double epoch = 59139.25;
	const std::vector<double> mjds{epoch + .25, epoch - .25, epoch};
	std::vector<double> jds;
	for (double t : mjds)
		jds.push_back(t + 2400000.5);
	const Vec3D r(1, 0, 0), v(0, .0172, 0);
	for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
		NBodyPropagator p(kind);
		p.setAllPlanetsEnabled(false);
		p.nonGrav().setModel(NgModel::MarsdenH2O);
		p.nonGrav().setEnabled(true);
		p.nonGrav().setParams(1e-9, 2e-10, -1e-10);
		const auto a = p.propagateAllWithStmAndSensMjd(r, v, epoch, mjds, p.nonGrav());
		const auto b = p.propagateAllWithStmAndSens(r, v, epoch + 2400000.5, jds, p.nonGrav());
		for (std::size_t i = 0; i < mjds.size(); ++i) {
			EXPECT_LT((a[i].R - b[i].R).magnitude(), 1e-14);
			EXPECT_LT((a[i].V - b[i].V).magnitude(), 1e-14);
			EXPECT_LT((a[i].Phi - b[i].Phi).norm(), 1e-12);
			EXPECT_LT((a[i].SensA - b[i].SensA).norm(), 1e-12);
			if (mjds[i] != epoch)
				EXPECT_GT(a[i].SensA.norm(), .01);
		}
	}
}

TEST_F(NBodyPropagatorTest, ZeroPropagationReturnsInitial)
{
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, 0.0172, 0.0);
	double epoch = 2460000.5;

	Vec3D result = propagator.propagatePosition(R, V, epoch, epoch);

	EXPECT_NEAR(result.x, R.x, 1e-10);
	EXPECT_NEAR(result.y, R.y, 1e-10);
	EXPECT_NEAR(result.z, R.z, 1e-10);
}

TEST_F(NBodyPropagatorTest, CircularOrbitPeriod)
{
	setDefaultEphemerisSource(nullptr);
	double v_circ = std::sqrt(Constants::GAUSS_K * Constants::GAUSS_K);
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, v_circ, 0.0);
	double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(false);

	double period = 2.0 * Constants::PI / Constants::GAUSS_K;
	Vec3D result = propagator.propagatePosition(R, V, epoch, epoch + period);

	EXPECT_NEAR(result.x, R.x, 0.001) << "x after full period";
	EXPECT_NEAR(result.y, R.y, 0.001) << "y after full period";
	EXPECT_NEAR(result.z, R.z, 0.001) << "z after full period";
}

TEST_F(NBodyPropagatorTest, QuarterPeriodCircular)
{
	double v_circ = std::sqrt(Constants::GAUSS_K * Constants::GAUSS_K);
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, v_circ, 0.0);
	double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(false);

	double quarter = 0.5 * Constants::PI / Constants::GAUSS_K;
	Vec3D result = propagator.propagatePosition(R, V, epoch, epoch + quarter);

	EXPECT_NEAR(result.x, 0.0, 0.01);
	EXPECT_NEAR(result.y, 1.0, 0.01);
}

TEST_F(NBodyPropagatorTest, PropagateAllBatch)
{
	double v_circ = std::sqrt(Constants::GAUSS_K * Constants::GAUSS_K);
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, v_circ, 0.0);
	double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(false);

	std::vector<double> jds = {epoch + 10.0, epoch + 20.0, epoch + 30.0};
	auto results = propagator.propagateAll(R, V, epoch, jds);

	ASSERT_EQ(results.size(), 3u);
	for (const auto& pos : results) {
		double dist = pos.magnitude();
		EXPECT_NEAR(dist, 1.0, 0.005) << "Circular orbit stays at 1 AU";
	}
}

TEST_F(NBodyPropagatorTest, PropagateAllMixedDirections)
{
	double v_circ = std::sqrt(Constants::GAUSS_K * Constants::GAUSS_K);
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, v_circ, 0.0);
	double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(false);

	std::vector<double> jds = {epoch - 5.0, epoch + 10.0, epoch - 10.0, epoch + 5.0};
	auto results = propagator.propagateAll(R, V, epoch, jds);

	ASSERT_EQ(results.size(), 4u);
	for (const auto& pos : results) {
		double dist = pos.magnitude();
		EXPECT_NEAR(dist, 1.0, 0.005);
	}
}

TEST_F(NBodyPropagatorTest, Ra15BatchMatchesIndependentTargets)
{
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K, 0.0);
	constexpr double epoch = 2460000.5;
	const std::vector<double> targets{epoch + 25.0, epoch - 7.0, epoch, epoch + 2.0, epoch - 19.0, epoch + 25.0};

	propagator.setIntegratorKind(IntegratorKind::RA15);
	propagator.setAllPlanetsEnabled(false);
	const auto batch = propagator.propagateAllWithVelocity(R, V, epoch, targets);

	ASSERT_EQ(batch.size(), targets.size());
	for (std::size_t index = 0; index < targets.size(); ++index) {
		NBodyPropagator independent(IntegratorKind::RA15);
		independent.setAllPlanetsEnabled(false);
		const auto single = independent.propagateAllWithVelocity(R, V, epoch, {targets[index]});
		ASSERT_EQ(single.size(), 1u);
		EXPECT_LT((batch[index].R - single[0].R).magnitude(), 2.0e-11) << "position at target " << index;
		EXPECT_LT((batch[index].V - single[0].V).magnitude(), 2.0e-12) << "velocity at target " << index;
	}
}

TEST_F(NBodyPropagatorTest, Ra15BatchStmMatchesIndependentTargets)
{
	const Vec3D R(1.2, 0.1, -0.02);
	const Vec3D V(-0.002, 0.015, 0.0003);
	constexpr double epoch = 2460000.5;
	const std::vector<double> targets{epoch + 3.0, epoch + 12.0};

	propagator.setIntegratorKind(IntegratorKind::RA15);
	propagator.setAllPlanetsEnabled(false);
	const auto batch = propagator.propagateAllWithStm(R, V, epoch, targets);

	ASSERT_EQ(batch.size(), targets.size());
	for (std::size_t index = 0; index < targets.size(); ++index) {
		NBodyPropagator independent(IntegratorKind::RA15);
		independent.setAllPlanetsEnabled(false);
		const auto single = independent.propagateAllWithStm(R, V, epoch, {targets[index]});
		ASSERT_EQ(single.size(), 1u);
		EXPECT_LT((batch[index].R - single[0].R).magnitude(), 2.0e-11);
		EXPECT_LT((batch[index].V - single[0].V).magnitude(), 2.0e-12);
		EXPECT_LT((batch[index].Phi - single[0].Phi).norm(), 2.0e-9);
	}
}

TEST_F(NBodyPropagatorTest, Ra15ChainedStmMatchesContinuousStm)
{
	const Vec3D R(1.2, 0.1, -0.02);
	const Vec3D V(-0.002, 0.015, 0.0003);
	constexpr double epoch = 2460000.5;
	const std::vector<double> targets{epoch + 3.0, epoch + 12.0, epoch + 40.0};

	propagator.setIntegratorKind(IntegratorKind::RA15);
	propagator.setAllPlanetsEnabled(false);

	IntegrationReport continuousReport;
	IntegrationReport chainedReport;
	const auto continuous =
	    propagator.propagateAllWithStm(R, V, epoch, targets, propagator.nonGrav(), &continuousReport);
	const auto chained =
	    propagator.propagateAllWithStmChained(R, V, epoch, targets, propagator.nonGrav(), &chainedReport);

	ASSERT_EQ(chained.size(), continuous.size());
	EXPECT_EQ(chainedReport.completedOutputs, targets.size());
	for (std::size_t index = 0; index < targets.size(); ++index) {
		EXPECT_LT((chained[index].R - continuous[index].R).magnitude(), 2.0e-11);
		EXPECT_LT((chained[index].V - continuous[index].V).magnitude(), 2.0e-12);
		EXPECT_LT((chained[index].Phi - continuous[index].Phi).norm(), 2.0e-9);
	}
}

TEST_F(NBodyPropagatorTest, InteriorPericenterMapPreservesTwoBodyInvariants)
{
	constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
	constexpr double kAuDayKmS = kAuKm / Constants::SECONDS_PER_DAY;
	const double radiusAu = orbital_core::hazard::planetRadiusKm(orbital_core::hazard::Planet::Earth) / kAuKm;
	const Vec3D r0(radiusAu, 0.0, 0.0);
	const Vec3D v0(-5.0 / kAuDayKmS, 14.0 / kAuDayKmS, 1.0 / kAuDayKmS);

	const auto mapped =
	    orbital_core::hazard::mapInteriorEncounterToPericenter(r0, v0, orbital_core::hazard::Planet::Earth);
	ASSERT_TRUE(mapped.valid);
	EXPECT_GT(mapped.deltaDays, 0.0);
	EXPECT_LT(mapped.relativePosition.magnitude(), radiusAu);
	EXPECT_NEAR(mapped.relativePosition.dot(mapped.relativeVelocity), 0.0, 1.0e-13);

	const double mu = orbital_core::hazard::planetMuKm3PerS2(orbital_core::hazard::Planet::Earth) *
	                  Constants::SECONDS_PER_DAY * Constants::SECONDS_PER_DAY / (kAuKm * kAuKm * kAuKm);
	const auto energy = [mu](const Vec3D& r, const Vec3D& v) { return 0.5 * v.dot(v) - mu / r.magnitude(); };
	const auto angularMomentum2 = [](const Vec3D& r, const Vec3D& v) { return r.cross(v).dot(r.cross(v)); };
	EXPECT_NEAR(energy(mapped.relativePosition, mapped.relativeVelocity), energy(r0, v0), 1.0e-13);
	EXPECT_NEAR(angularMomentum2(mapped.relativePosition, mapped.relativeVelocity), angularMomentum2(r0, v0), 1.0e-17);

	const auto propagateFixed = [mu, dt = mapped.deltaDays](Vec3D r, Vec3D v) {
		constexpr int kSteps = 4096;
		const double h = dt / kSteps;
		const auto acceleration = [mu](const Vec3D& position) {
			const double radius = position.magnitude();
			return position * (-mu / (radius * radius * radius));
		};
		for (int stepIndex = 0; stepIndex < kSteps; ++stepIndex) {
			const Vec3D k1r = v;
			const Vec3D k1v = acceleration(r);
			const Vec3D k2r = v + k1v * (0.5 * h);
			const Vec3D k2v = acceleration(r + k1r * (0.5 * h));
			const Vec3D k3r = v + k2v * (0.5 * h);
			const Vec3D k3v = acceleration(r + k2r * (0.5 * h));
			const Vec3D k4r = v + k3v * h;
			const Vec3D k4v = acceleration(r + k3r * h);
			r = r + (k1r + k2r * 2.0 + k3r * 2.0 + k4r) * (h / 6.0);
			v = v + (k1v + k2v * 2.0 + k3v * 2.0 + k4v) * (h / 6.0);
		}
		return std::pair<Vec3D, Vec3D>{r, v};
	};
	const auto fixedNominal = propagateFixed(r0, v0);
	EXPECT_LT((fixedNominal.first - mapped.relativePosition).magnitude(), 1.0e-13);
	EXPECT_LT((fixedNominal.second - mapped.relativeVelocity).magnitude(), 1.0e-12);

	Eigen::Matrix<double, 6, 6> finiteDifference;
	for (int column = 0; column < 6; ++column) {
		const double step = column < 3 ? 1.0e-10 : 1.0e-9;
		Vec3D rPlus = r0;
		Vec3D rMinus = r0;
		Vec3D vPlus = v0;
		Vec3D vMinus = v0;
		auto perturb = [column, step](Vec3D& r, Vec3D& v, double sign) {
			double* component = nullptr;
			switch (column) {
				case 0:
					component = &r.x;
					break;
				case 1:
					component = &r.y;
					break;
				case 2:
					component = &r.z;
					break;
				case 3:
					component = &v.x;
					break;
				case 4:
					component = &v.y;
					break;
				default:
					component = &v.z;
					break;
			}
			*component += sign * step;
		};
		perturb(rPlus, vPlus, 1.0);
		perturb(rMinus, vMinus, -1.0);
		const auto mappedPlus = propagateFixed(rPlus, vPlus);
		const auto mappedMinus = propagateFixed(rMinus, vMinus);
		const double inverseSpan = 0.5 / step;
		const Vec3D dr = (mappedPlus.first - mappedMinus.first) * inverseSpan;
		const Vec3D dv = (mappedPlus.second - mappedMinus.second) * inverseSpan;
		finiteDifference.col(column) << dr.x, dr.y, dr.z, dv.x, dv.y, dv.z;
	}
	const double stmScale = std::max(1.0, finiteDifference.norm());
	EXPECT_LT((mapped.stateTransition - finiteDifference).norm() / stmScale, 2.0e-5)
	    << "analytic:\n"
	    << mapped.stateTransition << "\nfinite difference:\n"
	    << finiteDifference;
}

TEST(InteriorPericenterTest, EventJacobianMatchesVariableTimeFiniteDifferences)
{
	namespace hz = orbital_core::hazard;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	constexpr double auDayKmS = auKm / Constants::SECONDS_PER_DAY;
	Eigen::Matrix<double, 6, 1> initial;
	initial << hz::planetRadiusKm(hz::Planet::Earth) / auKm, 0.0, 0.0, -5.0 / auDayKmS, 14.0 / auDayKmS, 1.0 / auDayKmS;
	const auto map = [](const Eigen::Matrix<double, 6, 1>& state) {
		return hz::mapInteriorEncounterToPericenter(Vec3D(state(0), state(1), state(2)),
		                                            Vec3D(state(3), state(4), state(5)), hz::Planet::Earth);
	};
	const auto center = map(initial);
	const auto derivative = center.eventJacobian(hz::Planet::Earth);
	ASSERT_TRUE(derivative);
	Eigen::Matrix<double, 6, 6> finiteDifference;
	for (int column = 0; column < 6; ++column) {
		const double step = column < 3 ? 1.0e-10 : 1.0e-9;
		auto plus = initial, minus = initial;
		plus(column) += step;
		minus(column) -= step;
		const auto a = map(plus), b = map(minus);
		ASSERT_TRUE(a.valid && b.valid);
		const Vec3D dr = (a.relativePosition - b.relativePosition) / (2.0 * step);
		const Vec3D dv = (a.relativeVelocity - b.relativeVelocity) / (2.0 * step);
		finiteDifference.col(column) << dr.x, dr.y, dr.z, dv.x, dv.y, dv.z;
		EXPECT_LT((derivative->col(column) - finiteDifference.col(column)).norm() / finiteDifference.col(column).norm(),
		          1.0e-6)
		    << column;
	}
	EXPECT_GT((center.stateTransition - finiteDifference).norm() / finiteDifference.norm(), 0.01);
	Eigen::Matrix<double, 1, 6> normal;
	normal << center.relativeVelocity.x, center.relativeVelocity.y, center.relativeVelocity.z,
	    center.relativePosition.x, center.relativePosition.y, center.relativePosition.z;
	EXPECT_LT((normal * *derivative).norm(), 1.0e-15);
	EXPECT_FALSE(hz::InteriorPericenterMapping{}.eventJacobian(hz::Planet::Earth));
}

TEST_F(NBodyPropagatorTest, PlanetDataGravitationalParametersMatchDe440Header)
{
	if (!de440 || !de440->isAvailable())
		GTEST_SKIP();

	const double sunGm = de440->gravitationalParameter(EphemerisPlanet::Sun);
	EXPECT_NEAR(Constants::GM_SUN, sunGm, sunGm * 2.0e-15);
	for (std::size_t index = 0; index < std::size(orbital_core::forces::detail::kPerturbers); ++index) {
		const auto& body = orbital_core::forces::detail::kPerturbers[index];
		const double expected = de440->gravitationalParameter(body.planet);
		EXPECT_NEAR(body.gm, expected, expected * 2.0e-15);
	}
	const double moonGm = de440->gravitationalParameter(EphemerisPlanet::Moon);
	EXPECT_NEAR(orbital_core::forces::detail::GM_MOON, moonGm, moonGm * 2.0e-15);

	constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
	constexpr double kAu3Day2ToKm3S2 =
	    (kAuKm * kAuKm * kAuKm) / (Constants::SECONDS_PER_DAY * Constants::SECONDS_PER_DAY);
	for (const auto planet : orbital_core::hazard::kAllPlanets) {
		const double expected =
		    de440->gravitationalParameter(orbital_core::hazard::toEphemerisPlanet(planet)) * kAu3Day2ToKm3S2;
		ASSERT_GT(expected, 0.0);
		EXPECT_NEAR(orbital_core::hazard::planetMuKm3PerS2(planet), expected, expected * 2.0e-15);
	}
}

TEST_F(NBodyPropagatorTest, IntegratorSelectionIsExplicit)
{
	EXPECT_EQ(propagator.integratorKind(), IntegratorKind::RA15);
	propagator.setIntegratorKind(IntegratorKind::IAS15);
	EXPECT_EQ(propagator.integratorKind(), IntegratorKind::IAS15);
	EXPECT_STREQ(propagator.integratorName(), "IAS15");
	EXPECT_DOUBLE_EQ(propagator.fixedAccuracyParameter(), 1.0e-9);

	propagator.setIntegratorKind(IntegratorKind::RA15);
	EXPECT_EQ(propagator.integratorKind(), IntegratorKind::RA15);
	EXPECT_STREQ(propagator.integratorName(), "RA15");
	EXPECT_DOUBLE_EQ(propagator.fixedAccuracyParameter(), 1.0e-12);
}

TEST_F(NBodyPropagatorTest, UnsupportedIntegratorIsRejected)
{
	EXPECT_THROW(NBodyPropagator(static_cast<IntegratorKind>(0)), std::invalid_argument);
}

TEST_F(NBodyPropagatorTest, EncounterClonePreservesIntegratorAndInitialStep)
{
	for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
		propagator.setIntegratorKind(kind);
		propagator.setStepSize(2.0);
		const auto clone = propagator.forEncounterMonitoring(*de440, 600, 0.2);
		EXPECT_EQ(clone.integratorKind(), kind);
		EXPECT_DOUBLE_EQ(clone.stepSize(), 2.0);
		EXPECT_DOUBLE_EQ(clone.fixedAccuracyParameter(), propagator.fixedAccuracyParameter());
	}
}

TEST_F(NBodyPropagatorTest, PlanetMaskDefault)
{
	EXPECT_EQ(propagator.planetMask(), FULL_BODY_MASK);
	EXPECT_TRUE(propagator.planetEnabled(EphemerisPlanet::Jupiter));
	EXPECT_TRUE(propagator.planetEnabled(EphemerisPlanet::Mercury));
	EXPECT_TRUE(propagator.planetEnabled(EphemerisPlanet::Moon));
}

TEST_F(NBodyPropagatorTest, DisablePlanet)
{
	propagator.setPlanetEnabled(EphemerisPlanet::Jupiter, false);
	EXPECT_FALSE(propagator.planetEnabled(EphemerisPlanet::Jupiter));
	EXPECT_TRUE(propagator.planetEnabled(EphemerisPlanet::Saturn));
	EXPECT_EQ(propagator.planetMask(), static_cast<uint16_t>(FULL_BODY_MASK & ~(1u << 4)));
}

TEST_F(NBodyPropagatorTest, SetAllPlanets)
{
	propagator.setAllPlanetsEnabled(false);
	EXPECT_EQ(propagator.planetMask(), 0x00);

	propagator.setAllPlanetsEnabled(true);
	EXPECT_EQ(propagator.planetMask(), FULL_BODY_MASK);
}

TEST_F(NBodyPropagatorTest, PerturbationsAffectResult)
{
	if (!de440 || !de440->isAvailable())
		GTEST_SKIP();

	Vec3D R(2.5, 0.0, 0.0);
	double v_circ = Constants::GAUSS_K / std::sqrt(2.5);
	Vec3D V(0.0, v_circ, 0.0);
	double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(true);
	Vec3D withPert = propagator.propagatePosition(R, V, epoch, epoch + 100.0);

	propagator.setAllPlanetsEnabled(false);
	Vec3D noPert = propagator.propagatePosition(R, V, epoch, epoch + 100.0);

	Vec3D diff = withPert - noPert;
	double diffMag = diff.magnitude();

	EXPECT_GT(diffMag, 1e-6) << "Perturbations should have measurable effect";
	EXPECT_LT(diffMag, 0.01) << "Perturbations shouldn't be huge over 100 days";
}

TEST_F(NBodyPropagatorTest, JupiterDominatesPerturbation)
{
	if (!de440 || !de440->isAvailable())
		GTEST_SKIP();

	Vec3D R(2.5, 0.0, 0.0);
	double v_circ = Constants::GAUSS_K / std::sqrt(2.5);
	Vec3D V(0.0, v_circ, 0.0);
	double epoch = 2460000.5;

	propagator.setPlanetMask(1 << 4);
	Vec3D jupOnly = propagator.propagatePosition(R, V, epoch, epoch + 100.0);

	propagator.setAllPlanetsEnabled(true);
	Vec3D allPlanets = propagator.propagatePosition(R, V, epoch, epoch + 100.0);

	propagator.setAllPlanetsEnabled(false);
	Vec3D noPert = propagator.propagatePosition(R, V, epoch, epoch + 100.0);

	double jupEffect = (jupOnly - noPert).magnitude();
	double totalEffect = (allPlanets - noPert).magnitude();

	EXPECT_GT(jupEffect, totalEffect * 0.3);
}

TEST_F(NBodyPropagatorTest, NgOverrideMatchesSharedState)
{
	if (!de440 || !de440->isAvailable())
		GTEST_SKIP();

	Vec3D R(2.5, 0.0, 0.0);
	double v_circ = Constants::GAUSS_K / std::sqrt(2.5);
	Vec3D V(0.0, v_circ, 0.0);
	double epoch = 2460000.5;

	std::vector<double> jds;
	for (int i = 0; i < 20; ++i)
		jds.push_back(epoch + i * 5.0);

	propagator.nonGrav().setModel(NgModel::MarsdenH2O);
	propagator.nonGrav().setEnabled(true);
	propagator.nonGrav().setParams(1e-9, 5e-10, -3e-10);
	propagator.nonGrav().setDT(2.0);
	auto shared = propagator.propagateAll(R, V, epoch, jds);

	NonGravModel override = propagator.nonGrav();
	propagator.nonGrav().setParams(0.0, 0.0, 0.0);
	propagator.nonGrav().setDT(0.0);
	propagator.nonGrav().setEnabled(false);
	auto viaOverride = propagator.propagateAll(R, V, epoch, jds, override);

	ASSERT_EQ(shared.size(), viaOverride.size());
	for (size_t i = 0; i < shared.size(); ++i) {
		EXPECT_DOUBLE_EQ(shared[i].x, viaOverride[i].x) << "x mismatch at jd[" << i << "]";
		EXPECT_DOUBLE_EQ(shared[i].y, viaOverride[i].y) << "y mismatch at jd[" << i << "]";
		EXPECT_DOUBLE_EQ(shared[i].z, viaOverride[i].z) << "z mismatch at jd[" << i << "]";
	}

	EXPECT_DOUBLE_EQ(propagator.nonGrav().A1(), 0.0);
	EXPECT_DOUBLE_EQ(propagator.nonGrav().A2(), 0.0);
	EXPECT_DOUBLE_EQ(propagator.nonGrav().A3(), 0.0);
	EXPECT_FALSE(propagator.nonGrav().enabled());
}

TEST_F(NBodyPropagatorTest, StmMatchesPropagateFD)
{
	if (!de440 || !de440->isAvailable())
		GTEST_SKIP();

	propagator.nonGrav().setEnabled(false);
	propagator.setAllPlanetsEnabled(true);

	Vec3D R0(2.3, 0.7, -0.4);
	Vec3D V0(-0.006, 0.012, 0.003);
	double epoch = 2460000.5;
	double target = epoch + 30.0;

	auto stm = propagator.propagateAllWithStm(R0, V0, epoch, {target});
	ASSERT_EQ(stm.size(), 1u);
	const auto& Phi = stm[0].Phi;

	const double h = 1e-6;
	auto propR = [&](const Vec3D& R, const Vec3D& V) { return propagator.propagatePosition(R, V, epoch, target); };
	Vec3D R_base = propR(R0, V0);

	auto setAxis = [](Vec3D& v, int dim, double val) { (dim == 0 ? v.x : dim == 1 ? v.y : v.z) = val; };

	for (int j = 0; j < 3; ++j) {
		Vec3D Rp = R0, Rm = R0;
		setAxis(Rp, j, (j == 0 ? R0.x : j == 1 ? R0.y : R0.z) + h);
		setAxis(Rm, j, (j == 0 ? R0.x : j == 1 ? R0.y : R0.z) - h);
		Vec3D Rpf = propR(Rp, V0);
		Vec3D Rmf = propR(Rm, V0);
		double fd0 = (Rpf.x - Rmf.x) / (2.0 * h);
		double fd1 = (Rpf.y - Rmf.y) / (2.0 * h);
		double fd2 = (Rpf.z - Rmf.z) / (2.0 * h);
		EXPECT_NEAR(Phi(0, j), fd0, 1e-5 * std::max(1.0, std::abs(fd0))) << "Phi(0," << j << ")";
		EXPECT_NEAR(Phi(1, j), fd1, 1e-5 * std::max(1.0, std::abs(fd1))) << "Phi(1," << j << ")";
		EXPECT_NEAR(Phi(2, j), fd2, 1e-5 * std::max(1.0, std::abs(fd2))) << "Phi(2," << j << ")";
	}

	for (int j = 0; j < 3; ++j) {
		Vec3D Vp = V0, Vm = V0;
		setAxis(Vp, j, (j == 0 ? V0.x : j == 1 ? V0.y : V0.z) + h);
		setAxis(Vm, j, (j == 0 ? V0.x : j == 1 ? V0.y : V0.z) - h);
		Vec3D Rpf = propR(R0, Vp);
		Vec3D Rmf = propR(R0, Vm);
		double fd0 = (Rpf.x - Rmf.x) / (2.0 * h);
		double fd1 = (Rpf.y - Rmf.y) / (2.0 * h);
		double fd2 = (Rpf.z - Rmf.z) / (2.0 * h);
		const int col = 3 + j;
		EXPECT_NEAR(Phi(0, col), fd0, 1e-5 * std::max(1.0, std::abs(fd0))) << "Phi(0," << col << ")";
		EXPECT_NEAR(Phi(1, col), fd1, 1e-5 * std::max(1.0, std::abs(fd1))) << "Phi(1," << col << ")";
		EXPECT_NEAR(Phi(2, col), fd2, 1e-5 * std::max(1.0, std::abs(fd2))) << "Phi(2," << col << ")";
	}

	EXPECT_NEAR(stm[0].R.x, R_base.x, 1e-10);
	EXPECT_NEAR(stm[0].R.y, R_base.y, 1e-10);
	EXPECT_NEAR(stm[0].R.z, R_base.z, 1e-10);
}

TEST_F(NBodyPropagatorTest, RelativisticVelocityPartialsMatchFullAcceleration)
{
	if (!de440 || !de440->isAvailable())
		GTEST_SKIP();

	Vec3D r(2.3, 0.7, -0.4);
	Vec3D v(-0.006, 0.012, 0.003);
	double jd = 2460000.5;

	propagator.nonGrav().setEnabled(false);
	propagator.setAllPlanetsEnabled(true);

	const auto analytic = propagator.accelerationPartialV(r, v, jd);
	EXPECT_GT(analytic.norm(), 0.0);
	Eigen::Matrix3d numerical;
	constexpr double h = 1e-3;
	for (int column = 0; column < 3; ++column) {
		auto plus = v, minus = v;
		double* positive = column == 0 ? &plus.x : column == 1 ? &plus.y : &plus.z;
		double* negative = column == 0 ? &minus.x : column == 1 ? &minus.y : &minus.z;
		*positive += h;
		*negative -= h;
		const Vec3D derivative =
		    (propagator.accelerationAt(r, plus, jd) - propagator.accelerationAt(r, minus, jd)) / (2 * h);
		numerical.col(column) << derivative.x, derivative.y, derivative.z;
	}
	EXPECT_LT((analytic - numerical).norm(), analytic.norm() * 2e-6);
}
