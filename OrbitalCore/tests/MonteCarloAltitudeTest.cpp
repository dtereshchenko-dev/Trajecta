#include <gtest/gtest.h>
#include "NBodyPropagator.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include <cmath>
#include <limits>

namespace {
	namespace hz = orbital_core::hazard;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	constexpr double epoch = 2460000.5;
	constexpr double mjd = 60000.0;

	class FixedEarth : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, {1, 0, 0}, {}}; }
		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }
		bool isAvailable() const override { return true; }
		std::string name() const override { return "altitude-test"; }
	};

	class MovingEarth final : public FixedEarth {
	   public:
		explicit MovingEarth(const Vec3D& velocity) : velocity_(velocity) {}
		EphemerisRecord queryJD(EphemerisPlanet planet, double jd) override
		{
			return queryMjd(planet, orbital_core::time::mjdFromJd(jd));
		}
		EphemerisRecord queryMjd(EphemerisPlanet, double t) override
		{
			return {orbital_core::time::jdFromMjd(t), Vec3D(1, 0, 0) + velocity_ * (t - mjd), velocity_};
		}
	   private:
		Vec3D velocity_;
	};

	Vec3D fromFixed(const Vec3D& r)
	{
		return {hz::icrfToBodyFixed(hz::Planet::Earth, {1, 0, 0}, epoch).dot(r),
		        hz::icrfToBodyFixed(hz::Planet::Earth, {0, 1, 0}, epoch).dot(r),
		        hz::icrfToBodyFixed(hz::Planet::Earth, {0, 0, 1}, epoch).dot(r)};
	}

	double height(const Vec3D& r)
	{
		return hz::referenceSurfaceSignedAltitudeKm(r, Constants::EARTH_EQ_RADIUS_KM, Constants::EARTH_FLATTENING);
	}

	class MonteCarloAltitudeTest : public ::testing::Test {
	   protected:
		FixedEarth eph;
		IEphemerisSource* previous{};
		orbital_core::OrbitSolution source;
		orbital_core::NBodySettings forces;
		hz::IntervalMcConfig cfg;

		void SetUp() override
		{
			previous = getDefaultEphemerisSource();
			setDefaultEphemerisSource(&eph);
			source.nominal.setEpoch(epoch);
			const Vec3D r = fromFixed({4600, 0, 4600});
			source.nominal.setR(Vec3D(1, 0, 0) + r * (1.0 / auKm));
			source.nominal.setV(r * (-10.0 / r.magnitude() * 86400.0 / auKm / Constants::GAUSS_K));
			source.covariance = orbital_core::StateCov6x6::Identity() * 1e-22;
			forces.enabled = true;
			forces.planetMask = 0;
			forces.stepSize = .002;
			cfg.startMjdTdb = mjd;
			cfg.endMjdTdb = mjd + .01;
			cfg.samples = 4;
			cfg.workers = 1;
		}

		void TearDown() override { setDefaultEphemerisSource(previous); }
		hz::IntervalMcResult run(double altitude)
		{
			cfg.altitudeKm = altitude;
			return hz::runIntervalMonteCarlo(source, forces, cfg);
		}
	};

	TEST_F(MonteCarloAltitudeTest, SameDrawsStopAtConfiguredHeightAndVelocityUsesRotatingFrame)
	{
		const auto ground = run(0);
		ASSERT_EQ(ground.entryClones, 4);
		for (double altitude : {50., 100.}) {
			const auto result = run(altitude);
			ASSERT_EQ(result.completed, 4);
			ASSERT_EQ(result.entryClones, 4);
			for (int i = 0; i < cfg.samples; ++i) {
				ASSERT_TRUE(result.clones[i].entry);
				ASSERT_TRUE(ground.clones[i].entry);
				const auto& entry = *result.clones[i].entry;
				const auto fixed = hz::icrfToBodyFixedMjd(hz::Planet::Earth, entry.relativePositionKm, entry.mjdTdb);
				EXPECT_NEAR(height(fixed), altitude, 5e-5);
				EXPECT_LT(entry.mjdTdb, ground.clones[i].entry->mjdTdb);
				const auto v = hz::bodyFixedVelocityMjd(hz::Planet::Earth, entry.relativePositionKm,
				                                         entry.relativeVelocityKmS, entry.mjdTdb);
				const double before = entry.mjdTdb - .1 / 86400, after = entry.mjdTdb + .1 / 86400;
				const auto finite = (hz::icrfToBodyFixedMjd(hz::Planet::Earth,
				    entry.relativePositionKm + entry.relativeVelocityKmS * ((after - entry.mjdTdb) * 86400), after) -
				    hz::icrfToBodyFixedMjd(hz::Planet::Earth,
				    entry.relativePositionKm + entry.relativeVelocityKmS * ((before - entry.mjdTdb) * 86400), before)) *
				    (1.0 / ((after - before) * 86400));
				EXPECT_LT((v - finite).magnitude(), 1e-5);
				EXPECT_LT(fixed.dot(v), 0);
			}
		}
	}

	TEST_F(MonteCarloAltitudeTest, EntryDoesNotRequireGroundImpact)
	{
		source.nominal.setR({1, -.0001, (Constants::EARTH_EQ_RADIUS_KM + 50) / auKm});
		source.nominal.setV({0, .02 / Constants::GAUSS_K, 0});
		EXPECT_EQ(run(0).entryClones, 0);
		EXPECT_EQ(run(100).entryClones, 4);
		source.nominal.setR({1, -.0001, (Constants::EARTH_EQ_RADIUS_KM + 150) / auKm});
		const auto misses = run(100);
		EXPECT_EQ(misses.completed, 4);
		EXPECT_EQ(misses.entryClones, 0);
	}

	TEST_F(MonteCarloAltitudeTest, MovingEarthUsesCrossingEpoch)
	{
		MovingEarth moving(source.nominal.getV() * (-3.0 * Constants::GAUSS_K));
		setDefaultEphemerisSource(&moving);
		cfg.samples = 16;
		const auto result = run(100);
		ASSERT_EQ(result.completed, 16);
		ASSERT_EQ(result.entryClones, 16);
		for (const auto& clone : result.clones) {
			ASSERT_TRUE(clone.entry);
			const auto fixed = hz::icrfToBodyFixedMjd(hz::Planet::Earth, clone.entry->relativePositionKm,
			                                           clone.entry->mjdTdb, &moving);
			EXPECT_NEAR(height(fixed), 100, 5e-5);
		}
	}

	TEST_F(MonteCarloAltitudeTest, InvalidHeightAndAlreadyInsideAreNotEntries)
	{
		EXPECT_THROW(run(-1), std::invalid_argument);
		EXPECT_THROW(run(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
		const auto result = run(1000);
		EXPECT_EQ(result.entryClones, 0);
		EXPECT_EQ(result.completed, 0);
		EXPECT_EQ(result.failedPropagation, 4);
	}

	TEST_F(MonteCarloAltitudeTest, BothIntegratorsAndWorkerCounts)
	{
		for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
			cfg.integrator = kind;
			cfg.workers = 1;
			const auto serial = run(100);
			cfg.workers = 4;
			const auto parallel = run(100);
			ASSERT_EQ(serial.completed, 4);
			ASSERT_EQ(parallel.completed, 4);
			for (int i = 0; i < 4; ++i) {
				ASSERT_TRUE(serial.clones[i].entry);
				ASSERT_TRUE(parallel.clones[i].entry);
				EXPECT_EQ(serial.clones[i].state, parallel.clones[i].state);
				EXPECT_DOUBLE_EQ(serial.clones[i].entry->mjdTdb, parallel.clones[i].entry->mjdTdb);
				const auto fixed = hz::icrfToBodyFixedMjd(hz::Planet::Earth,
				    serial.clones[i].entry->relativePositionKm, serial.clones[i].entry->mjdTdb);
				EXPECT_NEAR(height(fixed), 100, 5e-5);
			}
		}
	}
}  // namespace
