#include <gtest/gtest.h>
#include "NBodyPropagator.h"
#include "hazard/GaussianSampler.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/hazard/ObsMonteCarlo.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include <Eigen/Cholesky>
#include <array>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
	namespace hz = orbital_core::hazard;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	constexpr double epoch = 2460000.5;

	class FixedEarth : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, Vec3D(1, 0, 0), Vec3D()}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "altitude-test"; }
	};

	class MovingEarth final : public FixedEarth {
	   public:
		explicit MovingEarth(const Vec3D& velocity) : m_velocity(velocity) {}

		EphemerisRecord queryJD(EphemerisPlanet planet, double jd) override
		{
			return queryMjd(planet, orbital_core::time::mjdFromJd(jd));
		}

		EphemerisRecord queryMjd(EphemerisPlanet, double mjd) override
		{
			return {orbital_core::time::jdFromMjd(mjd),
			        Vec3D(1, 0, 0) + m_velocity * (mjd - orbital_core::time::mjdFromJd(epoch)), m_velocity};
		}

	   private:
		Vec3D m_velocity;
	};

	Vec3D fromFixed(const Vec3D& r, double jd)
	{
		return {hz::icrfToBodyFixed(hz::Planet::Earth, {1, 0, 0}, jd).dot(r),
		        hz::icrfToBodyFixed(hz::Planet::Earth, {0, 1, 0}, jd).dot(r),
		        hz::icrfToBodyFixed(hz::Planet::Earth, {0, 0, 1}, jd).dot(r)};
	}

	double height(const Vec3D& r)
	{
		return hz::referenceSurfaceSignedAltitudeKm(r, Constants::EARTH_EQ_RADIUS_KM, Constants::EARTH_FLATTENING);
	}

	class MonteCarloAltitudeTest : public ::testing::Test {
	   protected:
		FixedEarth m_eph;
		NBodyPropagator m_prop{IntegratorKind::RA15};
		Orbit m_orbit;
		hz::CloseApproachEvent m_event;
		orbital_core::StateCov6x6 m_cov = orbital_core::StateCov6x6::Identity() * 1e-22;

		void SetUp() override
		{
			m_prop.setAllPlanetsEnabled(false);
			m_prop.setUseSb441(false);
			m_prop.setStepSize(.002);
			m_orbit.setEpoch(epoch);
			const Vec3D r = fromFixed({4600, 0, 4600}, epoch);
			m_orbit.setR(Vec3D(1, 0, 0) + r * (1.0 / auKm));
			m_orbit.setV(r * (-10.0 / r.magnitude() * 86400.0 / auKm / Constants::GAUSS_K));
			m_event.planet = hz::Planet::Earth;
			m_event.jd = epoch + .002;
			m_event.rPlanet = {1, 0, 0};
			m_event.vInfKmS = 10;
		}

		hz::ObsMcResult run(double altitude, int clones = 4)
		{
			hz::ObsMcConfig cfg;
			cfg.numClones = clones;
			cfg.altitudeKm = altitude;
			return hz::runStateCovarianceMonteCarlo(m_orbit, m_cov, m_prop, {&m_event, 1}, cfg, &m_eph);
		}
	};

	TEST_F(MonteCarloAltitudeTest, SameDrawsStopAtConfiguredHeight)
	{
		const auto ground = run(0);
		ASSERT_EQ(ground.corridor.size(), 4);
		for (double altitude : {50., 100.}) {
			const auto result = run(altitude);
			EXPECT_EQ(result.earthFrameModel, hz::earthBodyFixedFrameModel());
			ASSERT_EQ(result.corridor.size(), 4);
			EXPECT_DOUBLE_EQ(result.altitudeKm, altitude);
			EXPECT_EQ(result.numFailedFit, 0);
			Eigen::Matrix<double, 6, 6> covariance = m_cov;
			const Eigen::Matrix<double, 6, 6> factor = covariance.llt().matrixL();
			for (const auto& point : result.corridor) {
				Eigen::Matrix<double, 6, 1> z;
				hz::sampling::fillStandardNormal(z, 42, point.cloneIndex);
				const auto nominalR = m_orbit.getR();
				const auto nominalV = m_orbit.getV() * Constants::GAUSS_K;
				Eigen::Matrix<double, 6, 1> nominal;
				nominal << nominalR.x, nominalR.y, nominalR.z, nominalV.x, nominalV.y, nominalV.z;
				const Eigen::Matrix<double, 6, 1> x = nominal + factor * z;
				const Vec3D r(x(0), x(1), x(2));
				const Vec3D v(x(3), x(4), x(5));
				const std::array targets{m_event.jd, m_event.jd + 1.0};
				std::vector<NBodyPropagator::StateVel> samples;
				const auto stopped = m_prop.propagateUntilSurface(
				    r, v, epoch, targets, samples, EphemerisPlanet::Earth, Constants::EARTH_EQ_RADIUS_KM, m_eph,
				    m_prop.nonGrav(), NBodyPropagator::SurfaceShape::ReferenceEllipsoid, altitude);
				ASSERT_TRUE(stopped.hitSurface);
				const auto expected =
				    hz::icrfToBodyFixed(hz::Planet::Earth, (stopped.R - Vec3D(1, 0, 0)) * auKm, stopped.jd);
				EXPECT_DOUBLE_EQ(point.tcaJd, stopped.jd);
				EXPECT_DOUBLE_EQ(point.epochMjdTdb, stopped.mjd);
				EXPECT_DOUBLE_EQ(point.rBodyFixedKm.x, expected.x);
				EXPECT_DOUBLE_EQ(point.rBodyFixedKm.y, expected.y);
				EXPECT_DOUBLE_EQ(point.rBodyFixedKm.z, expected.z);
				EXPECT_NEAR(height(point.rBodyFixedKm), altitude, 5e-5);
				EXPECT_LT(point.tcaJd, ground.corridor[point.cloneIndex].tcaJd);
				EXPECT_NEAR((point.epochMjdTdb - point.epochMjdUtc) * 86400, 69.184, .003);

				const double dt = .1 / 86400;
				const Vec3D rKm = (stopped.R - Vec3D(1, 0, 0)) * auKm;
				const Vec3D vKmDay = stopped.V * auKm;
				const double before = stopped.jd - dt, after = stopped.jd + dt;
				const auto finiteDifference =
				    (hz::icrfToBodyFixed(hz::Planet::Earth, rKm + vKmDay * (after - stopped.jd), after) -
				     hz::icrfToBodyFixed(hz::Planet::Earth, rKm + vKmDay * (before - stopped.jd), before)) *
				    (1.0 / ((after - before) * 86400));
				EXPECT_LT((point.vBodyFixedKmS - finiteDifference).magnitude(), 1e-5);
				EXPECT_LT(point.rBodyFixedKm.dot(point.vBodyFixedKmS), 0);
			}
		}
	}

	TEST_F(MonteCarloAltitudeTest, EntryDoesNotRequireGroundImpact)
	{
		m_orbit.setR({1, -.0001, (Constants::EARTH_EQ_RADIUS_KM + 50) / auKm});
		m_orbit.setV({0, .02 / Constants::GAUSS_K, 0});
		EXPECT_TRUE(run(0).corridor.empty());
		const auto entries = run(100);
		EXPECT_EQ(entries.corridor.size(), 4);
		EXPECT_EQ(entries.numFailedFit, 0);
		m_orbit.setR({1, -.0001, (Constants::EARTH_EQ_RADIUS_KM + 150) / auKm});
		const auto misses = run(100);
		EXPECT_EQ(misses.numConverged, 4);
		EXPECT_EQ(misses.numFailedFit, 0);
		EXPECT_TRUE(misses.corridor.empty());
	}

	TEST_F(MonteCarloAltitudeTest, MovingEarthUsesUnroundedCrossingEpoch)
	{
		MovingEarth eph(m_orbit.getV() * (-3.0 * Constants::GAUSS_K));
		hz::ObsMcConfig cfg;
		cfg.numClones = 16;
		cfg.altitudeKm = 100;
		const auto result = hz::runStateCovarianceMonteCarlo(m_orbit, m_cov, m_prop, {&m_event, 1}, cfg, &eph);
		ASSERT_EQ(result.numFailedFit, 0);
		ASSERT_EQ(result.corridor.size(), 16);
		bool hasSubJdPrecision = false;
		for (const auto& point : result.corridor) {
			EXPECT_NEAR(height(point.rBodyFixedKm), 100, 5e-5);
			hasSubJdPrecision |= point.epochMjdTdb != orbital_core::time::mjdFromJd(point.tcaJd);
		}
		EXPECT_TRUE(hasSubJdPrecision);
	}

	TEST_F(MonteCarloAltitudeTest, InvalidHeightAndAlreadyInsideAreNotEntries)
	{
		EXPECT_THROW(run(-1), std::invalid_argument);
		EXPECT_THROW(run(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
		EXPECT_THROW(run(std::numeric_limits<double>::infinity()), std::invalid_argument);
		const auto result = run(1000);
		EXPECT_TRUE(result.corridor.empty());
		EXPECT_EQ(result.numConverged, 0);
		EXPECT_EQ(result.numFailedFit, 4);
	}

	TEST_F(MonteCarloAltitudeTest, NonGravParameterCovarianceUsesSameBoundary)
	{
		NonGravModel ng;
		ng.setModel(NgModel::Yarkovsky);
		ng.setEnabled(true);
		orbital_core::FitParameterCovariance cov = orbital_core::FitParameterCovariance::Identity(7, 7) * 1e-30;
		hz::ObsMcConfig cfg;
		cfg.numClones = 4;
		cfg.altitudeKm = 100;
		const auto result = hz::runParameterCovarianceMonteCarlo(m_orbit, cov, ng, m_prop, {&m_event, 1}, cfg, &m_eph);
		ASSERT_EQ(result.corridor.size(), 4);
		for (const auto& point : result.corridor)
			EXPECT_NEAR(height(point.rBodyFixedKm), 100, 5e-5);
	}

	TEST_F(MonteCarloAltitudeTest, IndependentOfWorkerCount)
	{
#ifdef _OPENMP
		const int previous = omp_get_max_threads();
		omp_set_num_threads(1);
		const auto serial = run(100, 32);
		omp_set_num_threads(15);
		const auto parallel = run(100, 32);
		omp_set_num_threads(previous);
		ASSERT_EQ(serial.corridor.size(), 32);
		ASSERT_EQ(serial.corridor.size(), parallel.corridor.size());
		for (std::size_t i = 0; i < serial.corridor.size(); ++i) {
			const auto& a = serial.corridor[i];
			const auto& b = parallel.corridor[i];
			EXPECT_EQ(a.cloneIndex, b.cloneIndex);
			EXPECT_DOUBLE_EQ(a.tcaJd, b.tcaJd);
			EXPECT_DOUBLE_EQ(a.rBodyFixedKm.x, b.rBodyFixedKm.x);
			EXPECT_DOUBLE_EQ(a.rBodyFixedKm.y, b.rBodyFixedKm.y);
			EXPECT_DOUBLE_EQ(a.rBodyFixedKm.z, b.rBodyFixedKm.z);
			EXPECT_DOUBLE_EQ(a.vBodyFixedKmS.x, b.vBodyFixedKmS.x);
			EXPECT_DOUBLE_EQ(a.vBodyFixedKmS.y, b.vBodyFixedKmS.y);
			EXPECT_DOUBLE_EQ(a.vBodyFixedKmS.z, b.vBodyFixedKmS.z);
		}
#endif
	}

	TEST_F(MonteCarloAltitudeTest, BothIntegratorsUseEllipsoidalHeightAtEquatorAndPole)
	{
		for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
			m_prop.setIntegratorKind(kind);
			for (double latitude : {0., 45., 90.}) {
				SCOPED_TRACE(latitude);
				const double phi = latitude * Constants::DEG_TO_RAD;
				const double e2 = Constants::EARTH_FLATTENING * (2 - Constants::EARTH_FLATTENING);
				const double n = Constants::EARTH_EQ_RADIUS_KM / std::sqrt(1 - e2 * std::sin(phi) * std::sin(phi));
				const Vec3D r = fromFixed({(n + 200) * std::cos(phi), 0, (n * (1 - e2) + 200) * std::sin(phi)}, epoch);
				m_orbit.setR(Vec3D(1, 0, 0) + r * (1.0 / auKm));
				m_orbit.setV(r * (-10.0 / r.magnitude() * 86400.0 / auKm / Constants::GAUSS_K));
				const auto result = run(100);
				ASSERT_EQ(result.corridor.size(), 4);
				for (const auto& point : result.corridor) {
					EXPECT_NEAR(height(point.rBodyFixedKm), 100, 5e-5);
					if (latitude == 90)
						EXPECT_NEAR(point.rBodyFixedKm.magnitude(),
						            Constants::EARTH_EQ_RADIUS_KM * (1 - Constants::EARTH_FLATTENING) + 100, .001);
				}
			}
		}
	}
}  // namespace
