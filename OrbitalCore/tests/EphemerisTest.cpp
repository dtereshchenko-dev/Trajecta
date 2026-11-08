#include <gtest/gtest.h>

#include "orbital_core/Constants.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Ephemeris.h"
#include "orbital_core/FitSession.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/Types.h"
#include "obs/LightDeflection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <memory>
#include <limits>
#include <stop_token>
#include <string>
#include <vector>

namespace oc = orbital_core;

namespace {

	constexpr double kArcsec = Constants::PI / (180.0 * 3600.0);

	[[nodiscard]] bool resolveDataDir()
	{
		if (!getExecutableDir().empty())
			return true;
		const auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
		for (const auto& dir : {srcDir.parent_path() / "data", srcDir / "data", srcDir / "release", srcDir / "debug"}) {
			if (std::filesystem::exists(dir / "linux_p1550p2650.440")) {
				getExecutableDir() = dir;
				return true;
			}
		}
		return false;
	}

	// 2023 BU, JPL solution #23 (Horizons, heliocentric ecliptic J2000 ICRF, AU and
	// AU/day) at 2023-01-20 00:00 TDB, a week before its 9968 km geocentric pass.
	const Vec3D kBuPosition(-4.887497265974482E-01, 8.649413611927270E-01, 5.352552184023343E-03);
	const Vec3D kBuVelocity(-1.437612291038726E-02, -9.412434053643936E-03, -7.248763353981187E-04);
	constexpr double kBuEpochMjdTdb = 59964.0;

	[[nodiscard]] double separationArcsec(double ra1, double dec1, double ra2, double dec2) noexcept
	{
		const double sinDec = std::sin(0.5 * (dec2 - dec1));
		const double sinRa = std::sin(0.5 * (ra2 - ra1));
		const double h = sinDec * sinDec + std::cos(dec1) * std::cos(dec2) * sinRa * sinRa;
		return 2.0 * std::asin(std::sqrt(std::clamp(h, 0.0, 1.0))) / kArcsec;
	}

}  // namespace

TEST(TimeScaleTest, UtcTtAndTdbConversionsRoundTrip)
{
	const double mjdUtc = 59000.25;
	const double tt = oc::time::toMjdTt(mjdUtc, oc::time::Scale::UTC);
	EXPECT_NEAR((tt - mjdUtc) * 86400.0, 69.184, 1.0e-6);
	EXPECT_NEAR(oc::time::fromMjdTt(tt, oc::time::Scale::UTC), mjdUtc, 1.0e-11);

	const double tdb = oc::time::fromMjdTt(tt, oc::time::Scale::TDB);
	EXPECT_DOUBLE_EQ(oc::time::toMjdTdb(mjdUtc, oc::time::Scale::UTC), tdb);
	EXPECT_DOUBLE_EQ(oc::time::toMjdTdb(tdb, oc::time::Scale::TDB), tdb);
	EXPECT_DOUBLE_EQ(oc::time::fromMjdTdb(tdb, oc::time::Scale::TDB), tdb);
	EXPECT_NEAR(oc::time::fromMjdTdb(tdb, oc::time::Scale::UTC), mjdUtc, 1.0e-11);
	EXPECT_NEAR(oc::time::fromMjdTdb(tdb, oc::time::Scale::TT), tt, 1.0e-11);
	EXPECT_NEAR((tdb - tt) * 86400.0, Calculations::tdbMinusTtSeconds(tt), 2.0e-6);  // MJD ulp ~0.6 us
	EXPECT_LT(std::abs((tdb - tt) * 86400.0), 1.8e-3);
	EXPECT_NEAR(oc::time::toMjdTt(tdb, oc::time::Scale::TDB), tt, 1.0e-12);
	EXPECT_DOUBLE_EQ(oc::time::toMjdTt(tt, oc::time::Scale::TT), tt);

	const double beforeLeap = 57753.5;  // 2016-12-31 12:00 UTC, before the 37th second
	EXPECT_NEAR((oc::time::toMjdTt(beforeLeap, oc::time::Scale::UTC) - beforeLeap) * 86400.0, 68.184, 1.0e-6);
}

TEST(EphemerisGridTest, EndIsInclusiveAndInvalidGridsThrow)
{
	EXPECT_EQ(oc::ephemerisEpochCount(59000.0, 59014.0, 1.0), 15u);
	EXPECT_EQ(oc::ephemerisEpochCount(59000.0, 59000.0, 1.0), 1u);
	EXPECT_EQ(oc::ephemerisEpochCount(59000.0, 59001.0 - 5.0e-7, 0.5), 2u);
	EXPECT_THROW((void)oc::ephemerisEpochCount(59000.0, 59001.0, 0.0), std::invalid_argument);
	EXPECT_THROW((void)oc::ephemerisEpochCount(59001.0, 59000.0, 1.0), std::invalid_argument);
	EXPECT_GT(oc::ephemerisEpochCount(0.0, 1.0e6, 1.0e-3), oc::kMaxEphemerisEpochs);
}

TEST(EphemerisGridTest, FineStepsNeverAddEpochsBeyondTheEnd)
{
	const double start = 59000.0;
	for (const double step : {0.1, 1e-5, 1e-7, 1e-9}) {
		for (const std::size_t intervals : {1u, 2u, 100u}) {
			const double end = start + static_cast<double>(intervals) * step;
			EXPECT_EQ(oc::ephemerisEpochCount(start, end, step), intervals + 1);
			EXPECT_EQ(oc::ephemerisEpochCount(start, std::nextafter(end, start), step), intervals);
		}
	}
	const double ulp = std::nextafter(start, std::numeric_limits<double>::infinity()) - start;
	EXPECT_EQ(oc::ephemerisEpochCount(start, start, ulp / 2), 1u);
	EXPECT_EQ(oc::ephemerisEpochCount(start, start + 10 * ulp, ulp), 11u);
	EXPECT_THROW((void)oc::ephemerisEpochCount(start, start + 10 * ulp, ulp / 2), std::invalid_argument);
	EXPECT_EQ(oc::ephemerisEpochCount(0, oc::kMaxEphemerisEpochs - 1, 1), oc::kMaxEphemerisEpochs);
	EXPECT_EQ(oc::ephemerisEpochCount(0, oc::kMaxEphemerisEpochs, 1), oc::kMaxEphemerisEpochs + 1);
}

class EphemerisTest : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (!resolveDataDir())
			GTEST_SKIP() << "data/ directory not found";
		m_de440 = std::make_unique<DE440Source>();
		if (!m_de440->isAvailable())
			GTEST_SKIP() << "DE440 kernel not loadable";
		setDefaultEphemerisSource(m_de440.get());
		m_sb441 = std::make_unique<SpkReader>();
		if (m_sb441->isAvailable())
			setDefaultAsteroidSource(m_sb441.get());
	}

	void TearDown() override
	{
		setDefaultEphemerisSource(nullptr);
		setDefaultAsteroidSource(nullptr);
	}

	[[nodiscard]] oc::EphemerisRequest buRequest(oc::EphemerisCenter center, const char* station = "") const
	{
		oc::EphemerisRequest request;
		request.solution.position = kBuPosition;
		request.solution.velocity = kBuVelocity;
		request.solution.epochMjdTdb = kBuEpochMjdTdb;
		request.solution.forces.enabled = true;
		request.solution.forces.useSb441 = m_sb441 && m_sb441->isAvailable();
		request.solution.absoluteMagnitude = 29.69;
		request.center = center;
		request.station = station;
		return request;
	}

	std::unique_ptr<DE440Source> m_de440;
	std::unique_ptr<SpkReader> m_sb441;
};

TEST_F(EphemerisTest, ReproducesTheFitObservationModelExactly)
{
	const auto initial =
	    make_orbit_from_states_vectors(StateVectors(kBuPosition, kBuVelocity / Constants::GAUSS_K), 2459964.5);
	ASSERT_TRUE(initial);

	oc::FitConfig config;
	config.nbody.enabled = true;
	config.sigmaClip.enabled = false;
	const oc::FitSession session{config};
	const oc::FitRunOptions evaluate{
	    .selection = oc::FitSelectionMode::Fixed, .uncertainty = oc::FitUncertaintyMode::Skip, .evaluateOnly = true};

	for (const char* station : {"F51", "W68", "500"}) {
		oc::EphemerisRequest request;
		request.solution.forces = config.nbody;
		request.center = oc::EphemerisCenter::Station;
		request.station = station;
		request.startMjd = 59966.0;
		request.endMjd = 59971.5;
		request.stepDays = 0.25;
		const auto observe = [&](const Orbit& orbit) {
			request.solution.position = orbit.getR();
			request.solution.velocity = orbit.getV() * Constants::GAUSS_K;
			request.solution.epochMjdTdb = orbit.epochMjd();
			const auto ephemeris = oc::computeEphemeris(request);
			std::vector<ObsNode> rows(ephemeris.rows.size());
			for (std::size_t i = 0; i < rows.size(); ++i) {
				rows[i].JD = oc::time::jdFromMjd(ephemeris.rows[i].mjdUtc);
				const auto& row = ephemeris.rows[i];
				const double c = std::cos(Constants::EARTH_INCL), s = std::sin(Constants::EARTH_INCL);
				const Vec3D topo(row.position.x, row.position.y * c + row.position.z * s,
				                 -row.position.y * s + row.position.z * c);
				// The station offset is negligible for this correction at the test's tolerance.
				const Vec3D observer = m_de440->getPlanetPositionMjd(EphemerisPlanet::Earth, row.mjdTdb);
				const Vec3D corrected = oc::obs::catalogReferencedDirection(topo, topo + observer, observer);
				const Vec3D eq(corrected.x, corrected.y * c - corrected.z * s, corrected.y * s + corrected.z * c);
				rows[i].RA = std::atan2(eq.y, eq.x);
				if (rows[i].RA < 0.0)
					rows[i].RA += 2.0 * Constants::PI;
				rows[i].DEC = std::atan2(eq.z, std::hypot(eq.x, eq.y));
				rows[i].obsCode = station;
			}
			return rows;
		};

		// The session evaluates from its own fit epoch; predict from that state.
		const auto first = observe(*initial);
		ASSERT_EQ(first.size(), 23u);
		const auto rebased = session.optimize({.initialOrbit = initial, .observations = first, .options = evaluate});
		ASSERT_TRUE(rebased.fitEpochOrbit);
		const auto observations = observe(*rebased.fitEpochOrbit);
		const auto fit = session.optimize(
		    {.initialOrbit = rebased.fitEpochOrbit, .observations = observations, .options = evaluate});
		ASSERT_EQ(fit.residuals.count(), observations.size());
		double maxResidual = 0.0;
		for (std::size_t i = 0; i < observations.size(); ++i) {
			const auto& row = fit.residuals.data()[i];
			maxResidual = std::max({maxResidual, std::abs(row.dRA_arcsec), std::abs(row.dDEC_arcsec)});
			// The ephemeris uses the geometric light-time direction; synthetic
			// observations above use the catalogue-referenced convention.
			// The fit also carries the STM and uses a different step sequence.
			EXPECT_LT(std::abs(row.dRA_arcsec), 1.0e-4) << station << " row " << i;
			EXPECT_LT(std::abs(row.dDEC_arcsec), 1.0e-4) << station << " row " << i;
		}
		RecordProperty(std::string("max_residual_arcsec_") + station, std::to_string(maxResidual));
	}
}

TEST_F(EphemerisTest, LightTimeRatesDifferentiatePositionsAcrossCloseApproach)
{
	for (const auto center :
	     {oc::EphemerisCenter::Heliocentric, oc::EphemerisCenter::Geocentric, oc::EphemerisCenter::Station}) {
		for (bool lightTime : {false, true}) {
			auto request = buRequest(center, "W68");
			request.scale = oc::time::Scale::TT;
			request.lightTime = lightTime;
			const double epoch = oc::time::toMjdTt(59970.97843, oc::time::Scale::UTC);
			request.stepDays = 0.5 / 86400.0;
			request.startMjd = epoch - request.stepDays;
			request.endMjd = request.startMjd + 2 * request.stepDays;
			const auto rows = oc::computeEphemeris(request).rows;
			ASSERT_EQ(rows.size(), 3u);
			const double dt = rows[2].mjdTt - rows[0].mjdTt;
			const auto measured = (rows[2].position - rows[0].position) / dt;
			EXPECT_LT((rows[1].velocity - measured).magnitude() * Constants::AU_IN_METERS / 86400.0, 1e-3);
			const double raRate =
			    Helpers::wrapPi(rows[2].longitude - rows[0].longitude) * std::cos(rows[1].latitude) / dt;
			const double decRate = (rows[2].latitude - rows[0].latitude) / dt;
			EXPECT_NEAR(rows[1].longitudeRate, raRate, 5e-4 * kArcsec * 1440);
			EXPECT_NEAR(rows[1].latitudeRate, decRate, 5e-4 * kArcsec * 1440);
			if (!lightTime || center == oc::EphemerisCenter::Heliocentric) {
				EXPECT_NEAR(rows[1].distanceRateAuPerDay, (rows[2].distanceAu - rows[0].distanceAu) / dt, 1.0e-8);
				EXPECT_NEAR(rows[1].heliocentricRateAuPerDay, (rows[2].heliocentricAu - rows[0].heliocentricAu) / dt,
				            1.0e-8);
			}
			EXPECT_DOUBLE_EQ(rows[1].oneWayLightTimeDays, rows[1].distanceAu / Constants::AU_PER_DAY);
		}
	}
}

TEST_F(EphemerisTest, LightTimePhotometryGeometryAndRadialVelocitiesMatchReference)
{
	struct Reference {
		Vec3D position;
		Vec3D velocity;
		const char* station;
		std::array<double, 3> solarRange;
		std::array<double, 3> solarRateKmS;
		std::array<double, 3> observerRange;
		std::array<double, 3> observerRateKmS;
		std::array<double, 3> phaseDegrees;
	};

	const Vec3D chesleyPos(2.598867861350580, -1.390029528986849, -0.5576174308536911);
	const Vec3D chesleyVel(0.004807172863533565, 0.008695388694768491, -0.0005697612003221255);
	const std::array<Reference, 4> references{{
	    {chesleyPos,
	     chesleyVel,
	     "500",
	     {2.999534849595, 2.999776279113, 3.000017820655},
	     {0.4179589, 0.4181551, 0.4183466},
	     {3.10468362788685, 3.11818760404627, 3.13164466940049},
	     {23.4204475, 23.3412047, 23.2580206},
	     {18.5221, 18.4483, 18.3718}},
	    {chesleyPos,
	     chesleyVel,
	     "033",
	     {2.999534849562, 2.999776279079, 3.000017820620},
	     {0.4179589, 0.4181551, 0.4183466},
	     {3.10470776376887, 3.11821196625990, 3.13166925549799},
	     {23.6618044, 23.5815134, 23.4972472},
	     {18.5226, 18.4489, 18.3723}},
	    {Vec3D(-0.07852649149069046, -0.8197480519020646, 0.04189395153233909),
	     Vec3D(0.01987510249688846, 0.001322088445361402, 0.0003994960444223522),
	     "500",
	     {0.824595526480, 0.821425534501, 0.818291450626},
	     {-5.5184387, -5.4580785, -5.3940886},
	     {1.76617102756837, 1.76160557032962, 1.75702926329873},
	     {-7.8946213, -7.9141824, -7.9321729},
	     {14.2409, 14.4648, 14.6947}},
	    {Vec3D(0.8031176726254730, 0.9930978829545563, 0.2339101254563577),
	     Vec3D(-0.01393305158641565, 0.007535949457489370, -0.001390706422743106),
	     "500",
	     {1.298448154570, 1.295349430869, 1.292263006173},
	     {-5.3756423, -5.3547898, -5.3330552},
	     {0.40016966952235, 0.39962166129306, 0.39914726801918},
	     {-1.0140058, -0.8843426, -0.7591036},
	     {33.1581, 33.6420, 34.1299}},
	}};
	constexpr double auPerDayToKmS = 149597870.7 / 86400.0;
	for (const auto& reference : references) {
		auto request = buRequest(oc::EphemerisCenter::Station, reference.station);
		request.solution.position = reference.position;
		request.solution.velocity = reference.velocity;
		request.solution.epochMjdTdb = 61000.0;
		request.startMjd = 61000.0;
		request.endMjd = 61002.0;
		request.stepDays = 1.0;
		const auto rows = oc::computeEphemeris(request).rows;
		ASSERT_EQ(rows.size(), 3u);
		for (std::size_t i = 0; i < rows.size(); ++i) {
			EXPECT_NEAR(rows[i].heliocentricAu, reference.solarRange[i], 3.0e-9) << reference.station << " " << i;
			EXPECT_NEAR(rows[i].heliocentricRateAuPerDay * auPerDayToKmS, reference.solarRateKmS[i], 1.0e-5);
			EXPECT_NEAR(rows[i].distanceAu, reference.observerRange[i], 3.0e-9);
			EXPECT_NEAR(rows[i].distanceRateAuPerDay * auPerDayToKmS, reference.observerRateKmS[i], 1.0e-5);
			EXPECT_NEAR(rows[i].phase * Constants::RAD_TO_DEG, reference.phaseDegrees[i], 1.5e-4);
		}
	}
}

TEST_F(EphemerisTest, ChesleyObserverFieldsAgainstHorizonsAtSharedEpoch)
{
	// Horizons JPL#60, geocentric, 2025-11-21..23 00:00 UTC. The initial
	// heliocentric ecliptic state is Horizons at MJD 61000 TDB; Horizons uses
	// DE441 and SB441, while this test uses the locally available force model.
	struct Reference {
		double ra, dec, magnitude, solarDistance, solarRate, distance, distanceRate, lightMinutes, elongation, phase;
		double raBefore, raAfter, decBefore, decAfter;
	};
	constexpr std::array<Reference, 3> reference{{
	    {318.751209037, -26.800558241, 18.497, 2.999534849596, 0.4179589, 3.10468362788899, 23.4204475,
	     25.82086638, 74.7121, 18.5221, 318.751021552, 318.751396524, -26.800637075, -26.800479406},
	    {319.022563162, -26.686605948, 18.505, 2.999776279114, 0.4181551, 3.11818760404841, 23.3412047,
	     25.93317552, 73.9829, 18.4483, 319.022373775, 319.022752550, -26.686685379, -26.686526516},
	    {319.296625170, -26.571798263, 18.512, 3.000017820656, 0.4183466, 3.13164466940263, 23.2580206,
	     26.04509452, 73.2562, 18.3718, 319.296433923, 319.296816417, -26.571878286, -26.571718240},
	}};
	auto request = buRequest(oc::EphemerisCenter::Geocentric);
	request.solution.position = Vec3D(2.598867861350580, -1.390029528986849, -0.5576174308536911);
	request.solution.velocity = Vec3D(0.004807172863533565, 0.008695388694768491, -0.0005697612003221255);
	request.solution.epochMjdTdb = 61000.0;
	request.solution.absoluteMagnitude = 12.7;
	request.solution.slopeG = 0.15;
	request.startMjd = 61000.0;
	request.endMjd = 61002.0;
	request.stepDays = 1.0;
	const auto rows = oc::computeEphemeris(request).rows;
	ASSERT_EQ(rows.size(), reference.size());
	constexpr double kmPerSecond = 149597870.7 / 86400.0;
	constexpr double arcsecondsPerHour = 206264.80624709636 / 24.0;
	for (std::size_t i = 0; i < rows.size(); ++i) {
		const auto& row = rows[i];
		const auto& ref = reference[i];
		const double skyError = separationArcsec(row.longitude, row.latitude, ref.ra * Constants::DEG_TO_RAD,
		                                         ref.dec * Constants::DEG_TO_RAD);
		const double referenceRaRate = (ref.raAfter - ref.raBefore) * 30.0 * 3600.0 * std::cos(row.latitude);
		const double referenceDecRate = (ref.decAfter - ref.decBefore) * 30.0 * 3600.0;
		const double raRate = row.longitudeRate * arcsecondsPerHour;
		const double decRate = row.latitudeRate * arcsecondsPerHour;
		const double motion = std::hypot(raRate, decRate);
		const double referenceMotion = std::hypot(referenceRaRate, referenceDecRate);
		const double positionAngle = std::atan2(raRate, decRate) * Constants::RAD_TO_DEG;
		const double referencePositionAngle = std::atan2(referenceRaRate, referenceDecRate) * Constants::RAD_TO_DEG;
		ASSERT_TRUE(row.visualMagnitude);
		RecordProperty("sky_arcsec_" + std::to_string(i), std::to_string(skyError));
		RecordProperty("magnitude_" + std::to_string(i), std::to_string(*row.visualMagnitude - ref.magnitude));
		RecordProperty("solar_range_au_" + std::to_string(i), std::to_string(row.heliocentricAu - ref.solarDistance));
		RecordProperty("solar_rate_km_s_" + std::to_string(i),
		               std::to_string(row.heliocentricRateAuPerDay * kmPerSecond - ref.solarRate));
		RecordProperty("observer_range_au_" + std::to_string(i), std::to_string(row.distanceAu - ref.distance));
		RecordProperty("observer_rate_km_s_" + std::to_string(i),
		               std::to_string(row.distanceRateAuPerDay * kmPerSecond - ref.distanceRate));
		RecordProperty("light_time_min_" + std::to_string(i),
		               std::to_string(row.oneWayLightTimeDays * 1440.0 - ref.lightMinutes));
		RecordProperty("ra_rate_arcsec_h_" + std::to_string(i),
		               std::to_string(row.longitudeRate * arcsecondsPerHour - referenceRaRate));
		RecordProperty("dec_rate_arcsec_h_" + std::to_string(i),
		               std::to_string(row.latitudeRate * arcsecondsPerHour - referenceDecRate));
		RecordProperty("sky_motion_arcsec_h_" + std::to_string(i), std::to_string(motion - referenceMotion));
		RecordProperty("sky_motion_pa_deg_" + std::to_string(i), std::to_string(positionAngle - referencePositionAngle));
		RecordProperty("elongation_deg_" + std::to_string(i),
		               std::to_string(row.elongation * Constants::RAD_TO_DEG - ref.elongation));
		RecordProperty("phase_deg_" + std::to_string(i),
		               std::to_string(row.phase * Constants::RAD_TO_DEG - ref.phase));
		EXPECT_LT(skyError, 0.02);
		EXPECT_NEAR(*row.visualMagnitude, ref.magnitude, 0.01);
		EXPECT_NEAR(row.heliocentricAu, ref.solarDistance, 3e-9);
		EXPECT_NEAR(row.heliocentricRateAuPerDay * kmPerSecond, ref.solarRate, 1e-5);
		EXPECT_NEAR(row.distanceAu, ref.distance, 3e-9);
		EXPECT_NEAR(row.distanceRateAuPerDay * kmPerSecond, ref.distanceRate, 1e-5);
		EXPECT_NEAR(row.oneWayLightTimeDays * 1440.0, ref.lightMinutes, 3e-8);
		EXPECT_NEAR(row.elongation * Constants::RAD_TO_DEG, ref.elongation, 0.0001);
		// Horizons quantity 43 is the true phase angle, not quantity 24 (S-T-O).
		EXPECT_NEAR(row.phase * Constants::RAD_TO_DEG, ref.phase, 0.0001);
		EXPECT_NEAR(row.longitudeRate * arcsecondsPerHour, referenceRaRate, 0.001);
		EXPECT_NEAR(row.latitudeRate * arcsecondsPerHour, referenceDecRate, 0.001);
		EXPECT_NEAR(motion, referenceMotion, 0.001);
		EXPECT_NEAR(positionAngle, referencePositionAngle, 0.001);
	}
	request.center = oc::EphemerisCenter::Station;
	request.station = "033";
	const auto stationRows = oc::computeEphemeris(request).rows;
	ASSERT_EQ(stationRows.size(), 3u);
	constexpr std::array<double, 3> stationElongation{74.7107, 73.9815, 73.2548};
	for (std::size_t i = 0; i < stationRows.size(); ++i)
		EXPECT_NEAR(stationRows[i].elongation * Constants::RAD_TO_DEG, stationElongation[i], 0.0001);

	// Geometric equatorial J2000 vector at the same initial TDB epoch.
	request.center = oc::EphemerisCenter::Geocentric;
	request.scale = oc::time::Scale::TDB;
	request.startMjd = request.endMjd = 61000.0;
	request.lightTime = false;
	const auto geometric = oc::computeEphemeris(request).rows;
	ASSERT_EQ(geometric.size(), 1u);
	const Vec3D expectedPosition(2.083594703397822, -1.826976819216022, -1.399805990259024);
	const Vec3D expectedVelocity(0.01977218913068487, 0.00002906018281776681, -0.0006080047842801818);
	EXPECT_LT((geometric[0].position - expectedPosition).magnitude(), 3e-9);
	EXPECT_LT((geometric[0].velocity - expectedVelocity).magnitude(), 3e-10);
	const auto earth = m_de440->queryMjd(EphemerisPlanet::Earth, 61000.0);
	const double c = std::cos(Constants::EARTH_INCL), s = std::sin(Constants::EARTH_INCL);
	const Vec3D earthEq(earth.R.x, earth.R.y * c - earth.R.z * s, earth.R.y * s + earth.R.z * c);
	const double geometricElongation = std::acos(std::clamp(
	    (earthEq * -1.0).dot(geometric[0].position) / (earthEq.magnitude() * geometric[0].position.magnitude()), -1.0, 1.0));
	EXPECT_NEAR(geometric[0].elongation, geometricElongation, 1e-12);
}

TEST_F(EphemerisTest, SkyUncertaintyMatchesPerturbedEphemerides)
{
	for (const bool nBody : {false, true}) {
		for (const bool lightTime : {false, true}) {
		  for (const double epoch : {59966.0, 61000.0}) {
			auto request = buRequest(oc::EphemerisCenter::Station, "033");
			request.solution.forces.enabled = nBody;
			request.solution.forces.useSb441 = false;
			request.lightTime = lightTime;
			request.startMjd = request.endMjd = epoch;
			request.uncertainty = true;
			Eigen::Matrix<double, 6, 6> factor = Eigen::Matrix<double, 6, 6>::Zero();
			factor.diagonal() << 1e-8, 2e-8, 3e-8, 2e-10, 3e-10, 4e-10;
			factor(0, 4) = 3e-9;
			factor(2, 3) = -2e-9;
			request.solution.covariance = factor * factor.transpose();
			const auto nominal = oc::computeEphemeris(request).rows.at(0);
			ASSERT_TRUE(nominal.uncertainty);
			Eigen::Matrix<double, 4, 6> jacobian;
			request.uncertainty = false;
			for (int i = 0; i < 6; ++i) {
				auto plus = request, minus = request;
				const double step = i < 3 ? 1e-6 : 1e-8;
				Vec3D& p = i < 3 ? plus.solution.position : plus.solution.velocity;
				Vec3D& m = i < 3 ? minus.solution.position : minus.solution.velocity;
				double* pc[] = {&p.x, &p.y, &p.z};
				double* mc[] = {&m.x, &m.y, &m.z};
				*pc[i % 3] += step;
				*mc[i % 3] -= step;
				const auto a = oc::computeEphemeris(plus).rows.at(0);
				const auto b = oc::computeEphemeris(minus).rows.at(0);
				jacobian(0, i) = std::remainder(a.longitude - b.longitude, 2.0 * Constants::PI) *
				                 std::cos(nominal.latitude) / (2.0 * step);
				jacobian(1, i) = (a.latitude - b.latitude) / (2.0 * step);
				jacobian(2, i) = (a.distanceAu - b.distanceAu) / (2.0 * step);
				jacobian(3, i) = (a.distanceRateAuPerDay - b.distanceRateAuPerDay) / (2.0 * step);
			}
			const Eigen::Matrix4d expected = jacobian * request.solution.covariance * jacobian.transpose();
			constexpr double arcsec = 206264.80624709636;
			const auto& actual = *nominal.uncertainty;
			const double east = 3.0 * std::sqrt(expected(0, 0)) * arcsec;
			const double north = 3.0 * std::sqrt(expected(1, 1)) * arcsec;
			const double range = 3.0 * std::sqrt(expected(2, 2)) * (Constants::AU_IN_METERS / 1000.0);
			EXPECT_NEAR(actual.longitude3SigmaArcsec, east, east * 2e-4);
			EXPECT_NEAR(actual.latitude3SigmaArcsec, north, north * 2e-4);
			EXPECT_NEAR(actual.range3SigmaKm, range, range * 2e-4);
			const double rangeRate = 3.0 * std::sqrt(expected(3, 3)) * (Constants::AU_IN_METERS / 1000.0) / 86400.0;
			EXPECT_NEAR(actual.rangeRate3SigmaKmPerSecond, rangeRate, rangeRate * 3e-4);
			EXPECT_NEAR(actual.major3SigmaArcsec * actual.major3SigmaArcsec +
			                actual.minor3SigmaArcsec * actual.minor3SigmaArcsec,
			            east * east + north * north, (east * east + north * north) * 2e-4);
		  }
		}
	}
}

TEST_F(EphemerisTest, ChesleySkyUncertaintyAgreesWithHorizons)
{
	// SBDB JPL#60 covariance in [e, q, tp, node, peri, i], at JD 2458383.5 TDB.
	using Vector6 = Eigen::Matrix<double, 6, 1>;
	using Matrix6 = Eigen::Matrix<double, 6, 6>;
	const Vector6 elements(0.02466299754368284, 2.936562809597262, 2458725.988367522770,
	                       78.06106008819181, 182.5320169349131, 11.15274171364406);
	Matrix6 covariance;
	covariance << 5.341658843758314e-19, -2.413560534201591e-18, 2.137674205518474e-15, -4.814053735820984e-18, 4.73977827543846e-16, -1.06909165126936e-18,
	    -2.413560534201591e-18, 1.159221689362191e-17, -1.18727956270191e-14, 2.482937430756793e-17, -2.609689576634236e-15, 3.844903204503517e-18,
	    2.137674205518474e-15, -1.18727956270191e-14, 4.529949217506534e-11, 4.619217982682939e-14, 8.900935638769985e-12, 1.58267220373105e-15,
	    -4.814053735820984e-18, 2.482937430756793e-17, 4.619217982682939e-14, 1.575411787756052e-14, -7.885564600784617e-15, -4.784858482750135e-16,
	    4.73977827543846e-16, -2.609689576634236e-15, 8.900935638769985e-12, -7.885564600784617e-15, 1.776271207268168e-12, 7.212448314149736e-16,
	    -1.06909165126936e-18, 3.844903204503517e-18, 1.58267220373105e-15, -4.784858482750135e-16, 7.212448314149736e-16, 8.893618693010124e-16;
	const auto state = [](const Vector6& p) -> Vector6 {
		const double e = p[0], a = p[1] / (1.0 - e);
		const double mean = Constants::GAUSS_K * (2458383.5 - p[2]) / std::pow(a, 1.5);
		double eccentric = mean;
		for (int i = 0; i < 12; ++i)
			eccentric -= (eccentric - e * std::sin(eccentric) - mean) / (1.0 - e * std::cos(eccentric));
		const double anomaly = std::atan2(std::sqrt(1.0 - e * e) * std::sin(eccentric), std::cos(eccentric) - e);
		const oc::cometary::Vector6 cot(p[1], e, p[5] * Constants::DEG_TO_RAD, p[3] * Constants::DEG_TO_RAD,
		                               p[4] * Constants::DEG_TO_RAD, anomaly);
		return oc::cometary::toCartesian(cot).value();
	};
	const Vector6 nominal = state(elements);
	Matrix6 transform;
	const Vector6 steps(1e-5, 1e-4, 0.005, 0.001, 0.001, 0.001);
	for (int i = 0; i < 6; ++i) {
		Vector6 plus = elements, minus = elements;
		plus[i] += steps[i];
		minus[i] -= steps[i];
		transform.col(i) = (state(plus) - state(minus)) / (plus[i] - minus[i]);
	}
	auto request = buRequest(oc::EphemerisCenter::Geocentric);
	request.solution.position = Vec3D(nominal[0], nominal[1], nominal[2]);
	request.solution.velocity = Vec3D(nominal[3], nominal[4], nominal[5]);
	request.solution.epochMjdTdb = 58383.0;
	request.solution.covariance = transform * covariance * transform.transpose();
	request.scale = oc::time::Scale::TT;
	request.startMjd = request.endMjd = 58383.0;
	request.uncertainty = true;
	const auto row = oc::computeEphemeris(request).rows.at(0);
	ASSERT_TRUE(row.uncertainty);
	EXPECT_NEAR(row.uncertainty->major3SigmaArcsec, 0.00048, 0.00002);
	EXPECT_NEAR(row.uncertainty->minor3SigmaArcsec, 0.00023, 0.00002);
	EXPECT_NEAR(row.uncertainty->positionAngleDeg, -26.291, 2.0);
	EXPECT_NEAR(row.uncertainty->range3SigmaKm, 0.7874, 0.03);
	request.scale = oc::time::Scale::UTC;
	request.startMjd = request.endMjd = 61000.0;
	const auto propagated = oc::computeEphemeris(request).rows.at(0);
	ASSERT_TRUE(propagated.uncertainty);
	EXPECT_NEAR(propagated.uncertainty->major3SigmaArcsec, 0.00444, 0.0001);
	EXPECT_NEAR(propagated.uncertainty->minor3SigmaArcsec, 0.00032, 0.00002);
	EXPECT_NEAR(propagated.uncertainty->positionAngleDeg, 10.822, 2.0);
	EXPECT_NEAR(propagated.uncertainty->range3SigmaKm, 4.5566, 0.15);
	request.uncertainty = false;
	Eigen::Matrix<double, 1, 6> rangeRateJacobian;
	for (int i = 0; i < 6; ++i) {
		auto plus = request, minus = request;
		const double step = i < 3 ? 1e-6 : 1e-8;
		Vec3D& p = i < 3 ? plus.solution.position : plus.solution.velocity;
		Vec3D& m = i < 3 ? minus.solution.position : minus.solution.velocity;
		double* pc[] = {&p.x, &p.y, &p.z};
		double* mc[] = {&m.x, &m.y, &m.z};
		*pc[i % 3] += step;
		*mc[i % 3] -= step;
		rangeRateJacobian[i] = (oc::computeEphemeris(plus).rows.at(0).distanceRateAuPerDay -
		                        oc::computeEphemeris(minus).rows.at(0).distanceRateAuPerDay) / (2.0 * step);
	}
	const double numericSigma = 3.0 * std::sqrt((rangeRateJacobian * request.solution.covariance *
	                                             rangeRateJacobian.transpose())(0, 0)) *
	                            (Constants::AU_IN_METERS / 1000.0) / 86400.0;
	EXPECT_NEAR(propagated.uncertainty->rangeRate3SigmaKmPerSecond, numericSigma, numericSigma * 2e-3);
}

TEST_F(EphemerisTest, SkyUncertaintyIncludesNonGravitationalCrossCovariance)
{
	auto request = buRequest(oc::EphemerisCenter::Geocentric);
	request.solution.forces.useSb441 = false;
	request.solution.position = Vec3D(2.598867861350580, -1.390029528986849, -0.5576174308536911);
	request.solution.velocity = Vec3D(0.004807172863533565, 0.008695388694768491, -0.0005697612003221255);
	request.solution.epochMjdTdb = 61000.0;
	request.solution.nonGrav.setModel(NgModel::Yarkovsky);
	request.solution.nonGrav.setEnabled(true);
	request.solution.nonGrav.setFitParameterMask(NgParameterMask::A2);
	request.solution.nonGrav.setParams(0.0, -2e-12, 0.0);
	request.startMjd = request.endMjd = 61100.0;
	request.uncertainty = true;
	Eigen::Matrix<double, 7, 7> factor = Eigen::Matrix<double, 7, 7>::Zero();
	factor.diagonal() << 1e-8, 2e-8, 3e-8, 2e-10, 3e-10, 4e-10, 1e-12;
	factor(2, 6) = 2e-9;
	request.solution.covariance = factor * factor.transpose();
	const auto nominal = oc::computeEphemeris(request).rows.at(0);
	ASSERT_TRUE(nominal.uncertainty);
	Eigen::Matrix<double, 2, 7> jacobian;
	Eigen::Matrix<double, 1, 7> rangeRateJacobian;
	request.uncertainty = false;
	for (int i = 0; i < 7; ++i) {
		auto plus = request, minus = request;
		const double step = i < 3 ? 1e-6 : i < 6 ? 1e-8 : 1e-10;
		if (i == 6) {
			plus.solution.nonGrav.setFitParameterValue(0, request.solution.nonGrav.A2() + step);
			minus.solution.nonGrav.setFitParameterValue(0, request.solution.nonGrav.A2() - step);
		} else {
			Vec3D& p = i < 3 ? plus.solution.position : plus.solution.velocity;
			Vec3D& m = i < 3 ? minus.solution.position : minus.solution.velocity;
			double* pc[] = {&p.x, &p.y, &p.z};
			double* mc[] = {&m.x, &m.y, &m.z};
			*pc[i % 3] += step;
			*mc[i % 3] -= step;
		}
		const auto a = oc::computeEphemeris(plus).rows.at(0);
		const auto b = oc::computeEphemeris(minus).rows.at(0);
		jacobian(0, i) = std::remainder(a.longitude - b.longitude, 2.0 * Constants::PI) *
		                 std::cos(nominal.latitude) / (2.0 * step);
		jacobian(1, i) = (a.latitude - b.latitude) / (2.0 * step);
		rangeRateJacobian[i] = (a.distanceRateAuPerDay - b.distanceRateAuPerDay) / (2.0 * step);
	}
	const Eigen::Matrix2d expected = jacobian * request.solution.covariance * jacobian.transpose();
	constexpr double arcsec = 206264.80624709636;
	const double east = 3.0 * std::sqrt(expected(0, 0)) * arcsec;
	const double north = 3.0 * std::sqrt(expected(1, 1)) * arcsec;
	EXPECT_NEAR(nominal.uncertainty->longitude3SigmaArcsec, east, east * 0.002);
	EXPECT_NEAR(nominal.uncertainty->latitude3SigmaArcsec, north, north * 0.002);
	const double rangeRate = 3.0 * std::sqrt((rangeRateJacobian * request.solution.covariance *
	                                          rangeRateJacobian.transpose())(0, 0)) *
	                         (Constants::AU_IN_METERS / 1000.0) / 86400.0;
	EXPECT_NEAR(nominal.uncertainty->rangeRate3SigmaKmPerSecond, rangeRate, rangeRate * 0.002);
}

TEST_F(EphemerisTest, InvalidEphemerisCovarianceIsNotReportedAsZeroUncertainty)
{
	auto request = buRequest(oc::EphemerisCenter::Geocentric);
	request.startMjd = request.endMjd = kBuEpochMjdTdb;
	request.uncertainty = true;
	request.solution.covariance = Eigen::Matrix<double, 6, 6>::Identity();
	request.solution.covariance(0, 0) = -1.0;
	const auto result = oc::computeEphemeris(request);
	ASSERT_EQ(result.rows.size(), 1u);
	EXPECT_FALSE(result.uncertaintyComputed);
	EXPECT_FALSE(result.rows.front().uncertainty);
}

TEST_F(EphemerisTest, HyperbolicLongArcLightTimeWithNonGravitationalParameters)
{
	for (const bool lightTime : {false, true}) {
		for (const bool uncertainty : {false, true}) {
			auto request = buRequest(oc::EphemerisCenter::Geocentric);
			request.solution.position = {-1.4514193884395477, 1.5432737371263698, -0.18204353174280821};
			request.solution.velocity = {-0.0065919582586045824, -0.018174932418503386, -0.015876564630493396};
			request.solution.epochMjdTdb = 58793.873949228786;
			request.solution.forces.useSb441 = m_sb441 && m_sb441->isAvailable();
			auto& ng = request.solution.nonGrav;
			ng.setModel(NgModel::InverseSquare);
			ng.setEnabled(true);
			ng.setParams(4.1793572990108531e-08, 3.9732593271940552e-08, 2.2863821588803364e-09);
			ng.setDT(-61.518218831046198);
			ng.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3 |
			                       NgParameterMask::DT);
			request.solution.covariance = Eigen::MatrixXd::Identity(10, 10) * 1e-14;
			request.startMjd = 58793;
			request.endMjd = 131842;
			request.stepDays = 1.0;
			request.lightTime = lightTime;
			request.uncertainty = uncertainty;
			const auto result = oc::computeEphemeris(request);
			ASSERT_EQ(result.rows.size(), 73050u);
			EXPECT_EQ(result.uncertaintyComputed, uncertainty);
			if (uncertainty)
				EXPECT_TRUE(result.rows.back().uncertainty);
		}
	}
}

TEST_F(EphemerisTest, AgreesWithHorizonsThroughA2023BuCloseApproach)
{
	// Horizons astrometric RA/Dec (deg) at CA -6 h, -1 h, +1 h, +6 h (UT).
	constexpr std::array<double, 4> kMjdUtc{59970.77010, 59970.97843, 59971.06177, 59971.27010};
	constexpr std::array<std::array<double, 2>, 4> kGeocentric{{{130.489709204, 34.495515931},
	                                                            {102.499210956, 6.604643309},
	                                                            {275.339088143, 2.371867587},
	                                                            {250.009722532, 29.796542177}}};
	constexpr std::array<std::array<double, 2>, 4> kW68{{{131.829356529, 34.357187700},
	                                                     {114.955710028, 15.869356449},
	                                                     {272.078074652, 8.069984719},
	                                                     {253.629706269, 31.122175263}}};

	for (const bool topocentric : {false, true}) {
		const auto& expected = topocentric ? kW68 : kGeocentric;
		for (std::size_t i = 0; i < kMjdUtc.size(); ++i) {
			auto request = topocentric ? buRequest(oc::EphemerisCenter::Station, "W68")
			                           : buRequest(oc::EphemerisCenter::Geocentric);
			request.startMjd = request.endMjd = kMjdUtc[i];
			const auto result = oc::computeEphemeris(request);
			ASSERT_EQ(result.rows.size(), 1u);
			const auto& row = result.rows.front();
			const double error = separationArcsec(row.longitude, row.latitude, expected[i][0] * Constants::DEG_TO_RAD,
			                                      expected[i][1] * Constants::DEG_TO_RAD);
			EXPECT_LT(error, 0.5) << (topocentric ? "W68" : "geocentric") << " epoch " << i;
			RecordProperty((topocentric ? "w68_" : "geocentric_") + std::to_string(i) + "_arcsec",
			               std::to_string(error));
		}
	}
}

TEST_F(EphemerisTest, TimeScalesDescribeTheSameInstants)
{
	auto utc = buRequest(oc::EphemerisCenter::Station, "W68");
	utc.startMjd = utc.endMjd = 59970.5;
	auto tdb = utc;
	tdb.scale = oc::time::Scale::TDB;
	tdb.startMjd = tdb.endMjd =
	    oc::time::fromMjdTt(oc::time::toMjdTt(59970.5, oc::time::Scale::UTC), oc::time::Scale::TDB);
	const auto a = oc::computeEphemeris(utc).rows.at(0);
	const auto b = oc::computeEphemeris(tdb).rows.at(0);
	auto tt = utc;
	tt.scale = oc::time::Scale::TT;
	tt.startMjd = tt.endMjd = oc::time::fromMjdTdb(tdb.startMjd, oc::time::Scale::TT);
	const auto c = oc::computeEphemeris(tt).rows.at(0);
	EXPECT_DOUBLE_EQ(a.mjdTdb, b.mjdTdb);
	EXPECT_DOUBLE_EQ(a.mjdTdb, c.mjdTdb);
	EXPECT_LT((a.position - b.position).magnitude() * Constants::AU_IN_METERS, 1e-6);
	EXPECT_LT((a.position - c.position).magnitude() * Constants::AU_IN_METERS, 1e-6);
	EXPECT_NEAR(a.mjdTt, b.mjdTt, 1.0e-11);
	EXPECT_LT(separationArcsec(a.longitude, a.latitude, b.longitude, b.latitude), 1.0e-3);
	EXPECT_DOUBLE_EQ(b.mjd, tdb.startMjd);
	EXPECT_NEAR(b.mjdUtc, 59970.5, 1.0e-11);
}

TEST_F(EphemerisTest, TdbSolutionEpochDoesNotShiftTheInitialState)
{
	auto request = buRequest(oc::EphemerisCenter::Heliocentric);
	request.scale = oc::time::Scale::TDB;
	request.plane = oc::EphemerisPlane::Ecliptic;
	request.startMjd = request.endMjd = kBuEpochMjdTdb;
	const auto row = oc::computeEphemeris(request).rows.at(0);
	EXPECT_DOUBLE_EQ(row.mjdTdb, kBuEpochMjdTdb);
	EXPECT_EQ(row.position.x, kBuPosition.x);
	EXPECT_EQ(row.position.y, kBuPosition.y);
	EXPECT_EQ(row.position.z, kBuPosition.z);
	EXPECT_EQ(row.velocity.x, kBuVelocity.x);
	EXPECT_EQ(row.velocity.y, kBuVelocity.y);
	EXPECT_EQ(row.velocity.z, kBuVelocity.z);
}

TEST_F(EphemerisTest, LightTimeIsAstrometricAndHeliocentricIsGeometric)
{
	auto request = buRequest(oc::EphemerisCenter::Geocentric);
	request.startMjd = 59960.0;
	request.endMjd = 59960.0;
	const auto astrometric = oc::computeEphemeris(request).rows.at(0);
	request.lightTime = false;
	const auto geometric = oc::computeEphemeris(request).rows.at(0);
	const double tauDays = geometric.distanceAu / Constants::AU_PER_DAY;
	const double shift = (astrometric.position - geometric.position).magnitude();
	// The target moves about 30 km/s along its orbit during the delay.
	EXPECT_GT(shift, 0.5 * 0.0175 * tauDays);
	EXPECT_LT(shift, 2.0 * 0.0175 * tauDays);
	EXPECT_TRUE(astrometric.visualMagnitude.has_value());
	EXPECT_TRUE(std::isfinite(astrometric.phase));

	request.center = oc::EphemerisCenter::Heliocentric;
	request.lightTime = true;
	const auto helio = oc::computeEphemeris(request);
	EXPECT_FALSE(helio.lightTimeApplied);
	EXPECT_NEAR(helio.rows.at(0).distanceAu, helio.rows.at(0).heliocentricAu, 1.0e-15);
	EXPECT_TRUE(std::isnan(helio.rows.at(0).elongation));
	EXPECT_FALSE(helio.rows.at(0).visualMagnitude.has_value());
}

TEST_F(EphemerisTest, StationsAreValidatedAndGeocentreMatchesCode500)
{
	auto request = buRequest(oc::EphemerisCenter::Station, "XXX");
	request.startMjd = request.endMjd = 59960.0;
	EXPECT_THROW((void)oc::computeEphemeris(request), std::invalid_argument);
	request.station = "250";  // Hubble: no fixed geodetic position
	EXPECT_THROW((void)oc::computeEphemeris(request), std::invalid_argument);

	request.station = "500";
	const auto station500 = oc::computeEphemeris(request).rows.at(0);
	request.center = oc::EphemerisCenter::Geocentric;
	const auto geocentre = oc::computeEphemeris(request).rows.at(0);
	EXPECT_LT(separationArcsec(station500.longitude, station500.latitude, geocentre.longitude, geocentre.latitude),
	          1.0e-9);
}

TEST_F(EphemerisTest, SolutionWithoutNBodyUsesTwoBodyMotionAndCancellationReturnsNoRows)
{
	auto request = buRequest(oc::EphemerisCenter::Geocentric);
	request.startMjd = 59960.0;
	request.endMjd = 59961.0;
	request.solution.forces.enabled = false;
	request.solution.forces.useSb441 = false;
	const auto twoBody = oc::computeEphemeris(request);
	EXPECT_FALSE(twoBody.nBody);
	EXPECT_EQ(twoBody.rows.size(), 2u);

	std::stop_source stop;
	stop.request_stop();
	const auto cancelled = oc::computeEphemeris(request, {.cancel = stop.get_token()});
	EXPECT_TRUE(cancelled.cancelled);
	EXPECT_TRUE(cancelled.rows.empty());
}

TEST_F(EphemerisTest, FineGridIsStrictlyIncreasingAndBounded)
{
	auto request = buRequest(oc::EphemerisCenter::Heliocentric);
	request.solution.forces.enabled = false;
	request.solution.forces.useSb441 = false;
	request.scale = oc::time::Scale::TT;
	request.startMjd = 59960.0;
	request.stepDays = 1e-7;
	request.endMjd = request.startMjd + 2 * request.stepDays;
	const auto result = oc::computeEphemeris(request);
	ASSERT_EQ(result.rows.size(), 3u);
	for (std::size_t i = 0; i < result.rows.size(); ++i) {
		EXPECT_LE(result.rows[i].mjd, request.endMjd);
		if (i)
			EXPECT_GT(result.rows[i].mjd, result.rows[i - 1].mjd);
	}
}

TEST_F(EphemerisTest, CancellationDuringNBodyPropagationReturnsNoPartialRows)
{
	class CancellingSource final : public IEphemerisSource {
	   public:
		explicit CancellingSource(IEphemerisSource& source) : source(source) {}

		EphemerisRecord queryJD(EphemerisPlanet planet, double jd) override { return queryMjd(planet, jd - 2400000.5); }

		EphemerisRecord queryMjd(EphemerisPlanet planet, double mjd) override
		{
			if (++queries == 30)
				stop.request_stop();
			return source.queryMjd(planet, mjd);
		}

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet planet, double jd, int days) override
		{
			return source.fetchWindow(planet, jd, days);
		}

		bool isAvailable() const override { return true; }

		std::string name() const override { return "CancellingSource"; }

		IEphemerisSource& source;
		std::stop_source stop;
		int queries = 0;
	} source(*m_de440);

	setDefaultEphemerisSource(&source);
	auto request = buRequest(oc::EphemerisCenter::Geocentric);
	request.startMjd = 59960.0;
	request.endMjd = 60060.0;
	const auto result = oc::computeEphemeris(request, {.cancel = source.stop.get_token()});
	setDefaultEphemerisSource(m_de440.get());
	EXPECT_TRUE(source.stop.stop_requested());
	EXPECT_LT(source.queries, 100);
	EXPECT_TRUE(result.cancelled);
	EXPECT_TRUE(result.rows.empty());
	request.startMjd = request.endMjd = 59960.0;
	EXPECT_EQ(oc::computeEphemeris(request).rows.size(), 1u);
}
