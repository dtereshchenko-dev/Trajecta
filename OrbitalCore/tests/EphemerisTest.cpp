#include <gtest/gtest.h>

#include "orbital_core/Constants.h"
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
				if (rows[i].RA < 0.0) rows[i].RA += 2.0 * Constants::PI;
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
	for (const auto center : {oc::EphemerisCenter::Geocentric, oc::EphemerisCenter::Station}) {
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
