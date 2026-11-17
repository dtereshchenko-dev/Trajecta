#include <gtest/gtest.h>
#include "AdesParser.h"
#include "MPECParser.h"
#include "obs/Mpc80ColumnReader.h"
#include "fit/OrbitFitModel.h"
#include "fit/OrbitFitProblem.h"
#include "orbital_core/Residuals.h"
#include "orbital_core/Types.h"
#include "orbital_core/ObservationSet.h"
#include "orbital_core/DE440Source.h"
#include "fit/FitObservations.h"
#include "obs/RadarPrediction.h"
#include "earth_orientation/EarthPck.h"
#include "orbital_core/SolutionRefinement.h"
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <iterator>

namespace {
	const std::string delayDoppler = "99942         R2005 01 29.000000  19202850713  -    10251291   2380 251 JPLRS251";
	const std::string delayDopplerSigma =
	    "99942         r2005 01 29.000000C         4000           025        251 JPLRS251";
	const std::string doppler = "99942         R2005 01 27.979861               -    10084914   2380 251 JPLRS251";
	const std::string dopplerSigma = "99942         r2005 01 27.979861C                        0250       251 JPLRS251";
}  // namespace

TEST(RadarObservation, MpcApophisUnitsAndMissingDelay)
{
	MPECParser parser({delayDoppler, delayDopplerSigma, doppler, dopplerSigma});
	const auto records = parser.parse();
	ASSERT_EQ(parser.error(), MPECParser::Error::NoError);
	ASSERT_EQ(records.size(), 2u);
	ASSERT_TRUE(records[0].radar);
	const auto& r = *records[0].radar;
	ASSERT_TRUE(r.delaySeconds);
	ASSERT_TRUE(r.dopplerHz);
	EXPECT_DOUBLE_EQ(r.delaySeconds->value, 192.02850713);
	EXPECT_DOUBLE_EQ(r.delaySeconds->sigma, 4e-6);
	EXPECT_DOUBLE_EQ(r.dopplerHz->value, -102512.91);
	EXPECT_DOUBLE_EQ(r.dopplerHz->sigma, 0.25);
	EXPECT_DOUBLE_EQ(r.frequencyHz, 2380e6);
	EXPECT_EQ(r.transmitter, "251");
	EXPECT_EQ(records[0].code, "251");
	EXPECT_EQ(r.dimensions(), 2);
	ASSERT_TRUE(records[1].radar);
	EXPECT_FALSE(records[1].radar->delaySeconds);
	EXPECT_EQ(records[1].radar->dimensions(), 1);
	EXPECT_DOUBLE_EQ(records[1].radar->dopplerHz->value, -100849.14);
	EXPECT_DOUBLE_EQ(records[1].radar->dopplerHz->sigma, 0.25);
}

TEST(RadarObservation, MalformedRadarNeverSilentlyDrops)
{
	for (const auto& lines :
	     std::vector<std::vector<std::string>>{{doppler},
	                                           {dopplerSigma},
	                                           {doppler, delayDopplerSigma},
	                                           {doppler, std::string(dopplerSigma).replace(77, 3, "253")},
	                                           {doppler, std::string(dopplerSigma).replace(47, 15, 15, ' ')}}) {
		MPECParser parser(lines);
		EXPECT_TRUE(parser.parse().empty());
		EXPECT_NE(parser.error(), MPECParser::Error::NoError);
	}
}

TEST(RadarObservation, MpcRadarRecoversWholeUtcSecondsAndNormalizesMidnight)
{
	struct Epoch {
		const char* date;
		int mjd;
		double seconds;
	};

	for (const auto& epoch : {Epoch{"2005 10 30.570833", 53673, 49320}, Epoch{"2005 01 27.979861", 53397, 84660},
	                          Epoch{"2000 01 01.000000", 51544, 0}, Epoch{"2000 01 01.000012", 51544, 1},
	                          Epoch{"2000 01 01.000023", 51544, 2}, Epoch{"2000 01 01.999988", 51544, 86399},
	                          Epoch{"2000 01 01.999999", 51545, 0}, Epoch{"2000 01 31.999999", 51575, 0},
	                          Epoch{"2000 02 29.999999", 51604, 0}, Epoch{"2000 12 31.999999", 51910, 0}}) {
		for (const char kind : {'R', 'Q'}) {
			SCOPED_TRACE(epoch.date);
			SCOPED_TRACE(kind);
			auto first = doppler;
			auto second = dopplerSigma;
			first.replace(15, 17, epoch.date);
			second.replace(15, 17, epoch.date);
			first[14] = kind;
			second[14] = kind == 'R' ? 'r' : 'q';
			std::string error;
			const auto records = orbital_core::obs::Mpc80ColumnReader{}.read(first + '\n' + second + '\n', error);
			ASSERT_TRUE(error.empty()) << error;
			ASSERT_EQ(records.size(), 1u);
			const auto& record = records.front();
			ASSERT_TRUE(record.radar);
			EXPECT_TRUE(record.hasSplitUtc);
			EXPECT_EQ(record.utcMjdDay, epoch.mjd);
			EXPECT_DOUBLE_EQ(record.utcSecondsOfDay, epoch.seconds);
			EXPECT_DOUBLE_EQ(record.jd, 2400000.5 + (epoch.mjd + epoch.seconds / 86400.0));
			EXPECT_DOUBLE_EQ(record.timeAccuracyDays, 1.0 / 86400.0);
		}
	}
}

TEST(RadarObservation, MpcOpticalPreservesFractionalUtcSeconds)
{
	std::string error;
	const auto records = orbital_core::obs::Mpc80ColumnReader{}.read(
	    "01862          2005 10 30.57083302 45 33.68 +30 58 37.7                 19628690\n", error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(records.size(), 1u);
	EXPECT_FALSE(records[0].radar);
	EXPECT_TRUE(records[0].hasSplitUtc);
	EXPECT_EQ(records[0].utcMjdDay, 53673);
	EXPECT_NEAR(records[0].utcSecondsOfDay, 49319.9712, 1e-11);
	EXPECT_DOUBLE_EQ(records[0].timeAccuracyDays, 1e-6);
}

TEST(RadarObservation, AdesXmlAndPsvPreserveFractionalUtcSeconds)
{
	const std::string xml =
	    "<ades version=\"2022\"><radar><permID>99942</permID><trx>251</trx><rcv>251</rcv>"
	    "<obsTime>2005-10-30T13:41:59.9712Z</obsTime><com>1</com><frq>2380</frq>"
	    "<delay>192.02850713</delay><rmsDelay>4</rmsDelay></radar></ades>";
	const std::string psv =
	    "permID|trx|rcv|obsTime|com|frq|doppler|rmsDoppler\n"
	    "99942|251|251|2005-10-30T13:41:59.9712Z|1|2380|-10251.291|0.25\n";
	for (const auto& text : {xml, psv}) {
		std::string error;
		const auto records = orbital_core::ades::AdesReader{}.read(text, error);
		ASSERT_TRUE(error.empty()) << error;
		ASSERT_EQ(records.size(), 1u);
		ASSERT_TRUE(records[0].radar);
		EXPECT_TRUE(records[0].hasSplitUtc);
		EXPECT_EQ(records[0].utcMjdDay, 53673);
		EXPECT_DOUBLE_EQ(records[0].utcSecondsOfDay, 49319.9712);
		EXPECT_NEAR(records[0].timeAccuracyDays, 1e-4 / 86400.0, 1e-23);
	}
}

TEST(RadarObservation, AdesUnitsMatchMpc)
{
	orbital_core::ades::AdesReader reader;
	std::string error;
	const auto records = reader.read(
	    "<ades version=\"2022\"><radar><permID>99942</permID><trx>251</trx><rcv>251</rcv>"
	    "<obsTime>2005-01-29T00:00:00Z</obsTime><com>1</com><frq>2380</frq>"
	    "<delay>192.02850713</delay><rmsDelay>4</rmsDelay></radar></ades>",
	    error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(records.size(), 1u);
	ASSERT_TRUE(records[0].radar);
	EXPECT_DOUBLE_EQ(records[0].radar->delaySeconds->value, 192.02850713);
	EXPECT_DOUBLE_EQ(records[0].radar->delaySeconds->sigma, 4e-6);
	EXPECT_FALSE(records[0].radar->dopplerHz);
	const auto psv = reader.read(
	    "permID|trx|rcv|obsTime|com|frq|doppler|rmsDoppler\n"
	    "99942|251|251|2005-01-29T00:00:00Z|1|2380|-10251.291|0.25\n",
	    error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(psv.size(), 1u);
	EXPECT_DOUBLE_EQ(psv[0].radar->dopplerHz->value, -10251.291);
}

TEST(RadarObservation, CompactWhiteningAndOpticalRms)
{
	using namespace orbital_core;
	obs::PreparedObservation optical;
	optical.sigmaRa = optical.sigmaDec = Constants::ARCSEC_RAD;
	obs::PreparedObservation radar;
	radar.radar =
	    obs::RadarObservation{.dopplerHz = obs::RadarMeasurement{0, 0.25}, .frequencyHz = 2380e6, .transmitter = "251"};
	const auto orbit = make_orbit_from_states_vectors(StateVectors(Vec3D(1, 0.1, 0.02), Vec3D(0, 0.9, 0.1)), 2451545.0);
	OrbitOptimizerConfig config;
	config.useWeights = true;
	fit::OrbitFitModel model(orbit, {optical, radar, optical}, config);
	fit::OrbitFitProblem problem(model);
	EXPECT_EQ(problem.numResiduals(), 5);
	EXPECT_EQ(model.residualOffset(2), 3);
	std::vector<double> res{Constants::ARCSEC_RAD, 2 * Constants::ARCSEC_RAD, -0.5, 3 * Constants::ARCSEC_RAD,
	                        4 * Constants::ARCSEC_RAD};
	EXPECT_NEAR(model.opticalRms(res), std::sqrt(7.5), 1e-14);
	problem.weightResiduals(res, {});
	const std::vector<double> expected{1, 2, -2, 3, 4};
	for (std::size_t i = 0; i < res.size(); ++i)
		EXPECT_DOUBLE_EQ(res[i], expected[i]);
	EXPECT_DOUBLE_EQ(problem.sumSquares(res), 34);
	Eigen::MatrixXd jac = Eigen::MatrixXd::Ones(5, 6);
	problem.weightJacobian(jac, {});
	EXPECT_DOUBLE_EQ(jac(2, 0), 4);
	config.useWeights = false;
	EXPECT_THROW((fit::OrbitFitModel(orbit, {optical, radar}, config)), std::invalid_argument);
}

TEST(RadarObservation, RadarDoesNotChangeOpticalRms)
{
	Residual optical{};
	optical.dRA_arcsec = 3;
	optical.dDEC_arcsec = 4;
	Residual radar{};
	radar.radar.emplace();
	radar.dRA_arcsec = radar.dDEC_arcsec = std::numeric_limits<double>::quiet_NaN();
	Residuals residuals;
	residuals.add(optical);
	residuals.add(radar);
	EXPECT_DOUBLE_EQ(residuals.rms(), std::sqrt(12.5));
}

TEST(RadarObservation, SameStatePredictionsMatchIndependentReference)
{
	using namespace orbital_core;
	DE440Source eph;
	const auto directory = std::filesystem::path(__FILE__).parent_path() / "data";
	std::ifstream input(directory / "radar_cases.dat");
	std::ifstream expected(directory / "radar_expected.dat");
	ASSERT_TRUE(input && expected);
	int count;
	input >> count;
	ASSERT_EQ(count, 10);
	const auto readVector = [](std::istream& stream) {
		Vec3D v;
		stream >> v.x >> v.y >> v.z;
		return v;
	};
	const auto observer = [](const Vec3D& v, const char* code) {
		return obs::ObserverData{std::atan2(v.y, v.x) * Constants::RAD_TO_DEG,
		                         std::hypot(v.x, v.y) / Constants::EARTH_RADIUS_IN_AU,
		                         v.z / Constants::EARTH_RADIUS_IN_AU, code};
	};
	for (int i = 0; i < count; ++i) {
		SCOPED_TRACE(i);
		double epoch, frequency;
		obs::PreparedObservation record;
		input >> record.mjdTdb >> epoch >> frequency;
		record.mjdTdb = time::toMjdTdb(record.mjdTdb, time::Scale::TT);
		const auto position = readVector(input), velocity = readVector(input);
		record.observer = observer(readVector(input), "251");
		record.transmitter = observer(readVector(input), i % 2 ? "253" : "251");
		record.radar = obs::RadarObservation{.delaySeconds = obs::RadarMeasurement{1, 4e-6},
		                                     .dopplerHz = obs::RadarMeasurement{0, 0.1},
		                                     .frequencyHz = frequency,
		                                     .transmitter = record.transmitter.obsCode};
		int index;
		double delay, dopplerValue;
		expected >> index >> delay >> dopplerValue;
		readVector(expected);
		readVector(expected);
		ASSERT_TRUE(input && expected);
		const auto prediction = obs::predictRadar(record, position, velocity, eph);
		EXPECT_NEAR(prediction.delaySeconds, delay, 5e-10);
		EXPECT_NEAR(prediction.dopplerHz, dopplerValue, 0.001);
	}
}

TEST(RadarObservation, PckComparisonWithGrss)
{
	using namespace orbital_core;
	Vec3D position, velocity;
	if (!earthPckState(53398.0306038786, Vec3D(1, 0, 0), position, velocity))
		GTEST_SKIP() << "Earth PCK does not cover the reference epoch";
	DE440Source eph;
	const auto eqToEcl = [](const Vec3D& v) {
		const double c = std::cos(Constants::EARTH_INCL), s = std::sin(Constants::EARTH_INCL);
		return Vec3D(v.x, c * v.y + s * v.z, -s * v.y + c * v.z);
	};
	const Vec3D target(-0.509234620995513, 0.8741860351265314, 0.31192062720034325);
	const Vec3D targetVel(-0.014261107919596566, -0.005301452195241065, -0.00233848510505264);
	const Vec3D sun(0.004309465807282884, -0.00005377459807188601, -0.00013797382606907447);
	const Vec3D sunVel(3.11018104435471e-07, 6.638445944852413e-06, 2.8027487717581985e-06);
	const auto observer = [](const Vec3D& v, const char* code) {
		return obs::ObserverData{std::atan2(v.y, v.x) * Constants::RAD_TO_DEG,
		                         std::hypot(v.x, v.y) / Constants::EARTH_RADIUS_IN_AU,
		                         v.z / Constants::EARTH_RADIUS_IN_AU, code};
	};
	const Vec3D rx(1.5979452128291048e-05, -3.7198015198121182e-05, 1.3334e-05);
	const Vec3D tx(-1.5732979805588661e-05, -3.1025487683466637e-05, 2.4579600000000001e-05);
	ASSERT_TRUE(earthPckState(53398.0306038786, rx, position, velocity));
	const Vec3D grssPosition(1.2974700944412554e-05, 3.8352280704812e-05, 1.3326307858663178e-05);
	const Vec3D grssVelocity(-0.0002416309247645805, 8.170523254304e-05, 1.1308265004342927e-07);
	EXPECT_LT((position - grssPosition).magnitude(), 1e-13);
	EXPECT_LT((velocity - grssVelocity).magnitude(), 1e-13);

	struct Case {
		Vec3D station;
		double delay;
		double doppler;
	};

	for (const auto& item :
	     {Case{rx, 188.43312965399295, -103075.7351165945}, Case{tx, 188.4394916266295, -100518.14497411001}}) {
		obs::PreparedObservation record;
		record.mjdTdb = 53398.0306038786;
		record.observer = observer(rx, "251");
		record.transmitter = observer(item.station, "253");
		record.radar = obs::RadarObservation{.delaySeconds = obs::RadarMeasurement{1, 4e-6},
		                                     .dopplerHz = obs::RadarMeasurement{0, 0.1},
		                                     .frequencyHz = 2380e6,
		                                     .transmitter = record.transmitter.obsCode};
		const auto value = obs::predictRadar(record, eqToEcl(target - sun), eqToEcl(targetVel - sunVel), eph);
		const double transmit = record.mjdTdb - value.delaySeconds / Constants::SECONDS_PER_DAY;
		const double clock = Calculations::tdbMinusTtSeconds(time::fromMjdTdb(transmit, time::Scale::TT)) -
		                     Calculations::tdbMinusTtSeconds(time::fromMjdTdb(record.mjdTdb, time::Scale::TT));
		// GRSS and Trajecta use different atmospheric and relativistic corrections.
		EXPECT_NEAR(value.delaySeconds - clock, item.delay, 1e-7);
		EXPECT_NEAR(value.dopplerHz, item.doppler, 0.005);
	}
}

TEST(RadarObservation, SurfaceReturnWithoutRadiusUsesOneKilometre)
{
	using namespace orbital_core;

	struct RestoreEphemeris {
		IEphemerisSource* previous{getDefaultEphemerisSource()};

		~RestoreEphemeris() { setDefaultEphemerisSource(previous); }
	} restore;

	DE440Source ephemeris;
	ASSERT_TRUE(ephemeris.isAvailable());
	setDefaultEphemerisSource(&ephemeris);
	const auto path =
	    std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data/examples/NEO/Apollo1862.txt";
	const auto observations = loadObservations(path);
	std::vector<ObsNode> radar;
	std::ranges::copy_if(observations.all, std::back_inserter(radar),
	                     [](const auto& o) { return o.radar.has_value(); });
	ASSERT_EQ(radar.size(), 17u);
	EXPECT_TRUE(std::ranges::none_of(radar, [](const auto& o) { return o.radar->surfaceRadiusKm.has_value(); }));
	const auto prepared = fit::prepareFitObservations(std::span<const ObsNode>(radar), OrbitOptimizerConfig{});
	int surface = 0;
	int centre = 0;
	for (const auto& record : prepared) {
		ASSERT_TRUE(record.radar);
		if (record.radar->returnType == obs::RadarReturn::Surface) {
			++surface;
			ASSERT_TRUE(record.radar->surfaceRadiusKm);
			EXPECT_DOUBLE_EQ(*record.radar->surfaceRadiusKm, 1.0);
		} else {
			++centre;
			EXPECT_FALSE(record.radar->surfaceRadiusKm);
		}
	}
	EXPECT_EQ(surface, 12);
	EXPECT_EQ(centre, 5);
	EXPECT_TRUE(std::ranges::none_of(radar, [](const auto& o) { return o.radar->surfaceRadiusKm.has_value(); }));

	const auto found = std::ranges::find_if(prepared, [](const auto& record) {
		return record.radar->returnType == obs::RadarReturn::Surface && record.radar->delaySeconds;
	});
	ASSERT_NE(found, prepared.end());
	auto record = *found;
	const auto up = record.observerOffsetEclAu / record.observerOffsetEclAu.magnitude();
	const auto position = record.earthR + up * 0.05;
	const auto velocity = record.earthV + up * 0.001;
	const double surfaceDelay = obs::predictRadar(record, position, velocity, ephemeris).delaySeconds;
	record.radar->returnType = obs::RadarReturn::CenterOfMass;
	record.radar->surfaceRadiusKm.reset();
	const double centreDelay = obs::predictRadar(record, position, velocity, ephemeris).delaySeconds;
	EXPECT_NEAR(centreDelay - surfaceDelay, 2.0e3 / 299792458.0, 1e-9);
}
