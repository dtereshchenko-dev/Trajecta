#include <gtest/gtest.h>

#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "orbital_core/Exporter.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/ObservationSet.h"
#include "orbital_core/StateVectors.h"
#include "fit/FitObservations.h"
#include "obs/AstrometricResidual.h"
#include "obs/RwoReader.h"
#include "obs/ObservatoryCatalog.h"

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace orbital_core;

namespace {
	const std::string esaVersionThree =
	    "version =   3\nerrmod  = 'vfcc17'\nEND_OF_HEADER\n"
	    " 99942        O C   2004 03 15.10789       1.000E-05 04 06 08.080  1.435E-01    0.612 F   -0.247    0.251 +16 "
	    "55 "
	    "04.60  1.000E-01    0.612 F    0.140   -0.070                       o 691      0.43 1 0\n"
	    " 99942        S S   2020 12 18.97667       1.000E-05 11 30 15.530  1.474E-01    1.000 F    0.182   -0.700 -10 "
	    "46 "
	    "20.20  1.000E-01    1.000 F    0.135    0.887 19   R  1.00   0.01   L C51      1.13 1 1\n"
	    " 99942        S s   2020 12 18.97667 1      -5634.173400000000      -2466.265700000000      "
	    "-3038.392400000000 "
	    "C51\n"
	    " 99942        R c   2005 01 29 00:00:00  28784349.07929   0.59958   0.59958 F     0.00000     0.01857    -1   "
	    "-1  "
	    "       0.05  1\n"
	    " 99942        V c   2005 01 27 23:31:00    548781.78980   1.36040   1.36040 F     0.00000     0.22766    -1   "
	    "-1  "
	    "       0.17  1\n"
	    " 99942        R c   2021 03 11 08:20:00  17057650.31341   0.14990   0.14990 F     0.00000    -0.18525    -14  "
	    "-14 "
	    "       1.95  1\n";

	struct RwoRow {
		double sigmaRa, sigmaDec, biasRa, biasDec, residualRa, residualDec;
		double chi;
		char forceRa, forceDec;
		std::string magnitude;
		double magSigma{0.0};
		std::optional<double> magResidual;
		int astSelected{0};
		int magSelected{0};
	};

	struct RwoFile {
		double rms = -1.0;
		std::vector<RwoRow> rows;
	};

	RwoFile readRwo(const std::filesystem::path& path)
	{
		std::ifstream file(path);
		if (!file)
			throw std::runtime_error("Missing RWO file");
		RwoFile result;
		for (std::string line; std::getline(file, line);) {
			if (line.starts_with("RMSast"))
				result.rms = std::stod(line.substr(line.find('=') + 1));
			if (!line.starts_with(" 2020VV"))
				continue;
			const auto start = line.find("2020 ");
			if (start == std::string::npos)
				throw std::runtime_error("Missing RWO date");
			std::istringstream row(line.substr(start));
			double ignored;
			for (int i = 0; i < 8; ++i)
				row >> ignored;
			RwoRow parsed{};
			row >> parsed.sigmaRa >> parsed.forceRa >> parsed.biasRa >> parsed.residualRa;
			for (int i = 0; i < 4; ++i)
				row >> ignored;
			row >> parsed.sigmaDec >> parsed.forceDec >> parsed.biasDec >> parsed.residualDec;
			if (!row)
				throw std::runtime_error("Malformed RWO astrometry");
			std::istringstream tail(line);
			std::vector<std::string> tokens;
			for (std::string token; tail >> token;)
				tokens.push_back(token);
			parsed.chi = std::stod(tokens.at(tokens.size() - 3));
			parsed.astSelected = std::stoi(tokens.at(tokens.size() - 2));
			parsed.magSelected = std::stoi(tokens.back());
			parsed.magnitude = line.substr(156, 6);
			if (line.substr(163, 5).find_first_not_of(' ') != std::string::npos)
				parsed.magSigma = std::stod(line.substr(163, 5));
			if (line.substr(169, 6).find_first_not_of(' ') != std::string::npos)
				parsed.magResidual = std::stod(line.substr(169, 6));
			result.rows.push_back(parsed);
		}
		return result;
	}

	class ExporterRwoTest : public ::testing::Test {
	   protected:
		FitResult fit;
		std::vector<ObsNode> observations;
		std::filesystem::path root;
		bool ownsRoot = false;

		void SetUp() override
		{
			fit.orbit =
			    make_orbit_from_states_vectors(StateVectors(Vec3D(1.0, 0.1, 0.02), Vec3D(0.01, 0.9, 0.1)), 2459165.5);
			ASSERT_TRUE(fit.orbit);
			fit.fitEpochOrbit = fit.orbit;
			const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
			for (int attempt = 0; attempt < 100 && !ownsRoot; ++attempt) {
				root = std::filesystem::temp_directory_path() /
				       ("orbital_rwo_test_" + std::to_string(stamp) + "-" + std::to_string(attempt));
				ownsRoot = std::filesystem::create_directory(root);
			}
			ASSERT_TRUE(ownsRoot);
		}

		void TearDown() override
		{
			if (ownsRoot) {
				std::error_code error;
				std::filesystem::remove_all(root, error);
			}
		}

		Residual addObservation(const char* station, double sigma)
		{
			ObsNode observation;
			observation.JD = 2459165.5 + observations.size() * 0.001;
			observation.RA = 1.0;
			observation.DEC = 0.2;
			observation.obsCode = station;
			observation.note2 = 'C';
			observation.catalogCode = 'V';
			observation.hasObserverSigma = true;
			observation.rmsRaArcsec = sigma;
			observation.rmsDecArcsec = sigma;
			observations.push_back(observation);

			Residual row{};
			row.jd = observation.epochJdTdb();
			row.raObs = observation.RA;
			row.decObs = observation.DEC;
			row.obsCode = station;
			row.dRA_arcsec = sigma;
			row.dDEC_arcsec = -2.0 * sigma;
			row.total_arcsec = std::hypot(row.dRA_arcsec, row.dDEC_arcsec);
			obs::PreparedObservation prepared;
			prepared.sigmaRa = prepared.sigmaDec = sigma * Constants::ARCSEC_RAD;
			row.astrometricWeight = obs::residualWeightForReport(prepared, row.dRA_arcsec * Constants::ARCSEC_RAD,
			                                                     row.dDEC_arcsec * Constants::ARCSEC_RAD, true);
			return row;
		}

		RwoFile save()
		{
			ExportContext context;
			context.objectName = "2020VV";
			context.fitResult = &fit;
			context.observations = observations;
			context.rootDirectory = root;
			const auto report = Exporter::save(context, {ExportProduct::Residuals});
			return readRwo(report.directory / "weighted_residuals_rwo.txt");
		}
	};
}  // namespace

TEST_F(ExporterRwoTest, PreservesAllElevenEsaAssignmentsInsteadOfStationRules)
{
	const std::array<double, 11> sigmas{0.02, 0.02, 0.02, 0.08, 0.08, 0.08, 0.08, 0.07, 0.07, 0.08, 0.20};
	for (std::size_t i = 0; i < sigmas.size(); ++i)
		fit.residuals.add(addObservation((i < 3 || (i >= 7 && i < 10)) ? "309" : "G37", sigmas[i]));
	const auto result = save();
	ASSERT_EQ(result.rows.size(), sigmas.size());
	for (std::size_t i = 0; i < sigmas.size(); ++i) {
		EXPECT_DOUBLE_EQ(result.rows[i].sigmaRa, sigmas[i]);
		EXPECT_DOUBLE_EQ(result.rows[i].sigmaDec, sigmas[i]);
		EXPECT_EQ(result.rows[i].forceRa, 'F');
		EXPECT_NEAR(result.rows[i].chi, std::sqrt(5.0), 0.005);
	}
	EXPECT_NEAR(result.rms, std::sqrt(5.0 / 2.0), 5e-6);
}

TEST_F(ExporterRwoTest, RovingVvRoundTripPreservesPositionWithoutSatelliteSubstitution)
{
	fit.residuals.add(addObservation("247", 0.8));
	observations.back().note2 = 'V';
	observations.back().rovingObserver = obs::RovingObserver{15.123456789, -28.765432109, 1520.345};
	ExportContext context;
	context.objectName = "2020VV";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	const auto report = Exporter::save(context, {ExportProduct::Residuals});
	std::ifstream file(report.directory / "weighted_residuals_rwo.txt");
	const std::string text{std::istreambuf_iterator<char>(file), {}};
	EXPECT_NE(text.find(" O V "), std::string::npos);
	EXPECT_NE(text.find(" O v "), std::string::npos);
	EXPECT_EQ(text.find(" S S "), std::string::npos);
	EXPECT_EQ(text.find(" S s "), std::string::npos);
	std::string error;
	const auto records = obs::RwoReader{}.read(text, error);
	ASSERT_EQ(records.size(), 1u) << error;
	ASSERT_TRUE(records[0].rovingObserver);
	EXPECT_EQ(records[0].rovingObserver, observations[0].rovingObserver);
	EXPECT_FALSE(records[0].isSatellite);
	EXPECT_EQ(records[0].obsCode, "247");

	std::istringstream lines(text);
	std::string native, line;
	while (std::getline(lines, line))
		if (!line.starts_with("!"))
			native += line + '\n';
	const auto plain = obs::RwoReader{}.read(native, error);
	ASSERT_EQ(plain.size(), 1u) << error;
	ASSERT_TRUE(plain[0].rovingObserver);
	EXPECT_DOUBLE_EQ(plain[0].rovingObserver->longitudeDeg, 15.123457);
	EXPECT_DOUBLE_EQ(plain[0].rovingObserver->latitudeDeg, -28.765432);
	EXPECT_DOUBLE_EQ(plain[0].rovingObserver->heightMeters, 1520.3);
	EXPECT_EQ(plain[0].importedFit->selected, records[0].importedFit->selected);

	const auto continuation = native.find(" O v ");
	ASSERT_NE(continuation, std::string::npos);
	std::string mismatched = native;
	const auto station = mismatched.find("247", continuation);
	ASSERT_NE(station, std::string::npos);
	mismatched.replace(station, 3, "270");
	EXPECT_TRUE(obs::RwoReader{}.read(mismatched, error).empty());
	EXPECT_FALSE(error.empty());
	EXPECT_TRUE(obs::RwoReader{}.read(native.substr(0, native.rfind('\n', native.size() - 2) + 1), error).empty());
	EXPECT_FALSE(error.empty());
	std::string wrongPosition = text;
	const auto longitude = wrongPosition.find("15.123457");
	ASSERT_NE(longitude, std::string::npos);
	wrongPosition.replace(longitude, 9, "16.123457");
	EXPECT_TRUE(obs::RwoReader{}.read(wrongPosition, error).empty());
	EXPECT_FALSE(error.empty());
}

TEST_F(ExporterRwoTest, RovingExportRequiresTerrestrialCoordinates)
{
	fit.residuals.add(addObservation("247", 0.8));
	observations.back().isSatellite = true;
	observations.back().note2 = 'V';
	EXPECT_THROW(save(), std::invalid_argument);
}

TEST_F(ExporterRwoTest, WesternRovingLongitudePreservesPrecisionAndOriginalSign)
{
	fit.residuals.add(addObservation("247", 0.8));
	observations.back().note2 = 'V';
	observations.back().rovingObserver = obs::RovingObserver{-122.123456789, 48.6233, 506.0};
	ExportContext context;
	context.objectName = "2020VV";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	const auto report = Exporter::save(context, {ExportProduct::Residuals});
	std::ifstream file(report.directory / "weighted_residuals_rwo.txt");
	const std::string text{std::istreambuf_iterator<char>(file), {}};
	EXPECT_NE(text.find("237.876543"), std::string::npos);
	std::string error;
	const auto records = obs::RwoReader{}.read(text, error);
	ASSERT_EQ(records.size(), 1u) << error;
	EXPECT_EQ(records[0].rovingObserver, observations[0].rovingObserver);
	std::istringstream lines(text);
	std::string native, line;
	while (std::getline(lines, line))
		if (!line.starts_with("!"))
			native += line + '\n';
	const auto plain = obs::RwoReader{}.read(native, error);
	ASSERT_EQ(plain.size(), 1u) << error;
	ASSERT_TRUE(plain[0].rovingObserver);
	EXPECT_DOUBLE_EQ(plain[0].rovingObserver->longitudeDeg, 237.876543);
}

TEST_F(ExporterRwoTest, NativeVersionTwoRoundTripsRadarAndObserverGeometry)
{
	fit.residuals.add(addObservation("C51", 0.12));
	observations.back().isSatellite = true;
	observations.back().note2 = 'E';
	const Vec3D position(1200, -2500, 3100);
	observations.back().satelliteOffsetEclAu = Convertors::ecvToEcl(position / 149597870.7, Constants::EARTH_INCL);
	ObsNode radar;
	radar.hasSplitUtc = true;
	radar.utcMjdDay = 59165;
	radar.utcSecondsOfDay = 49319.9712;
	radar.JD = 2400000.5 + (radar.utcMjdDay + radar.utcSecondsOfDay / 86400.0);
	radar.obsCode = "251";
	radar.radar.emplace();
	radar.radar->transmitter = "253";
	radar.radar->frequencyHz = 8.56e9;
	radar.radar->delaySeconds = obs::RadarMeasurement{123.456, 2e-6};
	radar.radar->dopplerHz = obs::RadarMeasurement{-73.25, 0.05};
	observations.push_back(radar);
	Residual residual{};
	residual.radar.emplace();
	residual.radar->delaySeconds = Residual::RadarChannel{123.456, 123.456001, -1e-6, 2e-6};
	residual.radar->dopplerHz = Residual::RadarChannel{-73.25, -73.3, 0.05, 0.05};
	residual.radar->chiSquared = 1.25;
	fit.residuals.add(residual);
	fit.keptObsIndices = {0, 1};
	ExportContext context;
	context.objectName = "99942";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	const auto report = Exporter::save(context, {ExportProduct::Residuals});
	ASSERT_EQ(report.files.size(), 1u);
	EXPECT_FALSE(std::filesystem::exists(report.directory / "radar_residuals.csv"));
	std::ifstream input(report.files[0]);
	const std::string text(std::istreambuf_iterator<char>(input), {});
	EXPECT_TRUE(text.starts_with("version =   2"));
	EXPECT_NE(text.find(" S S "), std::string::npos);
	EXPECT_EQ(text.find(" S E "), std::string::npos);
	EXPECT_NE(text.find("Radar range/range rate (km or km/d)"), std::string::npos);
	std::istringstream headerLines(text);
	int headers = 0;
	for (std::string line; std::getline(headerLines, line);) {
		if (!line.starts_with("! Design"))
			continue;
		EXPECT_EQ(line.find("K T N"), 11u);
		EXPECT_EQ(line.find("YYYY"), 17u);
		if (++headers == 1) {
			EXPECT_EQ(line.find("HH MM"), 50u);
			EXPECT_EQ(line.find("sDD"), 103u);
		} else {
			EXPECT_EQ(line.find("TRX"), 99u);
			EXPECT_EQ(line.find("RCX"), 103u);
		}
	}
	EXPECT_EQ(headers, 2);
	obs::RwoReader reader;
	std::string error;
	const auto parsed = reader.read(text, error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(parsed.size(), 3u);
	EXPECT_TRUE(parsed[0].isSatellite);
	EXPECT_EQ(parsed[0].technique, obs::Technique::Occultation);
	EXPECT_NEAR(parsed[0].offsetX * 149597870.7, position.x, 1e-7);
	EXPECT_NEAR(parsed[0].offsetY * 149597870.7, position.y, 1e-7);
	EXPECT_NEAR(parsed[0].offsetZ * 149597870.7, position.z, 1e-7);
	ASSERT_TRUE(parsed[1].radar && parsed[2].radar);
	for (const auto i : {1, 2}) {
		EXPECT_TRUE(parsed[i].hasSplitUtc);
		EXPECT_EQ(parsed[i].utcMjdDay, radar.utcMjdDay);
		EXPECT_DOUBLE_EQ(parsed[i].utcSecondsOfDay, radar.utcSecondsOfDay);
		EXPECT_DOUBLE_EQ(parsed[i].jd, radar.JD);
	}
	EXPECT_EQ(parsed[1].radar->transmitter, "253");
	EXPECT_EQ(parsed[2].obsCode, "251");
	EXPECT_DOUBLE_EQ(parsed[1].radar->delaySeconds->value, 123.456);
	EXPECT_DOUBLE_EQ(parsed[2].radar->dopplerHz->value, -73.25);
	EXPECT_DOUBLE_EQ(parsed[2].radar->frequencyHz, 8.56e9);
	EXPECT_FALSE(parsed[2].radar->frequencyIsReference);
	std::string native;
	std::istringstream lines(text);
	for (std::string line; std::getline(lines, line);)
		if (!line.starts_with("! OC_RADAR"))
			native += line + '\n';
	const auto nativeParsed = reader.read(native, error);
	ASSERT_TRUE(error.empty()) << error;
	EXPECT_TRUE(nativeParsed[2].radar->frequencyIsReference);
	EXPECT_NEAR(nativeParsed[2].radar->dopplerHz->value * 8.56e9, -73.25, 5e-6);
	auto broken = text;
	const auto pos = broken.find(" S s ");
	ASSERT_NE(pos, std::string::npos);
	broken[pos + 3] = 'v';
	EXPECT_TRUE(reader.read(broken, error).empty());
	EXPECT_FALSE(error.empty());
}

TEST_F(ExporterRwoTest, GaiaAstrometryAndPhotocenterSurviveExportAndReload)
{
	struct RestoreEphemeris {
		IEphemerisSource* previous{getDefaultEphemerisSource()};

		~RestoreEphemeris() { setDefaultEphemerisSource(previous); }
	} restore;

	DE440Source ephemeris;
	ASSERT_TRUE(ephemeris.isAvailable());
	setDefaultEphemerisSource(&ephemeris);
	for (int i = 0; i < 3; ++i) {
		fit.residuals.add(addObservation("258", 0.12));
		auto& observation = observations.back();
		observation.isSatellite = true;
		observation.astrometricConvention = obs::AstrometricConvention::GaiaFpr;
		observation.satelliteOffsetEclAu = Vec3D(0.001, 0.002, 0.003);
	}
	fit.keptObsIndices = {0, 1, 2};
	ExportContext context;
	context.objectName = "28751";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	context.gaiaPhotocenterRadiusKm = 0.75;
	const auto report = Exporter::save(context, {ExportProduct::Residuals});
	ASSERT_EQ(report.files.size(), 1u);
	const auto loaded = loadObservations(report.files.front());
	ASSERT_EQ(loaded.all.size(), 3u);
	ASSERT_TRUE(loaded.gaiaPhotocenterRadiusKm);
	EXPECT_DOUBLE_EQ(*loaded.gaiaPhotocenterRadiusKm, 0.75);
	for (const auto& observation : loaded.all)
		EXPECT_EQ(observation.astrometricConvention, obs::AstrometricConvention::GaiaFpr);
}

TEST_F(ExporterRwoTest, ReadsEsaVersionThreeWithoutLosingRadarOrObserverPositions)
{
	obs::RwoReader reader;
	std::string error;
	const auto parsed = reader.read(esaVersionThree, error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(parsed.size(), 5u);
	EXPECT_EQ(parsed[0].objectName, "99942");
	EXPECT_EQ(parsed[0].obsCode, "691");
	EXPECT_NEAR(parsed[0].raRad, Convertors::raToRadians(4, 6, 8.080), 1e-15);
	EXPECT_NEAR(parsed[0].rmsRaArcsec, 0.612, 1e-15);
	ASSERT_TRUE(parsed[0].importedFit);
	EXPECT_TRUE(parsed[0].importedFit->selected);
	EXPECT_DOUBLE_EQ(parsed[0].importedFit->bias1, -0.247);
	ASSERT_TRUE(parsed[1].isSatellite);
	EXPECT_EQ(parsed[1].offsetUnits, obs::OffsetUnits::Km);
	EXPECT_DOUBLE_EQ(parsed[1].offsetX, -5634.1734);
	EXPECT_DOUBLE_EQ(parsed[1].offsetY, -2466.2657);
	EXPECT_DOUBLE_EQ(parsed[1].offsetZ, -3038.3924);
	EXPECT_EQ(parsed[1].obsCode, "C51");
	ASSERT_TRUE(parsed[2].radar && parsed[3].radar && parsed[4].radar);
	EXPECT_EQ(parsed[2].radar->transmitter, "251");
	EXPECT_EQ(parsed[2].obsCode, "251");
	EXPECT_EQ(parsed[4].radar->transmitter, "253");
	EXPECT_EQ(parsed[4].obsCode, "253");
	EXPECT_TRUE(parsed[3].radar->frequencyIsReference);
	constexpr double cKmSec = Constants::AU_PER_DAY * Constants::AU_IN_METERS / 1000.0 / 86400.0;
	EXPECT_NEAR(parsed[2].radar->delaySeconds->value * cKmSec / 2, 28784349.07929, 1e-8);
	EXPECT_NEAR(parsed[3].radar->dopplerHz->value * -cKmSec * 86400 / 2, 548781.78980, 1e-9);
	EXPECT_NEAR(parsed[3].radar->dopplerHz->sigma * cKmSec * 86400 / 2, 1.36040, 1e-12);
}

TEST_F(ExporterRwoTest, RejectsUnknownRadarAliasesAndMismatchedRwoVersions)
{
	obs::RwoReader reader;
	std::string error;
	auto unknownStation = esaVersionThree;
	const auto station = unknownStation.find("-14");
	ASSERT_NE(station, std::string::npos);
	unknownStation.replace(station, 3, "-99");
	EXPECT_TRUE(reader.read(unknownStation, error).empty());
	EXPECT_NE(error.find("-99"), std::string::npos);
	EXPECT_NE(error.find("RADCODE.dat"), std::string::npos);
	for (char version : {'2', '4'}) {
		auto invalid = esaVersionThree;
		invalid[12] = version;
		EXPECT_TRUE(reader.read(invalid, error).empty());
		EXPECT_FALSE(error.empty());
	}
}

TEST_F(ExporterRwoTest, RadarAliasesResolveToOneCanonicalCatalogRecord)
{
	const auto& stations = obs::radarObservatoryCatalog();
	std::ifstream source(getExecutableDir() / "RADCODE.dat");
	ASSERT_TRUE(source);
	const auto parsed = obs::readRadarCatalog(source);
	for (const auto& [alias, canonical] : parsed.aliases) {
		SCOPED_TRACE(alias);
		EXPECT_EQ(&obs::radarObservatory(alias), &stations.at(canonical));
		EXPECT_FALSE(stations.contains(alias));
	}
	EXPECT_EQ(&obs::radarObservatory("251"), &stations.at("251"));
	EXPECT_THROW((void)obs::radarObservatory("-999"), std::invalid_argument);
	EXPECT_DOUBLE_EQ(obs::observatoryCatalog().at("273").cosPhi, 0.0);
	EXPECT_GT(obs::radarObservatory("273").cosPhi, 0.0);
}

TEST_F(ExporterRwoTest, VersionedRadarCatalogDerivesIdentifiersFromStationRows)
{
	std::istringstream input(
	    "! RADCODE format=2\r\n"
	    "Q01 -101 Z01 123.5 123456.789 -234567.891 Test station\r\n"
	    "Q02 - Q02 42 200000 100000 No JPL identifier\n"
	    "Q03 - - 43 300000 -200000 No alternative identifiers\n"
	    "Q04 -102 - 44 200000 100000 No MPC identifier\n");
	const auto parsed = obs::readRadarCatalog(input);
	ASSERT_EQ(parsed.stations.size(), 4u);
	ASSERT_EQ(parsed.aliases.size(), 3u);
	EXPECT_EQ(parsed.aliases.at("-101"), "Q01");
	EXPECT_EQ(parsed.aliases.at("Z01"), "Q01");
	EXPECT_EQ(parsed.aliases.at("-102"), "Q04");
	EXPECT_FALSE(parsed.aliases.contains("-"));
	const auto& record = parsed.stations.at("Q01");
	EXPECT_EQ(record.name, "Test station");
	EXPECT_DOUBLE_EQ(record.longitude, 123.5);
	EXPECT_DOUBLE_EQ(record.cosPhi, 123456.789 * (1e-10 / Constants::EARTH_RADIUS_IN_AU));
	EXPECT_DOUBLE_EQ(record.sinPhi, -234567.891 * (1e-10 / Constants::EARTH_RADIUS_IN_AU));
}

TEST_F(ExporterRwoTest, VersionedRadarCatalogRejectsAmbiguousOrMalformedData)
{
	const std::string header = "! RADCODE format=2\n";
	const std::string valid = "Q01 -101 Z01 123 200000 100000 Station\n";
	for (const auto& text :
	     std::vector<std::string>{"",
	                              valid,
	                              "! RADCODE format=1\n" + valid,
	                              header,
	                              header + header + valid,
	                              header + valid + valid,
	                              header + "! ALIAS -1 Q01\n" + valid,
	                              header + valid + "Q02 -101 Z02 42 200000 100000 Duplicate JPL\n",
	                              header + valid + "Q02 -102 Z01 42 200000 100000 Duplicate MPC\n",
	                              header + valid + "Z01 -102 Z02 42 200000 100000 Canonical collision\n",
	                              header + valid + "Q02 -102 Q01 42 200000 100000 MPC collision\n",
	                              header + "Q1 -101 Z01 123 200000 100000 Invalid key\n",
	                              header + "Q01 101 Z01 123 200000 100000 Invalid JPL\n",
	                              header + "Q01 -0101 Z01 123 200000 100000 Invalid JPL\n",
	                              header + "Q01 -0 Z01 123 200000 100000 Invalid JPL\n",
	                              header + "Q01 -101 Z1 123 200000 100000 Invalid MPC\n",
	                              header + "Q01 -101 Z01 360 200000 100000 Longitude\n",
	                              header + "Q01 -101 Z01 -1 200000 100000 Longitude\n",
	                              header + "Q01 -101 Z01 123 -200000 100000 RhoXY\n",
	                              header + "Q01 -101 Z01 123 0 0 Origin\n",
	                              header + "Q01 -101 Z01 nan 200000 100000 Nonfinite\n",
	                              header + "Q01 -101 Z01 123 inf 100000 Nonfinite\n",
	                              header + "Q01 -101 Z01 123 200000 1e999 Overflow\n",
	                              header + "Q01 -101 Z01 123 200000 100000suffix\n",
	                              header + "Q01 -101 Z01 123 200000 100000\n",
	                              header + "Q01 123 200000 100000 Legacy row\n"}) {
		SCOPED_TRACE(text);
		std::istringstream input(text);
		EXPECT_THROW((void)obs::readRadarCatalog(input), std::runtime_error);
	}
}

TEST_F(ExporterRwoTest, RadarCoordinatesRetainJplPrecisionAndLegacyStations)
{
	const auto& arecibo = obs::radarObservatory("-1");
	const auto& goldstone = obs::radarObservatory("-14");
	EXPECT_DOUBLE_EQ(arecibo.longitude, 293.2473068);
	EXPECT_DOUBLE_EQ(arecibo.cosPhi, 404849.131831505 * (1e-10 / Constants::EARTH_RADIUS_IN_AU));
	EXPECT_DOUBLE_EQ(arecibo.sinPhi, 133339.192641330 * (1e-10 / Constants::EARTH_RADIUS_IN_AU));
	EXPECT_DOUBLE_EQ(goldstone.longitude, 243.1104618);
	EXPECT_DOUBLE_EQ(goldstone.cosPhi, 347865.711081849 * (1e-10 / Constants::EARTH_RADIUS_IN_AU));
	EXPECT_DOUBLE_EQ(goldstone.sinPhi, 245795.765861793 * (1e-10 / Constants::EARTH_RADIUS_IN_AU));
	for (const auto* code : {"257", "274", "275", "276"})
		EXPECT_TRUE(obs::radarObservatoryCatalog().contains(code));
	EXPECT_THROW((void)obs::radarObservatory("-"), std::invalid_argument);
}

TEST_F(ExporterRwoTest, RadarLoadingAndPreparationUseRadarNamespaceForBothStations)
{
	struct RestoreEphemeris {
		IEphemerisSource* previous{getDefaultEphemerisSource()};

		~RestoreEphemeris() { setDefaultEphemerisSource(previous); }
	} restore;

	DE440Source ephemeris;
	ASSERT_TRUE(ephemeris.isAvailable());
	setDefaultEphemerisSource(&ephemeris);
	const auto path = root / "stations.xml";
	{
		std::ofstream file(path);
		file << "<ades version=\"2022\">";
		for (int day : {1, 2, 3})
			file << "<optical><permID>99942</permID><stn>500</stn><mode>CCD</mode>"
			     << "<obsTime>2020-01-0" << day << "T00:00:00Z</obsTime><ra>30</ra><dec>10</dec></optical>";
		for (const auto* receiver : {"259", "273", "316", "317"})
			file << "<radar><permID>99942</permID><trx>263</trx><rcv>" << receiver << "</rcv>"
			     << "<obsTime>2020-01-04T00:00:00Z</obsTime><com>1</com><frq>8560</frq>"
			     << "<delay>100</delay><rmsDelay>1</rmsDelay></radar>";
		file << "</ades>";
	}
	const auto loaded = loadObservations(path);
	ASSERT_EQ(loaded.all.size(), 7u);
	ASSERT_EQ(loaded.iod.size(), 3u);
	std::vector<ObsNode> radar;
	for (const auto& record : loaded.all) {
		if (!record.radar)
			continue;
		EXPECT_EQ(record.radar->transmitter, "035");
		EXPECT_EQ(record.obsName, obs::radarObservatory(record.obsCode).name);
		EXPECT_DOUBLE_EQ(record.cosPhi, obs::radarObservatory(record.obsCode).cosPhi);
		radar.push_back(record);
	}
	ASSERT_EQ(radar.size(), 4u);
	const auto prepared = fit::prepareFitObservations(radar, OrbitOptimizerConfig{});
	ASSERT_EQ(prepared.size(), radar.size());
	for (std::size_t i = 0; i < radar.size(); ++i) {
		EXPECT_EQ(prepared[i].observer.obsCode, radar[i].obsCode);
		EXPECT_EQ(prepared[i].transmitter.obsCode, "035");
		EXPECT_DOUBLE_EQ(prepared[i].observer.cosPhi, radar[i].cosPhi);
		EXPECT_DOUBLE_EQ(prepared[i].transmitter.cosPhi, obs::radarObservatory("-35").cosPhi);
		EXPECT_GT(prepared[i].observerOffsetEclAu.magnitude(), 0.0);
	}
}

TEST_F(ExporterRwoTest, PhotometryPreservesAdesBandSigmaAndFullPrecisionH)
{
	for (const auto mag : {21.614, 19.1}) {
		auto row = addObservation("T12", 0.2);
		const bool ades = observations.size() == 1;
		auto& obs = observations.back();
		obs.MAG = row.mag = mag;
		obs.magBand = row.magBand = ades ? 'G' : 'V';
		obs.magString = ades ? "21.614" : "19.1";
		obs.rmsMag = ades ? 0.087 : 0.0;
		row.magSigma = photometry::resolveSigma(obs.magString, obs.rmsMag);
		row.helioDistAu = row.geoDistAu = 1.0;
		fit.residuals.add(row);
	}
	const auto result = save();
	ASSERT_EQ(result.rows.size(), 2u);
	EXPECT_EQ(result.rows[0].magnitude, "21.61G");
	EXPECT_EQ(result.rows[1].magnitude, "19.1 V");
	EXPECT_DOUBLE_EQ(result.rows[0].magSigma, 0.09);
	EXPECT_DOUBLE_EQ(result.rows[1].magSigma, 0.7);
	const auto solution = photometry::fit(fit.residuals.data());
	ASSERT_TRUE(solution.absoluteMagnitude);
	for (std::size_t i = 0; i < result.rows.size(); ++i) {
		const auto expected = photometry::evaluate(fit.residuals.data()[i], solution);
		ASSERT_TRUE(result.rows[i].magResidual);
		ASSERT_TRUE(expected.residual);
		EXPECT_NEAR(*result.rows[i].magResidual, *expected.residual, 0.0051);
	}

	ExportContext context;
	context.objectName = "2020VV";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	const auto report = Exporter::save(context, {ExportProduct::OrbitElements});
	std::ifstream oef(report.directory / "orbit_equinoctial_oef.txt");
	bool found = false;
	for (std::string line; std::getline(oef, line);) {
		if (line.starts_with(" MAG")) {
			EXPECT_NEAR(std::stod(line.substr(4)), *solution.absoluteMagnitude, 0.00051);
			found = true;
		}
	}
	EXPECT_TRUE(found);
}

TEST_F(ExporterRwoTest, CometResidualsUseCometLawAndDoNotClassifyUnknownMagnitudes)
{
	for (int i = 0; i < 3; ++i) {
		auto row = addObservation("309", 0.2);
		auto& obs = observations.back();
		row.helioDistAu = i == 1 ? 2.0 : 1.0;
		row.geoDistAu = 1.0;
		row.mag = obs.MAG = i == 2 ? 20.0 : 15.0 + 10.0 * std::log10(row.helioDistAu);
		row.magBand = obs.magBand = i == 2 ? 'V' : 'T';
		obs.magString = i == 1 ? "18.01" : i == 2 ? "20.0" : "15.0";
		row.magSigma = 0.5;
		if (i != 2)
			row.nuclearMagnitude = false;
		fit.residuals.add(row);
	}
	ExportContext context;
	context.objectName = "2020VV";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	context.cometPhotometry = true;
	const auto report = Exporter::save(context, {ExportProduct::Residuals});
	const auto result = readRwo(report.directory / "weighted_residuals_rwo.txt");
	ASSERT_EQ(result.rows.size(), 3u);
	ASSERT_TRUE(result.rows[0].magResidual);
	ASSERT_TRUE(result.rows[1].magResidual);
	EXPECT_NEAR(*result.rows[0].magResidual, 0.0, 0.0051);
	EXPECT_NEAR(*result.rows[1].magResidual, 0.0, 0.0051);
	EXPECT_FALSE(result.rows[2].magResidual);
}

TEST_F(ExporterRwoTest, AstrometricRejectionDoesNotDiscardPhotometryOrItsLeaveOneOutCheck)
{
	for (std::size_t i = 0; i < 6; ++i) {
		auto row = addObservation("309", 0.2);
		row.helioDistAu = row.geoDistAu = 1.0;
		row.magBand = observations.back().magBand = 'V';
		row.magSigma = 1.0;
		row.mag = observations.back().MAG = i < 5 ? 20.0 : 23.0;
		fit.allResiduals.add(row);
		if (i < 5) {
			fit.residuals.add(row);
			fit.keptObsIndices.push_back(i);
		}
	}
	const auto result = save();
	ASSERT_EQ(result.rows.size(), 6u);
	ASSERT_TRUE(result.rows[0].magResidual);
	ASSERT_TRUE(result.rows[5].magResidual);
	EXPECT_DOUBLE_EQ(*result.rows[0].magResidual, -0.5);
	EXPECT_DOUBLE_EQ(*result.rows[5].magResidual, 2.5);
	EXPECT_EQ(result.rows[0].astSelected, 1);
	EXPECT_EQ(result.rows[0].magSelected, 1);
	EXPECT_EQ(result.rows[5].astSelected, 0);
	EXPECT_EQ(result.rows[5].magSelected, 0);
}

TEST_F(ExporterRwoTest, CometPhotometryIncludesAstrometricallyRejectedMeasurements)
{
	for (std::size_t i = 0; i < 2; ++i) {
		auto row = addObservation("309", 0.2);
		row.helioDistAu = row.geoDistAu = 1.0;
		row.magBand = observations.back().magBand = 'T';
		row.nuclearMagnitude = false;
		row.magSigma = 1.0;
		row.mag = observations.back().MAG = i == 0 ? 15.0 : 17.0;
		observations.back().magString = i == 0 ? "15.0" : "17.0";
		fit.allResiduals.add(row);
		if (i == 0) {
			fit.residuals.add(row);
			fit.keptObsIndices.push_back(i);
		}
	}
	ExportContext context;
	context.objectName = "2020VV";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	context.cometPhotometry = true;
	const auto report = Exporter::save(context, {ExportProduct::Residuals});
	const auto result = readRwo(report.directory / "weighted_residuals_rwo.txt");
	ASSERT_EQ(result.rows.size(), 2u);
	ASSERT_TRUE(result.rows[0].magResidual);
	ASSERT_TRUE(result.rows[1].magResidual);
	EXPECT_NEAR(*result.rows[0].magResidual, -1.0, 0.0051);
	EXPECT_NEAR(*result.rows[1].magResidual, 1.0, 0.0051);
	EXPECT_EQ(result.rows[1].astSelected, 0);
	EXPECT_EQ(result.rows[1].magSelected, 1);
}

TEST_F(ExporterRwoTest, PrefersFinalSelectedMetricOverFullSetDiagnosticPass)
{
	auto selected = addObservation("309", 0.02);
	auto rejected = addObservation("G37", 0.08);
	auto diagnosticSelected = selected;
	diagnosticSelected.astrometricWeight->sigmaRaArcsec = 0.20;
	diagnosticSelected.astrometricWeight->sigmaDecArcsec = 0.20;
	diagnosticSelected.astrometricWeight->chiSquared = 100.0;
	rejected.astrometricWeight->chiSquared = 10000.0;
	fit.residuals.add(selected);
	fit.allResiduals.add(diagnosticSelected);
	fit.allResiduals.add(rejected);
	fit.keptObsIndices = {0};
	const auto result = save();
	ASSERT_EQ(result.rows.size(), 2u);
	EXPECT_DOUBLE_EQ(result.rows[0].sigmaRa, 0.02);
	EXPECT_DOUBLE_EQ(result.rows[1].sigmaRa, 0.08);
	EXPECT_NEAR(result.rms, std::sqrt(5.0 / 2.0), 5e-6);
}

TEST_F(ExporterRwoTest, DoesNotReapplyObserverSigmaOrReweightStoredMetric)
{
	auto row = addObservation("G37", 0.08);
	row.astrometricWeight = ResidualWeight{0.16, 0.24, 0.5, true, false, true};
	fit.residuals.add(row);
	const auto result = save();
	ASSERT_EQ(result.rows.size(), 1u);
	EXPECT_DOUBLE_EQ(result.rows[0].sigmaRa, 0.16);
	EXPECT_DOUBLE_EQ(result.rows[0].sigmaDec, 0.24);
	EXPECT_EQ(result.rows[0].forceRa, 'T');
	EXPECT_EQ(result.rows[0].forceDec, 'F');
	EXPECT_NEAR(result.rows[0].chi, std::sqrt(0.5), 0.005);
	EXPECT_DOUBLE_EQ(result.rms, 0.5);
}

TEST_F(ExporterRwoTest, ProjectsOnlyPrintedRaBiasToGreatCircle)
{
	fit.residuals.add(addObservation("309", 0.02));
	auto& observation = observations[0];
	observation.biasRaArcsec = 2.0;
	observation.biasDecArcsec = 0.3;
	observation.DEC = Constants::PI / 3.0 - 0.3 * Constants::ARCSEC_RAD;
	const auto result = save();
	ASSERT_EQ(result.rows.size(), 1u);
	EXPECT_DOUBLE_EQ(result.rows[0].biasRa, 1.0);
	EXPECT_DOUBLE_EQ(result.rows[0].biasDec, 0.3);
	EXPECT_DOUBLE_EQ(observation.biasRaArcsec, 2.0);
}

TEST_F(ExporterRwoTest, PostFitChiDoesNotReplaceMeasurementRmsOrSigmas)
{
	auto row = addObservation("309", 0.02);
	row.astrometricWeight->postFitChiSquared = 16.0;
	fit.residuals.add(row);
	const auto result = save();
	ASSERT_EQ(result.rows.size(), 1u);
	EXPECT_DOUBLE_EQ(result.rows[0].chi, 4.0);
	EXPECT_DOUBLE_EQ(result.rows[0].sigmaRa, 0.02);
	EXPECT_DOUBLE_EQ(result.rows[0].sigmaDec, 0.02);
	EXPECT_NEAR(result.rms, std::sqrt(5.0 / 2.0), 5e-6);
}

TEST_F(ExporterRwoTest, RefusesToInventMissingFitWeights)
{
	auto row = addObservation("309", 0.02);
	row.astrometricWeight.reset();
	fit.residuals.add(row);
	EXPECT_THROW((void)save(), std::runtime_error);
}

TEST_F(ExporterRwoTest, PreservesSmallSigmasAndOriginalAdesMetadata)
{
	auto residual = addObservation("F51", 0.0000123456789);
	fit.residuals.add(residual);
	auto& observation = observations.front();
	observation.astCat = "Gaia3E-original";
	observation.adesProgramCode = "abc";
	observation.notes = "quoted \"note\"";
	observation.hasRmsTime = true;
	observation.rmsTimeSec = 0.123456789;
	observation.rmsCorr = -0.1234567;
	observation.nuclearMagnitude = true;
	ExportContext context;
	context.objectName = "2020VV";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	const auto exported = Exporter::save(context, {ExportProduct::Residuals});
	std::ifstream input(exported.directory / "weighted_residuals_rwo.txt");
	const std::string text(std::istreambuf_iterator<char>(input), {});
	obs::RwoReader reader;
	std::string error;
	const auto restored = reader.read(text, error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(restored.size(), 1u);
	EXPECT_EQ(restored.front().astCat, observation.astCat);
	EXPECT_EQ(restored.front().notes, observation.notes);
	EXPECT_DOUBLE_EQ(restored.front().raRad, observation.RA);
	EXPECT_DOUBLE_EQ(restored.front().rmsRaArcsec, observation.rmsRaArcsec);
	EXPECT_DOUBLE_EQ(restored.front().rmsCorr, observation.rmsCorr);
	EXPECT_DOUBLE_EQ(restored.front().rmsTimeSec, observation.rmsTimeSec);
	ASSERT_TRUE(restored.front().nuclearMagnitude.has_value());
	EXPECT_TRUE(*restored.front().nuclearMagnitude);
	std::istringstream lines(text);
	std::string native;
	for (std::string line; std::getline(lines, line);)
		if (!line.starts_with("! OC_INPUT"))
			native += line + '\n';
	const auto standard = reader.read(native, error);
	ASSERT_TRUE(error.empty()) << error;
	EXPECT_GT(standard.front().rmsRaArcsec, 0);
	EXPECT_NEAR(standard.front().rmsRaArcsec, observation.rmsRaArcsec, 5e-8);
	auto corrupted = text;
	const auto version = corrupted.find("\"v\":1");
	ASSERT_NE(version, std::string::npos);
	corrupted[version + 4] = '2';
	EXPECT_TRUE(reader.read(corrupted, error).empty());
}

TEST_F(ExporterRwoTest, RejectsInvalidSnapshotRatherThanFallingBackToRules)
{
	auto row = addObservation("309", 0.02);
	row.astrometricWeight->sigmaRaArcsec = 0.0;
	fit.residuals.add(row);
	EXPECT_THROW((void)save(), std::runtime_error);
}

TEST_F(ExporterRwoTest, HistoricalEpochPreservesPhysicalInstantAndSubmittedScale)
{
	getExecutableDir() = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	for (const auto& [epoch, scale] :
	     {std::pair{2426323.692360, time::Scale::UT1}, std::pair{2438942.75, time::Scale::UTC}}) {
		observations.clear();
		fit.residuals = {};
		auto residual = addObservation("690", 1.0);
		auto& observation = observations.front();
		observation.JD = epoch;
		observation.timeScale = scale;
		residual.jd = observation.epochJdTdb();
		fit.residuals.add(residual);
		ExportContext context;
		context.objectName = "1862";
		context.fitResult = &fit;
		context.observations = observations;
		context.rootDirectory = root;
		const auto exported = Exporter::save(context, {ExportProduct::Residuals});
		std::ifstream input(exported.directory / "weighted_residuals_rwo.txt");
		const std::string text(std::istreambuf_iterator<char>(input), {});
		std::string error;
		const auto restored = obs::RwoReader{}.read(text, error);
		ASSERT_TRUE(error.empty()) << error;
		ASSERT_EQ(restored.size(), 1u);
		EXPECT_DOUBLE_EQ(restored[0].jd, epoch);
		EXPECT_EQ(restored[0].timeScale, scale);
		std::istringstream lines(text);
		std::string native;
		for (std::string line; std::getline(lines, line);)
			if (!line.starts_with("! OC_INPUT"))
				native += line + '\n';
		const auto standard = obs::RwoReader{}.read(native, error);
		ASSERT_TRUE(error.empty()) << error;
		ASSERT_EQ(standard.size(), 1u);
		EXPECT_EQ(standard[0].timeScale, time::Scale::UT1);
		const double expectedTt = time::toMjdTt(time::mjdFromJd(epoch), scale);
		const double standardTt = time::toMjdTt(time::mjdFromJd(standard[0].jd), standard[0].timeScale);
		EXPECT_NEAR((standardTt - expectedTt) * 86400., 0., .4321);
	}
}

TEST_F(ExporterRwoTest, EmptyEffectiveSelectionIsNotWrittenAsAllIncluded)
{
	fit.allResiduals.add(addObservation("F51", .2));
	const auto result = save();
	ASSERT_EQ(result.rows.size(), 1u);
	EXPECT_EQ(result.rows.front().astSelected, 0);
}

TEST_F(ExporterRwoTest, LongObjectNamesAndTinyPhotometricErrorsDoNotPreventSaving)
{
	auto row = addObservation("F51", .2);
	row.mag = observations.front().MAG = 18.0;
	observations.front().magString = "18.000";
	observations.front().rmsMag = row.magSigma = .0001;
	row.helioDistAu = row.geoDistAu = 1.0;
	fit.residuals.add(row);
	ExportContext context;
	context.objectName = "2I/Borisov";
	context.fitResult = &fit;
	context.observations = observations;
	context.rootDirectory = root;
	const auto exported = Exporter::save(context, {ExportProduct::Residuals});
	std::ifstream input(exported.directory / "weighted_residuals_rwo.txt");
	const std::string text(std::istreambuf_iterator<char>(input), {});
	obs::RwoReader reader;
	std::string error;
	const auto restored = reader.read(text, error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(restored.size(), 1u);
	EXPECT_EQ(restored.front().objectName, context.objectName);
	EXPECT_DOUBLE_EQ(restored.front().rmsMag, .0001);
}

TEST(ResidualWeightReport, PreservesCorrelatedTimingMetric)
{
	obs::PreparedObservation prepared;
	prepared.sigmaRa = 0.4 * Constants::ARCSEC_RAD;
	prepared.sigmaDec = 0.6 * Constants::ARCSEC_RAD;
	const double l00 = 0.5;
	const double l10 = 0.18 / l00;
	const double l11 = std::sqrt(0.52 - l10 * l10);
	const double arc = Constants::ARCSEC_RAD;
	const obs::AstrometricResidualTransform transform{1.0 / l00 / arc, 0.0, -l10 / (l00 * l11) / arc, 1.0 / l11 / arc,
	                                                  true};
	const auto result = obs::residualWeightForReport(prepared, 0.3 * arc, -0.2 * arc, true, transform);
	EXPECT_NEAR(result.sigmaRaArcsec, 0.5, 1e-14);
	EXPECT_NEAR(result.sigmaDecArcsec, std::sqrt(0.52), 1e-14);
	const double expected =
	    (0.52 * 0.3 * 0.3 + 2.0 * 0.18 * 0.3 * 0.2 + 0.25 * 0.2 * 0.2) / (0.25 * 0.52 - 0.18 * 0.18);
	EXPECT_NEAR(result.chiSquared, expected, 1e-14);
}

TEST(ResidualWeightReport, UnweightedModeIgnoresStationAndObserverSigmas)
{
	obs::PreparedObservation prepared;
	prepared.sigmaRa = prepared.sigmaDec = 0.02 * Constants::ARCSEC_RAD;
	prepared.forceWeightRa = prepared.forceWeightDec = true;
	const auto result =
	    obs::residualWeightForReport(prepared, 3.0 * Constants::ARCSEC_RAD, 4.0 * Constants::ARCSEC_RAD, false);
	EXPECT_DOUBLE_EQ(result.sigmaRaArcsec, 1.0);
	EXPECT_DOUBLE_EQ(result.sigmaDecArcsec, 1.0);
	EXPECT_NEAR(result.chiSquared, 25.0, 1e-14);
	EXPECT_FALSE(result.weighted);
	EXPECT_FALSE(result.forceWeightRa);
	EXPECT_FALSE(result.forceWeightDec);
}
