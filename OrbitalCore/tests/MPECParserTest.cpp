#include <gtest/gtest.h>
#include "MPECParser.h"
#include "obs/AstrometricSigma.h"
#include "obs/Mpc80ColumnReader.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Types.h"

#include <filesystem>

static const std::string VALID_LINE = "CK19K080 C2019 05 30.99860 19 35 09.15 +64 38 33.2 19.4 L L06";

static const std::string VALID_LINE_NEGATIVE_DEC = "CK19K080 C2019 06 01.89280 19 30 55.64 -65 29 49.6 19.1 V 033";

TEST(MPECParser, ParseSingleValidLine)
{
	std::vector<std::string> lines = {VALID_LINE};
	MPECParser parser(lines);
	auto result = parser.parse();

	ASSERT_EQ(parser.error(), MPECParser::Error::NoError);
	ASSERT_EQ(result.size(), 1u);

	const auto& data = result[0];
	EXPECT_EQ(data.objectName, "CK19K080");
	EXPECT_EQ(data.code, "L06");
	EXPECT_EQ(data.ra, "19 35 09.15");
	EXPECT_EQ(data.dec, "+64 38 33.2");
}

TEST(MPECParser, ParseNegativeDec)
{
	std::vector<std::string> lines = {VALID_LINE_NEGATIVE_DEC};
	MPECParser parser(lines);
	auto result = parser.parse();

	ASSERT_EQ(parser.error(), MPECParser::Error::NoError);
	ASSERT_EQ(result.size(), 1u);
	EXPECT_EQ(result[0].dec, "-65 29 49.6");
}

TEST(MPECParser, ParseMultipleLines)
{
	std::vector<std::string> lines = {VALID_LINE, VALID_LINE_NEGATIVE_DEC};
	MPECParser parser(lines);
	auto result = parser.parse();

	ASSERT_EQ(parser.error(), MPECParser::Error::NoError);
	EXPECT_EQ(result.size(), 2u);
}

TEST(MPECParser, JDConversion)
{
	std::vector<std::string> lines = {VALID_LINE};
	MPECParser parser(lines);
	auto result = parser.parse();

	ASSERT_EQ(result.size(), 1u);
	EXPECT_NEAR(result[0].JD, Convertors::convertGregorianToJD(2019, 5, 30.99860), 1e-6);
	EXPECT_TRUE(result[0].hasSplitUtc);
	EXPECT_EQ(result[0].utcMjdDay, 58633);
	EXPECT_DOUBLE_EQ(result[0].utcSecondsOfDay, 0.99860 * 86400.0);
}

TEST(MPECParser, RetainsMpcFractionalDayWithoutFullJdPrecisionLoss)
{
	const std::string line = "     K20V00V KC2020 11 15.25270 01 36 13.03 +13 45 15.5          22.5 Ro~4Mfd291";
	MPECParser parser({line});
	const auto result = parser.parse();

	ASSERT_EQ(parser.error(), MPECParser::Error::NoError);
	ASSERT_EQ(result.size(), 1u);
	EXPECT_EQ(result[0].utcMjdDay, 59168);
	EXPECT_DOUBLE_EQ(result[0].utcSecondsOfDay, 0.25270 * 86400.0);

	const double jdTdt =
	    2400000.5 + (static_cast<double>(result[0].utcMjdDay) + (result[0].utcSecondsOfDay + 69.184) / 86400.0);
	EXPECT_DOUBLE_EQ(jdTdt, 2459168.7535007410);
}

TEST(MPECParser, RetainsPackedCometDesignation)
{
	const std::string line = "    CK20F030* C2020 03 27.71352 07 00 39.56 -48 50 52.1          17   RLEG005C51";
	MPECParser parser({line});
	const auto result = parser.parse();
	ASSERT_EQ(parser.error(), MPECParser::Error::NoError);
	ASSERT_EQ(result.size(), 1u);
	EXPECT_EQ(result[0].classificationDesignation, "    CK20F030");
}

TEST(MPECParser, Magnitude)
{
	std::vector<std::string> lines = {VALID_LINE};
	MPECParser parser(lines);
	auto result = parser.parse();

	ASSERT_EQ(result.size(), 1u);
	EXPECT_EQ(result[0].magnitude, "19.4");
}

TEST(MPECParser, EmptyInput)
{
	std::vector<std::string> lines = {};
	MPECParser parser(lines);
	auto result = parser.parse();

	EXPECT_EQ(parser.error(), MPECParser::Error::InvalidData);
	EXPECT_TRUE(result.empty());
}

TEST(MPECParser, SkipsMalformedLinesAndKeepsValidObservations)
{
	std::vector<std::string> lines = {"CK19K080 C2019 05", VALID_LINE,
	                                  "CK19K080 C2019 05 30.99860 XX 35 09.15 +64 38 33.2 19.4 L L06"};
	MPECParser parser(lines);
	auto result = parser.parse();

	EXPECT_EQ(parser.error(), MPECParser::Error::NoError);
	ASSERT_EQ(result.size(), 1u);
	EXPECT_EQ(result[0].objectName, "CK19K080");
}

TEST(MPECParser, ErrorStrings)
{
	std::vector<std::string> lines = {};
	MPECParser p(lines);
	p.parse();
	EXPECT_EQ(p.errorString(), "Invalid data");
}

namespace {
	std::string buildMPC80(const std::string& number5, const std::string& provDes7, char discAsterisk,
	                       const std::string& notes2, const std::string& date17, const std::string& ra12,
	                       const std::string& dec12, const std::string& magBand6, char catalog,
	                       const std::string& obsCode3)
	{
		auto pad = [](const std::string& s, size_t n) {
			if (s.size() >= n)
				return s.substr(0, n);
			return s + std::string(n - s.size(), ' ');
		};
		std::string line;
		line += pad(number5, 5);
		line += pad(provDes7, 7);
		line += discAsterisk;
		line += pad(notes2, 2);
		line += pad(date17, 17);
		line += pad(ra12, 12);
		line += pad(dec12, 12);
		line += std::string(9, ' ');
		line += pad(magBand6, 6);
		line += catalog;
		line += std::string(5, ' ');
		line += pad(obsCode3, 3);
		return line;
	}
}  // namespace

TEST(MPECParser, MPC80Col_GaiaDR2)
{
	const std::string line = buildMPC80("     ", "K19Q080", ' ', " C", "2019 09 08.48596", "14 51 29.77", "-24 01 47.3",
	                                    " 20.6G", 'V', "I52");
	ASSERT_EQ(line.size(), 80u);

	MPECParser p({line});
	auto res = p.parse();
	ASSERT_EQ(p.error(), MPECParser::Error::NoError);
	ASSERT_EQ(res.size(), 1u);
	EXPECT_EQ(res[0].objectName, "K19Q080");
	EXPECT_EQ(res[0].ra, "14 51 29.77");
	EXPECT_EQ(res[0].dec, "-24 01 47.3");
	EXPECT_EQ(res[0].filter, "V");
	EXPECT_EQ(res[0].magBand, "G");
	EXPECT_EQ(res[0].catalog, "V");
	EXPECT_EQ(res[0].code, "I52");
}

TEST(MPECParser, MPC80Col_GaiaDR3_Positive_Dec)
{
	const std::string line = buildMPC80("     ", "K19Q080", ' ', " C", "2019 10 15.12345", "15 30 10.00", "+05 10 20.5",
	                                    " 20.0R", 'W', "F51");
	ASSERT_EQ(line.size(), 80u);
	MPECParser p({line});
	auto res = p.parse();
	ASSERT_EQ(p.error(), MPECParser::Error::NoError);
	EXPECT_EQ(res[0].filter, "W");
	EXPECT_EQ(res[0].magBand, "R");
	EXPECT_EQ(res[0].catalog, "W");
	EXPECT_EQ(res[0].code, "F51");
	EXPECT_EQ(res[0].dec, "+05 10 20.5");
}

TEST(MPECParser, MPC80Col_BlankCatalog)
{
	const std::string line = buildMPC80("00001", "       ", ' ', " C", "2000 01 01.00000", "12 00 00.00", "+00 00 00.0",
	                                    " 17.5V", ' ', "675");
	ASSERT_EQ(line.size(), 80u);
	MPECParser p({line});
	auto res = p.parse();
	ASSERT_EQ(p.error(), MPECParser::Error::NoError);
	EXPECT_EQ(res[0].filter, "");
	EXPECT_EQ(res[0].magBand, "V");
	EXPECT_EQ(res[0].catalog, "");
	EXPECT_EQ(res[0].code, "675");
	EXPECT_EQ(res[0].objectName, "00001");
}

TEST(MPECParser, MPC80Col_DateConversion)
{
	const std::string line = buildMPC80("     ", "K19Q080", ' ', " C", "2019 09 08.48596", "14 51 29.77", "-24 01 47.3",
	                                    " 20.6G", 'V', "I52");
	MPECParser p({line});
	auto res = p.parse();
	ASSERT_EQ(res.size(), 1u);
	EXPECT_NEAR(res[0].JD, Convertors::convertGregorianToJD(2019, 9, 8.48596), 1e-6);
}

TEST(MPECParser, DispatchesByLineLength)
{
	const std::string shortLine = "CK19K080 C2019 05 30.99860 19 35 09.15 +64 38 33.2 19.4 L L06";
	const std::string longLine = buildMPC80("     ", "K19Q080", ' ', " C", "2019 09 08.48596", "14 51 29.77",
	                                        "-24 01 47.3", " 20.6G", 'V', "I52");

	MPECParser p({shortLine, longLine});
	auto res = p.parse();
	ASSERT_EQ(p.error(), MPECParser::Error::NoError);
	ASSERT_EQ(res.size(), 2u);
	EXPECT_EQ(res[0].filter, "L");
	EXPECT_EQ(res[0].magBand, "L");
	EXPECT_EQ(res[0].catalog, "");
	EXPECT_EQ(res[0].code, "L06");
	EXPECT_EQ(res[1].filter, "V");
	EXPECT_EQ(res[1].magBand, "G");
	EXPECT_EQ(res[1].catalog, "V");
	EXPECT_EQ(res[1].code, "I52");
}

TEST(MpcObservationReader, SignedDeclinationIncludingNegativeZero)
{
	struct Case {
		const char* text;
		double expected;
	};

	const Case cases[] = {
	    {"-00 06 49.7", -Convertors::decToRadians(0, 6, 49.7)},
	    {"+00 06 49.7", Convertors::decToRadians(0, 6, 49.7)},
	    {"-15 01 33.6", -Convertors::decToRadians(15, 1, 33.6)},
	    {"+64 38 33.2", Convertors::decToRadians(64, 38, 33.2)},
	};
	for (const auto& test : cases) {
		const auto line = buildMPC80("     ", "K19Q080", ' ', " C", "2019 09 08.48596", "14 51 29.77", test.text,
		                             " 20.6G", 'V', "I52");
		std::string error;
		const auto observations = orbital_core::obs::Mpc80ColumnReader{}.read(line, error);
		ASSERT_EQ(observations.size(), 1u) << test.text << ": " << error;
		EXPECT_NEAR(observations[0].decRad, test.expected, 1e-12) << test.text;
	}
}

TEST(MPECParser, MPC80Col_WithDiscoveryAsteriskAndBatchRef)
{
	const std::string line = "0002IK19Q040* C2019 08 30.04398 08 26 37.82 +33 39 52.7          17.8 Rq~01XoL51";
	ASSERT_EQ(line.size(), 80u);

	MPECParser p({line});
	auto res = p.parse();
	ASSERT_EQ(p.error(), MPECParser::Error::NoError);
	ASSERT_EQ(res.size(), 1u);
	EXPECT_EQ(res[0].objectName, "0002I");
	EXPECT_EQ(res[0].ra, "08 26 37.82");
	EXPECT_EQ(res[0].dec, "+33 39 52.7");
	EXPECT_EQ(res[0].magnitude, "17.8");
	EXPECT_EQ(res[0].filter, "q");
	EXPECT_EQ(res[0].code, "L51");
	EXPECT_NEAR(res[0].JD, Convertors::convertGregorianToJD(2019, 8, 30.04398), 1e-6);
}

TEST(MpcObservationReader, UserVd3RecordsPreserveNote1TechniqueAndCatalog)
{
	const std::string first = "     K23V03D  C2023 11 04.51081003 51 52.068-00 11 31.16         21.24GV~8sxjF52";
	const std::string second = "     K23V03D KC2023 11 07.50185404 31 48.144+01 03 41.26         18.68GV~7RWOI52";
	ASSERT_EQ(first.size(), 80u);
	ASSERT_EQ(second.size(), 80u);
	std::string error;
	const auto observations = orbital_core::obs::Mpc80ColumnReader{}.read(first + '\n' + second, error);
	ASSERT_EQ(observations.size(), 2u) << error;
	EXPECT_EQ(observations[0].programCode, ' ');
	EXPECT_EQ(observations[1].programCode, 'K');
	EXPECT_EQ(observations[0].obsCode, "F52");
	EXPECT_EQ(observations[1].obsCode, "I52");
	EXPECT_LT(observations[0].decRad, 0.0);
	for (const auto& observation : observations) {
		EXPECT_EQ(observation.catalogCode, 'V');
		EXPECT_EQ(observation.magBand, 'G');
		EXPECT_EQ(orbital_core::obs::mpcNote2(observation.technique), 'C');
	}
}

TEST(MpcObservationReader, LegacyTokenizedRecordRetainsBothNotes)
{
	std::string error;
	const auto observations = orbital_core::obs::Mpc80ColumnReader{}.read(
	    "K23V03D KC2023 11 07.501854 04 31 48.144 +01 03 41.26 18.68 G I52", error);
	ASSERT_EQ(observations.size(), 1u) << error;
	EXPECT_EQ(observations[0].programCode, 'K');
	EXPECT_EQ(orbital_core::obs::mpcNote2(observations[0].technique), 'C');
}

namespace {
	const std::string satelliteParent =
	    "00001         S2022 01 26.26522 03 43 52.06 +19 11 10.9                L~5N2wC51";
	const std::string satelliteChild =
	    "00001         s2022 01 26.26522 1 + 3146.5614 + 5504.8439 + 2584.0547   ~5N2wC51";
}  // namespace

TEST(MpcObservationReader, SatellitePairUsesExplicitKilometerUnits)
{
	ASSERT_EQ(satelliteParent.size(), 80u);
	ASSERT_EQ(satelliteChild.size(), 80u);
	std::string error;
	const auto observations =
	    orbital_core::obs::Mpc80ColumnReader{}.read(satelliteParent + '\n' + satelliteChild, error);
	ASSERT_EQ(observations.size(), 1u) << error;
	EXPECT_TRUE(observations[0].isSatellite);
	EXPECT_EQ(observations[0].offsetUnits, orbital_core::obs::OffsetUnits::Km);
	EXPECT_DOUBLE_EQ(observations[0].offsetX, 3146.5614);
	EXPECT_DOUBLE_EQ(observations[0].offsetY, 5504.8439);
	EXPECT_DOUBLE_EQ(observations[0].offsetZ, 2584.0547);
	EXPECT_EQ(observations[0].obsCode, "C51");
	EXPECT_EQ(observations[0].catalogCode, 'L');
	EXPECT_TRUE(observations[0].magnitudeText.empty());

	const auto previousDir = getExecutableDir();
	const auto previousSource = getDefaultEphemerisSource();

	struct RestoreEnvironment {
		std::filesystem::path directory;
		IEphemerisSource* source;

		~RestoreEnvironment()
		{
			setDefaultEphemerisSource(source);
			getExecutableDir() = directory;
		}
	} restore{previousDir, previousSource};

	getExecutableDir() = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	DE440Source source;
	ASSERT_TRUE(source.isAvailable());
	setDefaultEphemerisSource(&source);
	const ObsNode node(observations[0], ObservatoryRecord{.code = "C51"});
	EXPECT_TRUE(node.isSatellite);
	EXPECT_EQ(node.note2, 'S');
	const Vec3D expectedOffset =
	    Convertors::ecvToEcl(Vec3D(3146.5614, 5504.8439, 2584.0547) / 149597870.7, Constants::EARTH_INCL);
	auto geocentricData = observations[0];
	geocentricData.isSatellite = false;
	const ObsNode geocenter(geocentricData, ObservatoryRecord{.code = "C51"});
	for (int i = 0; i < 3; ++i) {
		EXPECT_NEAR(node.satelliteOffsetEclAu[i], expectedOffset[i], 1e-16);
		EXPECT_NEAR(node.obsVector[i] - geocenter.obsVector[i], -expectedOffset[i], 1e-14);
	}
	const auto uncertainty =
	    orbital_core::obs::lookupUncertainty(node, node.JD, orbital_core::AstrometricWeightProfile::Vfcc17);
	EXPECT_TRUE(uncertainty.resolved);
	EXPECT_NEAR(uncertainty.sigmaRaRad / Constants::ARCSEC_RAD, 1.0, 1e-14);
	EXPECT_NEAR(uncertainty.sigmaDecRad / Constants::ARCSEC_RAD, 1.0, 1e-14);
	ObsNode genericSatellite = node;
	genericSatellite.obsCode = "999";
	const auto genericUncertainty =
	    orbital_core::obs::lookupUncertainty(genericSatellite, node.JD, orbital_core::AstrometricWeightProfile::Vfcc17);
	EXPECT_NEAR(genericUncertainty.sigmaRaRad / Constants::ARCSEC_RAD, 1.5, 1e-14);
	EXPECT_NEAR(genericUncertainty.sigmaDecRad / Constants::ARCSEC_RAD, 1.5, 1e-14);
}

TEST(MpcObservationReader, SatellitePairUsesExplicitAuUnits)
{
	std::string child = satelliteChild;
	child[32] = '2';
	const std::string payload = "+0.000021 -0.000038 +0.000015";
	child.replace(33, 36, 36, ' ');
	child.replace(33, payload.size(), payload);
	std::string error;
	const auto observations = orbital_core::obs::Mpc80ColumnReader{}.read(satelliteParent + '\n' + child, error);
	ASSERT_EQ(observations.size(), 1u) << error;
	EXPECT_EQ(observations[0].offsetUnits, orbital_core::obs::OffsetUnits::Au);
	EXPECT_DOUBLE_EQ(observations[0].offsetX, 0.000021);
	EXPECT_DOUBLE_EQ(observations[0].offsetY, -0.000038);
}

TEST(MpcObservationReader, RovingPairPreservesGeodeticPositionAndUsesGroundGeometry)
{
	std::string parent = satelliteParent;
	parent[14] = 'V';
	parent.replace(77, 3, "247");
	std::string child = parent;
	child[14] = 'v';
	child.replace(32, 45, 45, ' ');
	const std::string coordinates = " 15.250000 -28.500000 1500";
	child.replace(33, coordinates.size(), coordinates);
	std::string error;
	const auto records = orbital_core::obs::Mpc80ColumnReader{}.read(parent + '\n' + child, error);
	ASSERT_EQ(records.size(), 1u) << error;
	const auto& record = records[0];
	ASSERT_TRUE(record.rovingObserver);
	EXPECT_FALSE(record.isSatellite);
	EXPECT_DOUBLE_EQ(record.rovingObserver->longitudeDeg, 15.25);
	EXPECT_DOUBLE_EQ(record.rovingObserver->latitudeDeg, -28.5);
	EXPECT_DOUBLE_EQ(record.rovingObserver->heightMeters, 1500.0);
	const auto previousDir = getExecutableDir();
	const auto previousSource = getDefaultEphemerisSource();

	struct Restore {
		std::filesystem::path directory;
		IEphemerisSource* source;

		~Restore()
		{
			setDefaultEphemerisSource(source);
			getExecutableDir() = directory;
		}
	} restore{previousDir, previousSource};

	getExecutableDir() = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	DE440Source source;
	setDefaultEphemerisSource(&source);
	const ObsNode node(record, ObservatoryRecord{.code = "247"});
	ASSERT_TRUE(node.rovingObserver);
	EXPECT_FALSE(node.isSatellite);
	EXPECT_EQ(node.note2, 'V');
	const double latitude = -28.5 * Constants::DEG_TO_RAD;
	const double e2 = Constants::EARTH_FLATTENING * (2.0 - Constants::EARTH_FLATTENING);
	const double radius = Constants::EARTH_EQ_RADIUS_KM / std::sqrt(1.0 - e2 * std::pow(std::sin(latitude), 2));
	const double cosPhi = (radius + 1.5) * std::cos(latitude) / Constants::EARTH_PARALLAX_RADIUS_KM;
	const double sinPhi = (radius * (1.0 - e2) + 1.5) * std::sin(latitude) / Constants::EARTH_PARALLAX_RADIUS_KM;
	auto ground = record;
	ground.rovingObserver.reset();
	const ObsNode reference(ground,
	                        ObservatoryRecord{.code = "247", .longitude = 15.25, .cosPhi = cosPhi, .sinPhi = sinPhi});
	for (int i = 0; i < 3; ++i)
		EXPECT_DOUBLE_EQ(node.obsVector[i], reference.obsVector[i]);
}

TEST(MpcObservationReader, UnpairedOrMismatchedSatelliteIsNeverGroundObservation)
{
	std::string wrongEpoch = satelliteParent;
	wrongEpoch.replace(wrongEpoch.find("26.26522"), 8, "26.13406");
	std::string wrongStation = satelliteChild;
	wrongStation.replace(77, 3, "274");
	std::string wrongObject = satelliteChild;
	wrongObject.replace(0, 5, "00002");
	std::string badUnits = satelliteChild;
	badUnits[32] = '9';
	for (const auto& text : {satelliteParent, satelliteChild + '\n' + wrongEpoch, wrongEpoch + '\n' + satelliteChild,
	                         satelliteParent + '\n' + wrongStation, satelliteParent + '\n' + wrongObject,
	                         satelliteParent + '\n' + badUnits}) {
		std::string error;
		EXPECT_TRUE(orbital_core::obs::Mpc80ColumnReader{}.read(text, error).empty());
		EXPECT_FALSE(error.empty());
	}
}

TEST(MpcObservationReader, TechniqueCodesSurviveNeutralRecord)
{
	for (char code : std::string("cAHXx")) {
		std::string line = "     K23V03D  C2023 11 04.51081003 51 52.068-00 11 31.16         21.24GV~8sxjF52";
		line[14] = code;
		std::string error;
		const auto observations = orbital_core::obs::Mpc80ColumnReader{}.read(line, error);
		ASSERT_EQ(observations.size(), 1u) << error;
		EXPECT_EQ(orbital_core::obs::mpcNote2(observations[0].technique), code);
	}
}
