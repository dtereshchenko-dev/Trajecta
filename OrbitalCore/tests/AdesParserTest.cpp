#include <gtest/gtest.h>

#include "AdesParser.h"
#include "obs/Mpc80ColumnReader.h"
#include "obs/AstrometricSigma.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/ObservationSet.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Types.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <unordered_set>

namespace {
	std::string xml(const std::string& extra)
	{
		return "<ades version=\"2022\"><optical><provID>2020 VV</provID>"
		       "<stn>309</stn><mode>CCD</mode><obsTime>2020-11-12T04:38:36.067Z</obsTime>"
		       "<ra>10.123456</ra><dec>8.123456</dec><astCat>Gaia2</astCat>" +
		       extra + "</optical></ades>";
	}
}  // namespace

TEST(AdesParser, PermanentInterstellarIdentitySurvivesXmlAndPsvParsing)
{
	for (const auto& [permanentId, expected] : std::initializer_list<std::pair<std::string, bool>>{{"1I", true},
	                                                                                               {"2I", true},
	                                                                                               {"3I", true},
	                                                                                               {"42I", true},
	                                                                                               {"1P", false},
	                                                                                               {"99942", false},
	                                                                                               {"I", false},
	                                                                                               {"0I", false},
	                                                                                               {"2IA", false},
	                                                                                               {"K19Q04I", false},
	                                                                                               {"", false}}) {
		auto xmlText = xml("<permID>" + permanentId + "</permID><trkSub>arbitrary-tracklet</trkSub>");
		const std::string psv =
		    "# version=2022\n"
		    "permID|provID|trkSub|stn|mode|obsTime|ra|dec\n" +
		    permanentId + "|C/2019 Q4|arbitrary-tracklet|309|CCD|2020-11-12T04:38:36.067Z|10.123456|8.123456\n";
		for (const auto& text : {xmlText, psv}) {
			std::string error;
			const auto rows = orbital_core::ades::AdesReader{}.read(text, error);
			ASSERT_EQ(rows.size(), 1u) << error;
			EXPECT_EQ(rows[0].isInterstellar, expected) << permanentId;
			EXPECT_EQ(rows[0].objectName, "arbitrary-tracklet");
			EXPECT_EQ(rows[0].permanentDesignation, permanentId);
			if (!permanentId.empty())
				EXPECT_EQ(rows[0].classificationDesignation, permanentId);
		}
	}
}

TEST(AdesParser, ClassificationUsesDesignationRatherThanTracklet)
{
	const std::string psv =
	    "provID|trkSub|stn|mode|obsTime|ra|dec\n"
	    "C/2019 Q4|2I|309|CCD|2020-11-12T04:38:36.067Z|10.123456|8.123456\n";
	std::string error;
	const auto rows = orbital_core::ades::AdesReader{}.read(psv, error);
	ASSERT_EQ(rows.size(), 1u) << error;
	EXPECT_EQ(rows[0].objectName, "2I");
	EXPECT_EQ(rows[0].classificationDesignation, "C/2019 Q4");
}

TEST(AdesParser, TrackletNameDoesNotDetermineInterstellarClassification)
{
	for (const std::string name : {"2I", "0002I", "interstellar-candidate"}) {
		std::string error;
		const auto rows = orbital_core::ades::AdesReader{}.read(xml("<trkSub>" + name + "</trkSub>"), error);
		ASSERT_EQ(rows.size(), 1u) << error;
		EXPECT_FALSE(rows[0].isInterstellar);
	}
}

TEST(AdesParser, Wgs84PositionIsTerrestrialInXmlAndPsv)
{
	const std::string psv =
	    "provID|stn|mode|obsTime|ra|dec|sys|ctr|pos1|pos2|pos3\n"
	    "2020 VV|247|CCD|2020-11-12T04:38:36.067Z|10.123456|8.123456|WGS84|399|15.25|-28.5|1570.3\n";
	for (const auto& text :
	     {xml("<sys>WGS84</sys><ctr>399</ctr><pos1>15.25</pos1><pos2>-28.5</pos2><pos3>1570.3</pos3>"), psv}) {
		std::string error;
		const auto records = orbital_core::ades::AdesReader{}.read(text, error);
		ASSERT_EQ(records.size(), 1u) << error;
		ASSERT_TRUE(records[0].rovingObserver);
		EXPECT_FALSE(records[0].isSatellite);
		EXPECT_DOUBLE_EQ(records[0].rovingObserver->longitudeDeg, 15.25);
		EXPECT_DOUBLE_EQ(records[0].rovingObserver->latitudeDeg, -28.5);
		EXPECT_DOUBLE_EQ(records[0].rovingObserver->heightMeters, 1570.3);
	}
}

TEST(AdesParser, InvalidOrUnsupportedObserverFrameIsNotTreatedAsIcrf)
{
	for (const auto& extra : {"<sys>WGS84</sys><pos1>15</pos1><pos2>91</pos2><pos3>0</pos3>",
	                          "<sys>ITRF</sys><pos1>1</pos1><pos2>2</pos2><pos3>3</pos3>",
	                          "<sys>ICRF_KM</sys><ctr>10</ctr><pos1>1</pos1><pos2>2</pos2><pos3>3</pos3>",
	                          "<sys>WGS84</sys><pos1>15</pos1><pos2>20</pos2>"}) {
		std::string error;
		EXPECT_TRUE(orbital_core::ades::AdesReader{}.read(xml(extra), error).empty());
		EXPECT_FALSE(error.empty());
	}
}

TEST(AdesParser, OpticalXmlAndPsvPreserveFractionalUtcSeconds)
{
	const std::string psv =
	    "provID|stn|mode|obsTime|ra|dec\n"
	    "2020 VV|309|CCD|2020-11-12T04:38:36.067Z|10.123456|8.123456\n";
	for (const auto& text : {xml(""), psv}) {
		std::string error;
		const auto records = orbital_core::ades::AdesReader{}.read(text, error);
		ASSERT_TRUE(error.empty()) << error;
		ASSERT_EQ(records.size(), 1u);
		EXPECT_FALSE(records[0].radar);
		EXPECT_TRUE(records[0].hasSplitUtc);
		EXPECT_EQ(records[0].utcMjdDay, 59165);
		EXPECT_DOUBLE_EQ(records[0].utcSecondsOfDay, 16716.067);
		EXPECT_NEAR(records[0].timeAccuracyDays, 1e-3 / 86400.0, 1e-22);
	}
}

TEST(AdesParser, InterstellarFlagSurvivesObservationConstructionAndCopy)
{
	getExecutableDir() = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	DE440Source source;
	ASSERT_TRUE(source.isAvailable());
	setDefaultEphemerisSource(&source);

	struct ResetSource {
		~ResetSource() { setDefaultEphemerisSource(nullptr); }
	} reset;

	std::string error;
	const auto parsed = orbital_core::ades::AdesReader{}.read(xml("<permID>2I</permID><trkSub>0002I</trkSub>"), error);
	ASSERT_EQ(parsed.size(), 1u) << error;
	const ObsNode node(parsed[0], ObservatoryRecord{.code = "309"});
	const auto copy = node;
	EXPECT_TRUE(node.isInterstellar);
	EXPECT_TRUE(copy.isInterstellar);
}

TEST(AdesParser, PreservesOriginalCatalogAndUsesSharedCodesInXmlAndPsv)
{
	for (const auto& [name, code] : std::initializer_list<std::pair<std::string, char>>{{"GaIa3E", 'X'},
	                                                                                    {"ZZCAT", '7'},
	                                                                                    {"APASS", '8'},
	                                                                                    {"UBSC", '5'},
	                                                                                    {"AKARI", '!'},
	                                                                                    {"AG", '@'},
	                                                                                    {"WISE", '#'},
	                                                                                    {"LSST2502", '$'},
	                                                                                    {"IRASPSC", '^'},
	                                                                                    {"SDSS9", '&'},
	                                                                                    {"SDSS16", '*'},
	                                                                                    {"FutureCatalog_DR9", ' '},
	                                                                                    {"UNK", ' '},
	                                                                                    {"", ' '}}) {
		auto xmlText = xml("");
		const auto begin = xmlText.find("Gaia2");
		xmlText.replace(begin, 5, name);
		const auto psv =
		    "# version=2022\nprovID|stn|mode|obsTime|ra|dec|astCat\n"
		    "2020 VV|309|CCD|2020-11-12T04:38:36.067Z|10.123456|8.123456|" +
		    name + "\n";
		for (const auto& text : {xmlText, psv}) {
			std::string error;
			const auto rows = orbital_core::ades::AdesReader{}.read(text, error);
			ASSERT_EQ(rows.size(), 1u) << error;
			EXPECT_EQ(rows[0].catalogCode, code) << name;
			EXPECT_EQ(rows[0].astCat, name);
		}
	}
	auto escaped = xml("");
	escaped.replace(escaped.find("Gaia2"), 5, "Custom&amp;Catalog");
	std::string error;
	const auto rows = orbital_core::ades::AdesReader{}.read(escaped, error);
	ASSERT_EQ(rows.size(), 1u) << error;
	EXPECT_EQ(rows[0].astCat, "Custom&Catalog");
	EXPECT_EQ(rows[0].catalogCode, ' ');
}

TEST(AdesParser, OriginalCatalogSurvivesObservationConstructionAndCopy)
{
	getExecutableDir() = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	DE440Source source;
	ASSERT_TRUE(source.isAvailable());
	setDefaultEphemerisSource(&source);

	struct ResetSource {
		~ResetSource() { setDefaultEphemerisSource(nullptr); }
	} reset;

	for (const std::string name : {"GaIa2", "FutureCatalog_DR9", "UNK", ""}) {
		auto text = xml("");
		text.replace(text.find("Gaia2"), 5, name);
		std::string error;
		auto parsed = orbital_core::ades::AdesReader{}.read(text, error);
		ASSERT_EQ(parsed.size(), 1u) << error;
		const ObsNode node(parsed[0], ObservatoryRecord{.code = "309"});
		parsed.clear();
		auto copy = node;
		auto moved = std::move(copy);
		EXPECT_EQ(node.astCat, name);
		EXPECT_EQ(moved.astCat, name);
		EXPECT_EQ(node.catalogCode, name == "GaIa2" ? 'V' : ' ');
	}
}

TEST(AdesParser, ProgramIsDecodedIndependentlyOfNotesInXmlAndPsv)
{
	orbital_core::ades::AdesReader reader;
	const std::string psv =
	    "# version=2022\n"
	    "provID|stn|mode|obsTime|ra|dec|astCat|notes|prog\n"
	    "2020 VV|309|CCD|2020-11-12T04:38:36.067Z|10.123456|8.123456|Gaia2|Kst|0F\n";
	for (const auto& text : {xml("<prog>0F</prog><notes>Kst</notes>"), psv}) {
		std::string error;
		const auto observations = reader.read(text, error);
		ASSERT_EQ(observations.size(), 1u) << error;
		EXPECT_EQ(observations[0].programCode, '&');
		EXPECT_EQ(observations[0].adesProgramCode, "0F");
		EXPECT_EQ(observations[0].notes, "Kst");
	}
}

TEST(AdesParser, PhotometricSigmaAndPrecisionSurviveXmlAndPsvParsing)
{
	const std::string psv =
	    "# version=2022\n"
	    "provID|stn|mode|obsTime|ra|dec|mag|band|rmsMag|nucMag\n"
	    "2020 VV|309|CCD|2020-11-12T04:38:36.067Z|10.123456|8.123456|21.614|G|0.087|1\n";
	for (const auto& text : {xml("<mag>21.614</mag><band>G</band><rmsMag>0.087</rmsMag><nucMag>1</nucMag>"), psv}) {
		std::string error;
		const auto observations = orbital_core::ades::AdesReader{}.read(text, error);
		ASSERT_EQ(observations.size(), 1u) << error;
		EXPECT_EQ(observations[0].magnitudeText, "21.614");
		EXPECT_EQ(observations[0].magBand, 'G');
		EXPECT_DOUBLE_EQ(observations[0].rmsMag, 0.087);
		ASSERT_TRUE(observations[0].nuclearMagnitude.has_value());
		EXPECT_TRUE(*observations[0].nuclearMagnitude);
		EXPECT_FALSE(observations[0].hasObserverSigma);
	}
}

TEST(AdesParser, InvalidPhotometricSigmaLeavesFallbackAvailable)
{
	for (const std::string value : {"", "0", "-0.1", "nan", "inf"}) {
		std::string error;
		const auto observations = orbital_core::ades::AdesReader{}.read(
		    xml("<mag>19.1</mag><band>V</band><rmsMag>" + value + "</rmsMag>"), error);
		ASSERT_EQ(observations.size(), 1u) << error;
		EXPECT_DOUBLE_EQ(observations[0].rmsMag, 0.0);
	}
}

TEST(AdesParser, Base62ProgramsIncludingNonMpc80Codes)
{
	orbital_core::ades::AdesReader reader;
	for (const auto& [prog, expected] : std::initializer_list<std::pair<std::string, char>>{{"00", '0'},
	                                                                                        {"01", '1'},
	                                                                                        {"05", '5'},
	                                                                                        {"0C", '#'},
	                                                                                        {"0F", '&'},
	                                                                                        {"0W", '|'},
	                                                                                        {"1V", 'z'},
	                                                                                        {"1W", ' '},
	                                                                                        {"zz", ' '}}) {
		std::string error;
		const auto observations = reader.read(xml("<prog>" + prog + "</prog>"), error);
		ASSERT_EQ(observations.size(), 1u) << error;
		EXPECT_EQ(observations[0].programCode, expected) << prog;
		EXPECT_EQ(observations[0].adesProgramCode, prog);
	}
}

TEST(AdesParser, NotesAloneNeverBecomeProgram)
{
	std::string error;
	const auto observations = orbital_core::ades::AdesReader{}.read(xml("<prog/><notes>K&amp;st</notes>"), error);
	ASSERT_EQ(observations.size(), 1u) << error;
	EXPECT_EQ(observations[0].programCode, ' ');
	EXPECT_TRUE(observations[0].adesProgramCode.empty());
	EXPECT_EQ(observations[0].notes, "K&st");
}

TEST(AdesParser, InvalidProgramIsNotSilentlyTruncated)
{
	for (const std::string prog : {"0", "0F0", "?F"}) {
		std::string error;
		EXPECT_TRUE(orbital_core::ades::AdesReader{}.read(xml("<prog>" + prog + "</prog>"), error).empty());
		EXPECT_NE(error.find("prog"), std::string::npos);
	}
}

TEST(AdesParser, ProgramReachesWeightLookupAndReportedRmsSurvives)
{
	const auto dataDir = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	getExecutableDir() = dataDir;
	DE440Source source;
	ASSERT_TRUE(source.isAvailable());
	setDefaultEphemerisSource(&source);

	struct ResetSource {
		~ResetSource() { setDefaultEphemerisSource(nullptr); }
	} reset;

	std::string error;
	const auto observations = orbital_core::ades::AdesReader{}.read(
	    xml("<prog>0F</prog><notes>K</notes><rmsRA>0.12</rmsRA><rmsDec>0.13</rmsDec>"
	        "<mag>21.614</mag><band>G</band><rmsMag>0.087</rmsMag>"),
	    error);
	ASSERT_EQ(observations.size(), 1u) << error;
	const ObsNode node(observations[0], ObservatoryRecord{.code = "309"});
	EXPECT_EQ(node.programCode, '&');
	EXPECT_EQ(node.note1, '&');
	EXPECT_EQ(node.notes, "K");
	EXPECT_EQ(node.adesProgramCode, "0F");
	EXPECT_DOUBLE_EQ(node.rmsRaArcsec, 0.12);
	EXPECT_DOUBLE_EQ(node.rmsDecArcsec, 0.13);
	EXPECT_DOUBLE_EQ(node.rmsMag, 0.087);
	EXPECT_DOUBLE_EQ(node.MAG, 21.614);
	EXPECT_EQ(node.magString, "21.614");
	EXPECT_EQ(node.magBand, 'G');
	const auto uncertainty =
	    orbital_core::obs::lookupUncertainty(node, node.JD, orbital_core::AstrometricWeightProfile::Vfcc17);
	EXPECT_NEAR(uncertainty.sigmaRaRad / Constants::ARCSEC_RAD, 0.2, 1e-14);
	EXPECT_NEAR(uncertainty.sigmaDecRad / Constants::ARCSEC_RAD, 0.2, 1e-14);
}

TEST(AdesParser, AllSupportedModesMapToWeightingTechniques)
{
	for (const auto& [mode, code] : std::initializer_list<std::pair<std::string, char>>{{"CCD", 'C'},
	                                                                                    {"TDI", 'C'},
	                                                                                    {"CMO", 'B'},
	                                                                                    {"VID", 'n'},
	                                                                                    {"PHO", 'P'},
	                                                                                    {"PMT", 'H'},
	                                                                                    {"ENC", 'e'},
	                                                                                    {"MIC", 'M'},
	                                                                                    {"MER", 'T'},
	                                                                                    {"OCC", 'E'}}) {
		std::string text = xml("");
		text.replace(text.find("<mode>CCD</mode>"), 16, "<mode>" + mode + "</mode>");
		std::string error;
		const auto observations = orbital_core::ades::AdesReader{}.read(text, error);
		ASSERT_EQ(observations.size(), 1u) << error;
		EXPECT_EQ(orbital_core::obs::mpcNote2(observations[0].technique), code);
	}
}

TEST(AdesParser, Vv2020MatchesAll79ReferenceSigmas)
{
	const auto root = std::filesystem::path(__FILE__).parent_path();
	getExecutableDir() = root.parent_path().parent_path() / "data";
	std::ifstream input(root / "data/2020VV_ADES_2022.xml");
	std::ifstream expected(root / "data/2020VV_vfcc17_sigmas.txt");
	ASSERT_TRUE(input.good());
	ASSERT_TRUE(expected.good());
	const std::string text{std::istreambuf_iterator<char>(input), {}};
	std::string error;
	const auto observations = orbital_core::ades::AdesReader{}.read(text, error);
	ASSERT_EQ(observations.size(), 79u) << error;
	int programs = 0, reported = 0, station309 = 0;
	for (const auto& data : observations) {
		double sigmaRa = 0.0, sigmaDec = 0.0;
		ASSERT_TRUE(static_cast<bool>(expected >> sigmaRa >> sigmaDec));
		if (!data.adesProgramCode.empty())
			++programs;
		if (data.obsCode == "309") {
			++station309;
			EXPECT_EQ(data.programCode, '&');
		}
		if (data.hasObserverSigma) {
			++reported;
			EXPECT_NEAR(data.rmsRaArcsec, sigmaRa, 1e-14);
			EXPECT_NEAR(data.rmsDecArcsec, sigmaDec, 1e-14);
		} else {
			ObsNode node;
			node.obsCode = data.obsCode;
			node.catalogCode = data.catalogCode;
			node.note2 = orbital_core::obs::mpcNote2(data.technique);
			node.programCode = data.programCode;
			const auto uncertainty =
			    orbital_core::obs::lookupUncertainty(node, data.jd, orbital_core::AstrometricWeightProfile::Vfcc17);
			EXPECT_NEAR(uncertainty.sigmaRaRad / Constants::ARCSEC_RAD, sigmaRa, 1e-14) << data.obsCode;
			EXPECT_NEAR(uncertainty.sigmaDecRad / Constants::ARCSEC_RAD, sigmaDec, 1e-14) << data.obsCode;
		}
	}
	EXPECT_EQ(programs, 28);
	EXPECT_EQ(reported, 18);
	EXPECT_EQ(station309, 6);
}

namespace {
	std::string occultationXml(const std::string& offsets, const std::string& rms)
	{
		return "<ades version=\"2022\"><occultation><permID>99942</permID><trkSub>99942</trkSub>"
		       "<mode>OCC</mode><stn>275</stn><sys>ICRF_KM</sys><ctr>399</ctr><pos1>-1585.585</pos1>"
		       "<pos2>-1216.481</pos2><pos3>3799.302</pos3><obsTime>2021-03-07T05:55:57.56Z</obsTime>"
		       "<raStar>140.344694356</raStar><decStar>-6.671805427</decStar>" +
		       offsets + rms + "<astCat>Gaia3E</astCat><shapeOcc>0</shapeOcc></occultation></ades>";
	}
}  // namespace

TEST(AdesParser, XmlKeepsMixedRecordOrderAndRejectsUnterminatedBlocks)
{
	const auto optical = xml("");
	const std::string radar =
	    "<radar><permID>99942</permID><trx>251</trx><rcv>251</rcv>"
	    "<obsTime>2005-10-30T13:41:59.9712Z</obsTime><com>1</com><frq>2380</frq>"
	    "<delay>192.02850713</delay><rmsDelay>4</rmsDelay></radar>";
	const auto occultation = occultationXml("<deltaRA>0</deltaRA><deltaDec>0</deltaDec>", "");
	std::string error;
	const auto records = orbital_core::ades::AdesReader{}.read(optical + radar + occultation + optical, error);
	ASSERT_TRUE(error.empty()) << error;
	ASSERT_EQ(records.size(), 4u);
	EXPECT_FALSE(records[0].radar);
	ASSERT_TRUE(records[1].radar);
	EXPECT_EQ(records[2].technique, orbital_core::obs::Technique::Occultation);
	EXPECT_FALSE(records[3].radar);
	EXPECT_EQ(records[0].jd, records[3].jd);

	for (const std::string tag : {"optical", "radar", "occultation"}) {
		error.clear();
		EXPECT_TRUE(orbital_core::ades::AdesReader{}.read(optical + "<" + tag + ">", error).empty());
		EXPECT_EQ(error, "ADES XML: unterminated <" + tag + "> block");
	}
}

TEST(AdesParser, OccultationOffsetsConvertConsistentlyInXmlAndPsv)
{
	const std::string psv =
	    "# version=2022\n"
	    "permID|trkSub|mode|stn|obsTime|ra|dec|astCat\n"
	    "99942|K04M04N|CCD|691|2004-03-15T02:35:21.696Z|61.53367|16.91794|USNOB1\n"
	    "permID|trkSub|mode|stn|sys|ctr|pos1|pos2|pos3|obsTime|raStar|decStar|deltaRA|deltaDec|rmsRA|rmsDec|rmsCorr|"
	    "astCat\n"
	    "99942|99942|OCC|275|ICRF_KM|399|-1585.585|-1216.481|3799.302|2021-03-07T05:55:57.56Z|"
	    "140.344694356|-6.671805427|0.0000|0.0000|0.0008|0.0007|-0.755|Gaia3E\n";
	const auto xmlText = occultationXml("<deltaRA>0.0000</deltaRA><deltaDec>0.0000</deltaDec>",
	                                    "<rmsRA>0.0008</rmsRA><rmsDec>0.0007</rmsDec><rmsCorr>-0.755</rmsCorr>");
	for (const auto& [text, rows] : {std::pair{xmlText, 1u}, std::pair{psv, 2u}}) {
		std::string error;
		const auto parsed = orbital_core::ades::AdesReader{}.read(text, error);
		ASSERT_EQ(parsed.size(), rows) << error;
		const auto& occultation = parsed.back();
		EXPECT_EQ(occultation.objectName, "99942");
		EXPECT_EQ(occultation.obsCode, "275");
		EXPECT_EQ(occultation.technique, orbital_core::obs::Technique::Occultation);
		EXPECT_DOUBLE_EQ(occultation.raRad, 140.344694356 * Constants::DEG_TO_RAD);
		EXPECT_DOUBLE_EQ(occultation.decRad, -6.671805427 * Constants::DEG_TO_RAD);
		EXPECT_TRUE(occultation.isSatellite);
		EXPECT_EQ(occultation.offsetUnits, orbital_core::obs::OffsetUnits::Km);
		EXPECT_DOUBLE_EQ(occultation.offsetX, -1585.585);
		EXPECT_DOUBLE_EQ(occultation.offsetY, -1216.481);
		EXPECT_DOUBLE_EQ(occultation.offsetZ, 3799.302);
		EXPECT_TRUE(occultation.hasObserverSigma);
		EXPECT_DOUBLE_EQ(occultation.rmsRaArcsec, 0.0008);
		EXPECT_DOUBLE_EQ(occultation.rmsDecArcsec, 0.0007);
		EXPECT_DOUBLE_EQ(occultation.rmsCorr, -0.755);
		EXPECT_DOUBLE_EQ(occultation.raAccuracyArcsec, 1.0e-4);
		EXPECT_DOUBLE_EQ(occultation.decAccuracyArcsec, 1.0e-4);
		if (rows == 2u)
			EXPECT_EQ(parsed.front().technique, orbital_core::obs::Technique::Ccd);
	}
}

TEST(AdesParser, OccultationOffsetsRespectPrecisionFloor)
{
	std::string error;
	const auto offset = orbital_core::ades::AdesReader{}.read(
	    occultationXml("<deltaRA>0.0123</deltaRA><deltaDec>-0.0456</deltaDec>",
	                   "<rmsRA>0.00005</rmsRA><rmsDec>0.00004</rmsDec><rmsCorr>0.3</rmsCorr>"),
	    error);
	ASSERT_EQ(offset.size(), 1u) << error;
	const double dec = (-6.671805427 - 0.0456 / 3600.0) * Constants::DEG_TO_RAD;
	EXPECT_NEAR(offset[0].decRad, dec, 1e-15);
	EXPECT_NEAR(offset[0].raRad, 140.344694356 * Constants::DEG_TO_RAD + 0.0123 * Constants::ARCSEC_RAD / std::cos(dec),
	            1e-15);
	EXPECT_DOUBLE_EQ(offset[0].rmsRaArcsec, 1.0e-4);
	EXPECT_DOUBLE_EQ(offset[0].rmsDecArcsec, 1.0e-4);
	EXPECT_DOUBLE_EQ(offset[0].rmsCorr, 0.0);

	for (const std::string& invalid :
	     {occultationXml("<deltaRA>0.0</deltaRA>", ""), occultationXml("<dist>0.01</dist><pa>45</pa>", ""),
	      occultationXml("<deltaRA>0.0</deltaRA><deltaDec>0.0</deltaDec><ra>1.0</ra>", "")}) {
		EXPECT_TRUE(orbital_core::ades::AdesReader{}.read(invalid, error).empty());
		EXPECT_FALSE(error.empty());
	}
}

TEST(AdesParser, DeprecatedRecordIsReplacedAsInTheEightyColumnRecord)
{
	std::string error;
	const auto parsed = orbital_core::ades::AdesReader{}.read(xml("<deprecated>X</deprecated>"), error);
	ASSERT_EQ(parsed.size(), 1u) << error;
	EXPECT_EQ(orbital_core::obs::mpcNote2(parsed[0].technique), 'X');
}

TEST(AdesParser, OccultationKeepsItsTechniqueAndOwnSigma)
{
	getExecutableDir() = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	DE440Source source;
	ASSERT_TRUE(source.isAvailable());
	setDefaultEphemerisSource(&source);

	struct ResetSource {
		~ResetSource() { setDefaultEphemerisSource(nullptr); }
	} reset;

	const auto node = [](const std::string& text, const char* station) {
		std::string error;
		const auto parsed = orbital_core::ades::AdesReader{}.read(text, error);
		EXPECT_EQ(parsed.size(), 1u) << error;
		return ObsNode(parsed.at(0), ObservatoryRecord{.code = station});
	};
	const auto occultation = node(occultationXml("<deltaRA>0.0000</deltaRA><deltaDec>0.0000</deltaDec>",
	                                             "<rmsRA>0.0008</rmsRA><rmsDec>0.0007</rmsDec>"),
	                              "275");
	EXPECT_EQ(occultation.note2, 'E');
	EXPECT_TRUE(occultation.isSatellite);
	EXPECT_TRUE(orbital_core::obs::usesObserverSigma(occultation, false));

	const auto satellite = node(xml("<sys>ICRF_KM</sys><ctr>399</ctr><pos1>1</pos1><pos2>2</pos2><pos3>3</pos3>"
	                                "<rmsRA>0.1</rmsRA><rmsDec>0.1</rmsDec>"),
	                            "309");
	EXPECT_EQ(satellite.note2, 'S');
	EXPECT_FALSE(orbital_core::obs::usesObserverSigma(satellite, false));
	EXPECT_TRUE(orbital_core::obs::usesObserverSigma(satellite, true));

	const auto deprecated = node(xml("<deprecated>X</deprecated><rmsRA>0.1</rmsRA><rmsDec>0.1</rmsDec>"), "309");
	EXPECT_EQ(deprecated.note2, 'X');
	EXPECT_FALSE(orbital_core::obs::usesObserverSigma(deprecated, true));
}
