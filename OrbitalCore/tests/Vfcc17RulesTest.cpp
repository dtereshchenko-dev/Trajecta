#include <gtest/gtest.h>

#include "obs/Vfcc17Rules.h"
#include "obs/AstrometricSigma.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/FitConfig.h"
#include "orbital_core/ObservationSet.h"
#include "orbital_core/Types.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {
	using orbital_core::obs::Vfcc17Rules;

	std::string rule(const char* station, const char* tech, const char* cat, const char* prog, double ra, double dec,
	                 const char* start = "", const char* end = "")
	{
		std::ostringstream os;
		os << std::left << std::setw(3) << station << " t=" << std::setw(5) << tech << " c=" << std::setw(10) << cat
		   << " p=" << std::setw(3) << prog << " > " << std::setw(10) << start << " < " << std::setw(10) << end << " @ "
		   << std::right << std::fixed << std::setprecision(2) << std::setw(5) << ra << ", " << std::setw(5) << dec
		   << '\n';
		return os.str();
	}

	double jd(int year, int month = 1, double day = 1.0) { return Convertors::convertGregorianToJD(year, month, day); }

	double ra(const Vfcc17Rules& rules, const char* station, char tech = 'C', char cat = 'V', char program = ' ',
	          double epoch = jd(2020))
	{
		return rules.lookup(station, tech, cat, program, epoch).sigmaRaRad / Constants::ARCSEC_RAD;
	}

	class Vfcc17RulesTest : public ::testing::Test {
	   protected:
		void SetUp() override
		{
			std::ifstream stream(std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() /
			                     "data/weights/vfcc17.rules");
			ASSERT_TRUE(stream.good());
			const std::string text{std::istreambuf_iterator<char>(stream), {}};
			ASSERT_TRUE(rules.load(text, error)) << error;
		}

		Vfcc17Rules rules;
		std::string error;
	};
}  // namespace

TEST(AstrometricWeightProfileTest, DefaultsToVfcc17)
{
	using namespace orbital_core;
	EXPECT_EQ(FitConfig{}.weightProfile, AstrometricWeightProfile::Vfcc17);
	EXPECT_EQ(ObservationPreparationOptions{}.weightProfile, AstrometricWeightProfile::Vfcc17);
	EXPECT_EQ(static_cast<int>(AstrometricWeightProfile::Vfcc17), 1);
	EXPECT_EQ(static_cast<int>(AstrometricWeightProfile::Uniform), 2);
}

TEST(AstrometricWeightProfileTest, UniformUsesHalfArcsecondWithoutSigmaFile)
{
	ObsNode observation;
	observation.obsCode = "999";
	const auto sigma =
	    orbital_core::obs::lookupUncertainty(observation, jd(2020), orbital_core::AstrometricWeightProfile::Uniform);
	EXPECT_TRUE(sigma.resolved);
	EXPECT_DOUBLE_EQ(sigma.sigmaRaRad, 0.5 * Constants::ARCSEC_RAD);
	EXPECT_DOUBLE_EQ(sigma.sigmaDecRad, 0.5 * Constants::ARCSEC_RAD);
}

TEST_F(Vfcc17RulesTest, LoadsEveryProfileRule)
{
	EXPECT_EQ(rules.size(), 100u);
	EXPECT_DOUBLE_EQ(ra(rules, "999", 'C', 'V'), 0.6);
	EXPECT_DOUBLE_EQ(ra(rules, "999", 'C', ' '), 1.0);
	EXPECT_DOUBLE_EQ(ra(rules, "999", 'B', ' '), 1.0);
	EXPECT_DOUBLE_EQ(ra(rules, "595", 'B', 'X'), 1.0);
	EXPECT_DOUBLE_EQ(rules.lookup("595", 'B', 'X', ' ', jd(2023, 11, 7.93)).sigmaDecRad / Constants::ARCSEC_RAD, 1.0);
	EXPECT_DOUBLE_EQ(ra(rules, "F51"), 0.2);
	EXPECT_DOUBLE_EQ(ra(rules, "G96"), 0.5);
}

TEST_F(Vfcc17RulesTest, TechniquesAndHistoricalDateBoundaries)
{
	for (char t : std::string("PAN")) {
		EXPECT_DOUBLE_EQ(ra(rules, "999", t, ' ', ' ', jd(1880)), 10.0);
		EXPECT_DOUBLE_EQ(ra(rules, "999", t, ' ', ' ', jd(1890)), 10.0);
		EXPECT_DOUBLE_EQ(ra(rules, "999", t, ' ', ' ', jd(1890) + 0.001), 5.0);
		EXPECT_DOUBLE_EQ(ra(rules, "999", t, ' ', ' ', jd(1950)), 5.0);
		EXPECT_DOUBLE_EQ(ra(rules, "999", t, ' ', ' ', jd(1950) + 0.001), 2.5);
	}
	for (const auto& [tech, sigma] : std::initializer_list<std::pair<char, double>>{{'E', 0.2},
	                                                                                {'H', 0.4},
	                                                                                {'T', 0.5},
	                                                                                {'e', 0.75},
	                                                                                {'M', 3.0},
	                                                                                {'S', 1.5},
	                                                                                {'B', 1.0},
	                                                                                {'V', 1.0},
	                                                                                {'n', 1.0},
	                                                                                {'x', 999.0},
	                                                                                {'X', 999.0}}) {
		EXPECT_DOUBLE_EQ(ra(rules, "999", tech), sigma) << tech;
	}
}

TEST_F(Vfcc17RulesTest, CatalogAndProgramSets)
{
	for (char tech : std::string("cC"))
		for (char catalog : std::string("UVXW")) {
			EXPECT_DOUBLE_EQ(ra(rules, "309", tech, catalog, '&'), 0.2);
			EXPECT_DOUBLE_EQ(ra(rules, "309", tech, catalog, '%'), 0.2);
			EXPECT_DOUBLE_EQ(ra(rules, "309", tech, catalog, ' '), 0.6);
		}
	EXPECT_DOUBLE_EQ(ra(rules, "309", 'C', 't', '&'), 0.3);
	EXPECT_DOUBLE_EQ(ra(rules, "568", 'C', 'V', '2'), 0.1);
	EXPECT_DOUBLE_EQ(ra(rules, "568", 'C', 'o', '2'), 0.5);
	EXPECT_DOUBLE_EQ(ra(rules, "568", 'C', 'V', '^'), 0.2);
	EXPECT_DOUBLE_EQ(ra(rules, "568", 'C', 'V'), 0.5);
	EXPECT_DOUBLE_EQ(ra(rules, "T09", 'C', 'V', '0'), 0.1);
	EXPECT_DOUBLE_EQ(ra(rules, "T14", 'C', 'V', '7'), 0.1);
	EXPECT_DOUBLE_EQ(ra(rules, "G83", 'C', 'V', '2'), 0.2);
	EXPECT_DOUBLE_EQ(ra(rules, "Z18", 'C', 'V', '1'), 0.2);
	EXPECT_DOUBLE_EQ(ra(rules, "C65", 'C', 'V', '3'), 0.2);
	EXPECT_DOUBLE_EQ(ra(rules, "094", 'C', 'V', '9'), 0.5);
}

TEST_F(Vfcc17RulesTest, StationProgramDateBoundary)
{
	const double boundary = jd(2022, 10, 6);
	EXPECT_DOUBLE_EQ(ra(rules, "T12", 'C', 'V', '1', boundary), 0.1);
	EXPECT_DOUBLE_EQ(ra(rules, "T12", 'C', 'V', '1', boundary + 0.001), 0.2);
	EXPECT_DOUBLE_EQ(ra(rules, "T12", 'C', 'V', ' ', boundary + 0.001), 0.5);
}

TEST(Vfcc17RuleParser, PriorityAndFirstRowWins)
{
	Vfcc17Rules rules;
	std::string error;
	ASSERT_TRUE(rules.load(rule("ALL", "C", "V", "&", 0.9, 0.8) + rule("309", "C", "*", "", 0.7, 0.7) +
	                           rule("309", "C", "*", "&", 0.4, 0.4) + rule("309", "C", "V", "", 0.3, 0.3) +
	                           rule("309", "C", "V", "&", 0.2, 0.1) + rule("309", "C", "V", "&", 0.6, 0.6),
	                       error))
	    << error;
	EXPECT_DOUBLE_EQ(ra(rules, "309", 'C', 'V', '&'), 0.2);
	EXPECT_DOUBLE_EQ(ra(rules, "309", 'C', 'V', '2'), 0.3);
	EXPECT_DOUBLE_EQ(ra(rules, "309", 'C', 't', '&'), 0.4);
	EXPECT_DOUBLE_EQ(ra(rules, "309", 'C', 't', '2'), 0.7);
	EXPECT_DOUBLE_EQ(ra(rules, "999", 'C', 'V', '&'), 0.9);
	EXPECT_DOUBLE_EQ(rules.lookup("309", 'C', 'V', '&', jd(2020)).sigmaDecRad / Constants::ARCSEC_RAD, 0.1);
}

TEST(Vfcc17RuleParser, ProgramPunctuationIsNotCommentOrWildcard)
{
	for (char program : std::string("!#|*^<>@")) {
		Vfcc17Rules rules;
		std::string error;
		const std::string p(1, program);
		ASSERT_TRUE(rules.load("! comment\r\n# comment\r\n" + rule("309", "C", "V", p.c_str(), 0.2, 0.2), error))
		    << error;
		EXPECT_DOUBLE_EQ(ra(rules, "309", 'C', 'V', program), 0.2);
		EXPECT_DOUBLE_EQ(ra(rules, "309", 'C', 'V', ' '), 1.0);
	}
}

TEST(Vfcc17RuleParser, TechniqueWildcardAndAllFallback)
{
	Vfcc17Rules rules;
	std::string error;
	ASSERT_TRUE(rules.load(rule("309", "*", "*", "&", 0.3, 0.3) + rule("ALL", "*", "*", "", 0.7, 0.7), error)) << error;
	EXPECT_DOUBLE_EQ(ra(rules, "309", 'E', 'q', '&'), 0.3);
	EXPECT_DOUBLE_EQ(ra(rules, "999", 'E', 'q', ' '), 0.7);
}

TEST(Vfcc17RuleParser, MalformedFileIsNotPartiallyApplied)
{
	const auto good = rule("ALL", "C", "*", "", 0.5, 0.5);
	for (const auto& bad : {std::string("309 * 0.2 0.2\n"), rule("309", "C", "V", "", -0.2, 0.2),
	                        rule("309", "C", "V", "", 0.2, 0.2, "2021-02-29"),
	                        rule("309", "C", "V", "", 0.2, 0.2, "2022-01-01", "2020-01-01")}) {
		Vfcc17Rules rules;
		std::string error;
		ASSERT_TRUE(rules.load(good, error));
		EXPECT_FALSE(rules.load(good + bad, error));
		EXPECT_EQ(rules.size(), 0u);
		EXPECT_NE(error.find("line 2"), std::string::npos);
		EXPECT_DOUBLE_EQ(ra(rules, "309"), 1.0);
	}
}
