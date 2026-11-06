#include <gtest/gtest.h>

#include "orbital_core/Photometry.h"
#include "orbital_core/FitResult.h"
#include "AdesParser.h"
#include "obs/Mpc80ColumnReader.h"

#include <array>
#include <cmath>
#include <limits>

namespace photo = orbital_core::photometry;

namespace {
	Residual magnitude(double mag, char band = 'V', double sigma = 0.7)
	{
		Residual row{};
		row.mag = mag;
		row.magBand = band;
		row.magSigma = sigma;
		row.helioDistAu = row.geoDistAu = 1.0;
		return row;
	}
}  // namespace

TEST(Photometry, ResolvesReportedAndPrecisionSigmasInOrder)
{
	EXPECT_DOUBLE_EQ(photo::resolveSigma("21.614", 0.087), 0.087);
	EXPECT_DOUBLE_EQ(photo::resolveSigma("21.614", 0.0), 0.5);
	EXPECT_DOUBLE_EQ(photo::resolveSigma("21"), 1.0);
	EXPECT_DOUBLE_EQ(photo::resolveSigma("21.6"), 0.7);
	EXPECT_DOUBLE_EQ(photo::resolveSigma(" 21.60 "), 0.5);
	EXPECT_DOUBLE_EQ(photo::resolveSigma("21.614"), 0.5);
	EXPECT_DOUBLE_EQ(photo::resolveSigma(""), 0.0);
	EXPECT_DOUBLE_EQ(photo::resolveSigma("21.6", -1.0), 0.7);
	EXPECT_DOUBLE_EQ(photo::resolveSigma("21.6", std::numeric_limits<double>::quiet_NaN()), 0.7);
}

TEST(Photometry, FormatsBandSeparatelyWithoutDiscardingIntegerDigits)
{
	EXPECT_EQ(photo::formatMagnitude(21.614, 'G', "21.614"), "21.61G");
	EXPECT_EQ(photo::formatMagnitude(19.1, 'V', "19.1"), "19.1 V");
	EXPECT_EQ(photo::formatMagnitude(9.1, 'V', "9.1"), " 9.1 V");
	EXPECT_EQ(photo::formatMagnitude(19.0, 'R', "19"), "19   R");
	EXPECT_EQ(photo::formatMagnitude(21.60, ' ', "21.60"), "21.60 ");
	EXPECT_EQ(photo::formatMagnitude(0.0, 'V', ""), "      ");
}

TEST(Photometry, HAndResidualsUseOneSigmaWeightedSolution)
{
	const std::array rows{magnitude(20.0, 'G', 0.22), magnitude(21.0, 'r', 0.7)};
	const double expected = (20.28 / 0.22 + 21.23 / 0.7) / (1.0 / 0.22 + 1.0 / 0.7);
	const auto solution = photo::fit(rows);
	ASSERT_TRUE(solution.absoluteMagnitude);
	EXPECT_NEAR(*solution.absoluteMagnitude, expected, 1e-13);
	double weightedResidual = 0.0;
	for (const auto& row : rows) {
		const auto result = photo::evaluate(row, solution);
		ASSERT_TRUE(result.residual);
		EXPECT_TRUE(result.selected);
		EXPECT_DOUBLE_EQ(result.sigma, row.magSigma);
		weightedResidual += *result.residual / result.sigma;
	}
	EXPECT_NEAR(weightedResidual, 0.0, 1e-12);
}

TEST(Photometry, UnsupportedBandsAreNotSilentlyInterpretedAsV)
{
	for (const char band : {'?', 'W', 'v'}) {
		SCOPED_TRACE(band);
		const std::array rows{magnitude(20.0), magnitude(25.0, band)};
		const auto solution = photo::fit(rows);
		ASSERT_TRUE(solution.absoluteMagnitude);
		EXPECT_DOUBLE_EQ(*solution.absoluteMagnitude, 20.0);
		EXPECT_EQ(solution.count, 1u);
		EXPECT_FALSE(photo::evaluate(rows.back(), solution).residual);
	}
}

TEST(Photometry, HgPhaseFunctionIsUnchanged)
{
	auto row = magnitude(21.0, 'V');
	row.helioDistAu = 1.2;
	row.geoDistAu = 0.04;
	row.phaseAngleRad = 0.6;
	const double t = std::tan(row.phaseAngleRad / 2.0);
	const double expected =
	    row.mag - 5.0 * std::log10(1.2 * 0.04) +
	    2.5 * std::log10(0.85 * std::exp(-3.33 * std::pow(t, 0.63)) + 0.15 * std::exp(-1.87 * std::pow(t, 1.22)));
	const auto solution = photo::fit(std::span<const Residual>(&row, 1));
	ASSERT_TRUE(solution.absoluteMagnitude);
	EXPECT_NEAR(*solution.absoluteMagnitude, expected, 1e-13);
}

TEST(Photometry, MissingAndInvalidSamplesDoNotBiasH)
{
	auto invalid = magnitude(21.0);
	invalid.geoDistAu = std::numeric_limits<double>::quiet_NaN();
	const std::array rows{magnitude(20.0), magnitude(0.0), invalid};
	const auto solution = photo::fit(rows);
	ASSERT_TRUE(solution.absoluteMagnitude);
	EXPECT_DOUBLE_EQ(*solution.absoluteMagnitude, 20.0);
	EXPECT_EQ(solution.count, 1u);
	EXPECT_FALSE(photo::evaluate(invalid, solution).residual);
	EXPECT_FALSE(photo::fit({}).absoluteMagnitude);
}

TEST(Photometry, LeaveOneOutFlagsDoNotSilentlyChangeTheMean)
{
	const std::array rows{magnitude(20.0), magnitude(20.0), magnitude(20.0), magnitude(20.0), magnitude(25.0)};
	const auto solution = photo::fit(rows);
	ASSERT_TRUE(solution.absoluteMagnitude);
	EXPECT_NEAR(*solution.absoluteMagnitude, 21.0, 1e-13);
	EXPECT_FALSE(photo::evaluate(rows.back(), solution).selected);
	EXPECT_TRUE(photo::evaluate(rows.front(), solution).selected);
}

TEST(Photometry, Mpc80AndAdesWithoutRmsMagHaveTheSameFallback)
{
	const std::string mpc = "     K14A00A* C2014 01 01.26257 05 32 35.55 +13 59 45.0          19.1 Vq~0yn5G96";
	std::string error;
	const auto legacy = orbital_core::obs::Mpc80ColumnReader{}.read(mpc, error);
	ASSERT_EQ(legacy.size(), 1u) << error;
	const auto ades = orbital_core::ades::AdesReader{}.read(
	    "# version=2022\nprovID|stn|mode|obsTime|ra|dec|mag|band\n"
	    "2014 AA|G96|CCD|2014-01-01T06:18:06.048Z|83.148125|13.995833|19.1|V\n",
	    error);
	ASSERT_EQ(ades.size(), 1u) << error;
	EXPECT_EQ(legacy[0].magnitudeText, ades[0].magnitudeText);
	EXPECT_EQ(legacy[0].magBand, ades[0].magBand);
	EXPECT_DOUBLE_EQ(legacy[0].rmsMag, 0.0);
	EXPECT_DOUBLE_EQ(ades[0].rmsMag, 0.0);
	EXPECT_DOUBLE_EQ(photo::resolveSigma(legacy[0].magnitudeText), 0.7);
	EXPECT_DOUBLE_EQ(photo::resolveSigma(ades[0].magnitudeText), 0.7);
}

TEST(Photometry, CometMagnitudesUseSeparateTypesAndFixedSlope)
{
	auto total = magnitude(12.0, 'T', 0.5);
	total.nuclearMagnitude = false;
	total.helioDistAu = 2.0;
	total.geoDistAu = 0.5;
	auto nuclear = magnitude(16.0, 'N', 1.0);
	nuclear.nuclearMagnitude = true;
	const std::array rows{total, nuclear, magnitude(9.0)};
	const auto result = photo::fitComet(rows);
	ASSERT_TRUE(result.total.absoluteMagnitude);
	ASSERT_TRUE(result.nuclear.absoluteMagnitude);
	EXPECT_NEAR(*result.total.absoluteMagnitude, 12.0 - 5.0 * std::log10(0.5) - 10.0 * std::log10(2.0), 1e-12);
	EXPECT_DOUBLE_EQ(*result.nuclear.absoluteMagnitude, 16.0);
	EXPECT_EQ(result.total.count, 1u);
	EXPECT_EQ(result.nuclear.count, 1u);
	EXPECT_EQ(result.unknownTypeCount, 1u);
	const auto predicted = photo::cometMagnitude(*result.total.absoluteMagnitude, result.total.slope, 2.0, 0.5);
	ASSERT_TRUE(predicted);
	EXPECT_NEAR(*predicted, total.mag, 1e-12);
	const auto evaluation = photo::evaluateComet(total, result);
	ASSERT_TRUE(evaluation.residual);
	EXPECT_NEAR(*evaluation.residual, 0.0, 1e-12);
	EXPECT_FALSE(photo::evaluateComet(rows.back(), result).residual);
}

TEST(Photometry, CometMagnitudeTypeIsReadFromMpcAndAdes)
{
	std::string mpc = "     K14A00A* C2014 01 01.26257 05 32 35.55 +13 59 45.0          19.1 Vq~0yn5G96";
	mpc[70] = 'N';
	std::string error;
	const auto legacy = orbital_core::obs::Mpc80ColumnReader{}.read(mpc, error);
	ASSERT_EQ(legacy.size(), 1u) << error;
	ASSERT_TRUE(legacy[0].nuclearMagnitude.has_value());
	EXPECT_TRUE(*legacy[0].nuclearMagnitude);
	const auto ades = orbital_core::ades::AdesReader{}.read(
	    "# version=2022\nprovID|stn|mode|obsTime|ra|dec|mag|band|nucMag\n"
	    "2014 AA|G96|CCD|2014-01-01T06:18:06.048Z|83.148125|13.995833|19.1|Pr|0\n",
	    error);
	ASSERT_EQ(ades.size(), 1u) << error;
	EXPECT_FALSE(ades[0].nuclearMagnitude.has_value());
	EXPECT_EQ(ades[0].magBand, 'r');
}

TEST(Photometry, NonDefaultCoefficientsAreUsedInReductionAndPrediction)
{
	orbital_core::PhotometryOptions options;
	options.asteroidG = 0.42;
	options.cometK1 = 8.0;
	options.cometK2 = 5.0;
	options.cometNuclearPhaseCoefficient = 0.03;
	auto total = magnitude(0.0, 'T');
	total.nuclearMagnitude = false;
	total.helioDistAu = 2.0;
	total.geoDistAu = 0.5;
	total.phaseAngleRad = std::acos(-1.0) / 3.0;
	total.mag = 12.0 + 5.0 * std::log10(0.5) + 8.0 * std::log10(2.0);
	auto nuclear = total;
	nuclear.magBand = 'N';
	nuclear.nuclearMagnitude = true;
	nuclear.mag = 16.0 + 5.0 * std::log10(0.5) + 5.0 * std::log10(2.0) + 0.03 * 60.0;
	const std::array rows{total, nuclear};
	const auto solution = photo::fitComet(rows, options);
	ASSERT_TRUE(solution.total.absoluteMagnitude);
	ASSERT_TRUE(solution.nuclear.absoluteMagnitude);
	EXPECT_NEAR(*solution.total.absoluteMagnitude, 12.0, 1e-12);
	EXPECT_NEAR(*solution.nuclear.absoluteMagnitude, 16.0, 1e-12);
	ASSERT_TRUE(photo::evaluateComet(nuclear, solution).residual);
	EXPECT_NEAR(*photo::evaluateComet(nuclear, solution).residual, 0.0, 1e-12);
	auto asteroid = magnitude(20.0);
	asteroid.phaseAngleRad = 0.6;
	const auto h = photo::fit(std::span<const Residual>(&asteroid, 1), options.asteroidG);
	ASSERT_TRUE(h.absoluteMagnitude);
	const auto prediction =
	    photo::visualMagnitude(*h.absoluteMagnitude, options.asteroidG, 1.0, 1.0, asteroid.phaseAngleRad);
	ASSERT_TRUE(prediction);
	EXPECT_NEAR(*prediction, asteroid.mag, 1e-12);
}
