#include <gtest/gtest.h>
#include "AstrometryBias.h"
#include "orbital_core/Types.h"
#include "orbital_core/ObsNode.h"
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace {

	constexpr int kExpectedNpix = 49152;
	constexpr int kRowColumns = 104;
	constexpr double kBiasReferenceJd = 2'451'544.5;

	std::filesystem::path findBiasDat()
	{
		namespace fs = std::filesystem;
		const fs::path srcDir = fs::path{__FILE__}.parent_path().parent_path();
		const fs::path candidates[] = {
		    srcDir.parent_path() / "data" / "debias" / "efcc18.dat",
		    srcDir / "data" / "debias" / "efcc18.dat",
		    fs::current_path() / "debias" / "efcc18.dat",
		    srcDir.parent_path() / "data" / "debias_2018" / "bias.dat",
		    srcDir / "data" / "debias_2018" / "bias.dat",
		    srcDir / "dependencies" / "debias_2018" / "bias.dat",
		    fs::current_path() / ".." / "data" / "debias_2018" / "bias.dat",
		    fs::current_path() / ".." / "dependencies" / "debias_2018" / "bias.dat",
		};

		for (const auto& p : candidates) {
			if (fs::exists(p))
				return p;
		}

		return {};
	}

	std::filesystem::path findTilesDat()
	{
		namespace fs = std::filesystem;
		const fs::path srcDir = fs::path{__FILE__}.parent_path().parent_path();
		const fs::path candidates[] = {
		    srcDir.parent_path() / "data" / "debias" / "efcc18_tiles.dat",
		    srcDir / "data" / "debias" / "efcc18_tiles.dat",
		    srcDir.parent_path() / "data" / "debias_2018" / "tiles.dat",
		    srcDir / "data" / "debias_2018" / "tiles.dat",
		    srcDir / "dependencies" / "debias_2018" / "tiles.dat",
		    fs::current_path() / ".." / "data" / "debias_2018" / "tiles.dat",
		    fs::current_path() / "dependencies" / "debias_2018" / "tiles.dat",
		};

		for (const auto& p : candidates)
			if (fs::exists(p))
				return p;

		return {};
	}

	std::filesystem::path findFcct14Dat()
	{
		namespace fs = std::filesystem;
		const fs::path srcDir = fs::path{__FILE__}.parent_path().parent_path();
		const fs::path candidates[] = {
		    srcDir.parent_path() / "data" / "debias" / "fcct14.dat",
		    srcDir / "data" / "debias" / "fcct14.dat",
		    fs::current_path() / "debias" / "fcct14.dat",
		    srcDir.parent_path() / "data" / "debias_fcct14" / "xcat.bias_fcct14",
		    srcDir / "data" / "debias_fcct14" / "xcat.bias_fcct14",
		    fs::current_path() / ".." / "data" / "debias_fcct14" / "xcat.bias_fcct14",
		};

		for (const auto& p : candidates)
			if (fs::exists(p))
				return p;
		return {};
	}

	double hmsToRadians(int hour, int minute, double second)
	{
		return (hour + minute / 60.0 + second / 3600.0) * Constants::PI / 12.0;
	}

	double dmsToRadians(int degree, int minute, double second)
	{
		const double sign = degree < 0 ? -1.0 : 1.0;
		return sign * (std::abs(degree) + minute / 60.0 + second / 3600.0) * Constants::PI / 180.0;
	}
}  // namespace

TEST(AstrometryBiasPaths, PrefersCurrentFileAndFallsBackToLegacy)
{
	namespace fs = std::filesystem;
	const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
	const auto base = fs::temp_directory_path() / ("trajecta-debias-" + std::to_string(stamp));
	const auto currentDir = base / "debias";
	const auto legacyDir = base / "debias_2018";
	const auto current = currentDir / "efcc18.dat";
	const auto legacy = legacyDir / "bias.dat";

	struct Cleanup {
		fs::path current, legacy, currentDir, legacyDir, base;

		~Cleanup()
		{
			std::error_code ec;
			fs::remove(current, ec);
			fs::remove(legacy, ec);
			fs::remove(currentDir, ec);
			fs::remove(legacyDir, ec);
			fs::remove(base, ec);
		}
	} cleanup{current, legacy, currentDir, legacyDir, base};

	ASSERT_TRUE(fs::create_directories(currentDir));
	ASSERT_TRUE(fs::create_directories(legacyDir));
	const auto pathFor = [&](orbital_core::AstrometricWeightProfile profile) {
		return Calculations::astrometryBiasPathForProfile(base, profile);
	};
	EXPECT_EQ(pathFor(orbital_core::AstrometricWeightProfile::Vfcc17), current);
	{
		std::ofstream file(legacy);
		ASSERT_TRUE(file);
	}
	EXPECT_EQ(pathFor(orbital_core::AstrometricWeightProfile::Vfcc17), legacy);
	{
		std::ofstream file(current);
		ASSERT_TRUE(file);
	}
	EXPECT_EQ(pathFor(orbital_core::AstrometricWeightProfile::Vfcc17), current);
	EXPECT_EQ(pathFor(orbital_core::AstrometricWeightProfile::Uniform), current);
}

TEST(HEALPix, NorthPoleAndSouthPole)
{
	const int northPix = AstrometryBias::ang2pixRing(1e-10, 0.0);
	EXPECT_GE(northPix, 0);
	EXPECT_LT(northPix, 4);

	const int southPix = AstrometryBias::ang2pixRing(Constants::PI - 1e-10, 0.0);
	EXPECT_GE(southPix, kExpectedNpix - 4);
	EXPECT_LT(southPix, kExpectedNpix);
}

TEST(HEALPix, PixelsInRange)
{
	for (double theta = 0.01; theta < Constants::PI; theta += 0.1) {
		for (double phi = 0.0; phi < 2 * Constants::PI; phi += 0.25) {
			const int p = AstrometryBias::ang2pixRing(theta, phi);
			EXPECT_GE(p, 0);
			EXPECT_LT(p, kExpectedNpix);
		}
	}
}

TEST(HEALPix, PhiWrapsModulo2Pi)
{
	for (double theta = 0.3; theta < Constants::PI; theta += 0.5) {
		for (double phi = 0.1; phi < 2 * Constants::PI; phi += 0.7) {
			EXPECT_EQ(AstrometryBias::ang2pixRing(theta, phi),
			          AstrometryBias::ang2pixRing(theta, phi + 2 * Constants::PI));
			EXPECT_EQ(AstrometryBias::ang2pixRing(theta, phi),
			          AstrometryBias::ang2pixRing(theta, phi - 2 * Constants::PI));
		}
	}
}

TEST(HEALPix, RoundTripAllPixelsViaTilesDat)
{
	const auto tilesPath = findTilesDat();
	if (tilesPath.empty())
		GTEST_SKIP() << "tiles.dat not found";

	std::ifstream in(tilesPath);
	ASSERT_TRUE(in) << "cannot open " << tilesPath;

	int checked = 0;
	int mismatches = 0;
	std::string line;
	while (std::getline(in, line)) {
		const auto first = line.find_first_not_of(" \t\r\n");
		if (first == std::string::npos)
			continue;
		if (line[first] == '!')
			continue;

		std::istringstream iss(line);
		int pixExpected = -1;
		double raRad = 0.0, decRad = 0.0;
		iss >> pixExpected >> raRad >> decRad;
		if (!iss)
			continue;

		const double theta = Constants::PI * 0.5 - decRad;
		const int pixActual = AstrometryBias::ang2pixRing(theta, raRad);
		if (pixActual != pixExpected)
			++mismatches;
		++checked;
	}

	EXPECT_EQ(checked, kExpectedNpix) << "expected 49152 tile rows";
	EXPECT_EQ(mismatches, 0) << "HEALPix ang2pix disagrees with tiles.dat centres";
}

TEST(AstrometryBias, LoadAndCoverMainCatalogs)
{
	const auto p = findBiasDat();
	if (p.empty())
		GTEST_SKIP() << "EFCC18 table not found";

	AstrometryBias bias;
	ASSERT_NO_THROW(bias.load(p));
	EXPECT_TRUE(bias.loaded());

	for (char c : {'q', 'L', 't', 'o', 'Q', 'U', 'v', 'm', 'u', 'r', 'c', 'd'}) {
		const auto corr = bias.correction(Constants::PI, 0.1, c, 2458736.0);
		EXPECT_TRUE(std::isfinite(corr.dRA_arcsec));
		EXPECT_TRUE(std::isfinite(corr.dDEC_arcsec));
	}
}

TEST(AstrometryBias, UnknownCatalogReturnsZero)
{
	const auto p = findBiasDat();
	if (p.empty())
		GTEST_SKIP() << "EFCC18 table not found";

	AstrometryBias bias;
	bias.load(p);

	for (char c : {'V', 'W', ' ', '\0', 'Z'}) {
		const auto corr = bias.correction(1.0, 0.5, c, 2458700.0);
		EXPECT_DOUBLE_EQ(corr.dRA_arcsec, 0.0) << "catalog='" << c << "'";
		EXPECT_DOUBLE_EQ(corr.dDEC_arcsec, 0.0) << "catalog='" << c << "'";
	}
}

TEST(AstrometryBias, CorrectionHasExpectedMagnitude)
{
	const auto p = findBiasDat();
	if (p.empty())
		GTEST_SKIP() << "EFCC18 table not found";

	AstrometryBias bias;
	bias.load(p);

	const double raRad = 10.0 * Constants::PI / 12.0;
	const double decRad = 0.0;
	const auto corr = bias.correction(raRad, decRad, 'q', 2458800.0);

	EXPECT_LT(std::fabs(corr.dRA_arcsec), 2.0) << "RA bias unreasonably large";
	EXPECT_LT(std::fabs(corr.dDEC_arcsec), 2.0) << "DEC bias unreasonably large";
}

TEST(AstrometryBias, LoadParsesExactValuesForKnownPixels)
{
	const auto p = findBiasDat();
	if (p.empty())
		GTEST_SKIP() << "EFCC18 table not found";

	AstrometryBias bias;
	ASSERT_NO_THROW(bias.load(p));

	const double ra = 3.129321;
	const double dec = 0.0;

	const auto cA = bias.correction(ra, dec, 'a', kBiasReferenceJd);
	EXPECT_NEAR(cA.dRA_arcsec, -0.591, 1e-5);
	EXPECT_NEAR(cA.dDEC_arcsec, -0.068, 1e-5);

	const auto cB = bias.correction(ra, dec, 'b', kBiasReferenceJd);
	EXPECT_NEAR(cB.dRA_arcsec, -0.578, 1e-5);
	EXPECT_NEAR(cB.dDEC_arcsec, -0.057, 1e-5);

	const auto cY = bias.correction(ra, dec, 'Y', kBiasReferenceJd);
	EXPECT_NEAR(cY.dRA_arcsec, 0.002, 1e-5);
	EXPECT_NEAR(cY.dDEC_arcsec, 0.003, 1e-5);
}

TEST(AstrometryBias, Fcct14MatchesPublishedRwoBiases)
{
	const auto p = findFcct14Dat();
	if (p.empty())
		GTEST_SKIP() << "FCCT14 table not found";

	AstrometryBias bias;
	ASSERT_NO_THROW(bias.load(p));

	struct Reference {
		double jd;
		double ra;
		double dec;
		char catalog;
		double raCosDecArcsec;
		double decArcsec;
	};

	const Reference references[] = {
	    {2459159.594130, hmsToRadians(0, 49, 57.680), dmsToRadians(9, 21, 46.20), 'u', -0.027, 0.133},
	    {2459167.614300, hmsToRadians(1, 31, 50.950), dmsToRadians(13, 21, 58.90), 'u', -0.095, -0.099},
	    {2459163.755940, hmsToRadians(1, 14, 21.230), dmsToRadians(11, 47, 23.90), 'o', -0.312, 0.391},
	};

	for (const auto& ref : references) {
		const auto correction = bias.correction(ref.ra, ref.dec, ref.catalog, ref.jd);
		EXPECT_NEAR(correction.dRA_arcsec * std::cos(ref.dec), ref.raCosDecArcsec, 5.1e-4);
		EXPECT_NEAR(correction.dDEC_arcsec, ref.decArcsec, 5.1e-4);
	}
}

TEST(AstrometryBias, Fcct14DoesNotDebiasTycho2OrAct)
{
	const auto p = findFcct14Dat();
	if (p.empty())
		GTEST_SKIP() << "FCCT14 table not found";

	AstrometryBias bias;
	ASSERT_NO_THROW(bias.load(p));
	for (const char catalog : {'g', 'l'}) {
		const auto correction = bias.correction(1.0, 0.2, catalog, 2459000.0);
		EXPECT_DOUBLE_EQ(correction.dRA_arcsec, 0.0);
		EXPECT_DOUBLE_EQ(correction.dDEC_arcsec, 0.0);
	}
}

TEST(AstrometryBias, PreservesDoublePrecisionAndSubmittedCoordinates)
{
	const auto path = findFcct14Dat();
	if (path.empty())
		GTEST_SKIP() << "FCCT14 table not found";
	ObsNode observation;
	observation.RA = 29.66388 * Constants::PI / 180.0;
	observation.DEC = 15.54522 * Constants::PI / 180.0;
	observation.JD = 2400000.5 + 59175.17818;
	observation.catalogCode = 'o';
	const double ra = observation.RA;
	const double dec = observation.DEC;
	const double years = (observation.JD - kBiasReferenceJd) / 365.25;
	std::vector<ObsNode> observations{observation};
	ASSERT_EQ(Calculations::debiasObservations(observations, path), 1u);
	const auto& corrected = observations.front();
	ASSERT_TRUE(corrected.submittedCoordinates.has_value());
	EXPECT_EQ(corrected.submittedRaRad(), ra);
	EXPECT_EQ(corrected.submittedDecRad(), dec);
	EXPECT_NEAR(corrected.biasRaArcsec, (-0.188 - 5.73 * years / 1000.0) / std::cos(dec), 1e-15);
	EXPECT_NEAR(corrected.biasDecArcsec, 0.098 + 7.93 * years / 1000.0, 1e-15);
	const double correctedRa = corrected.RA;
	const double correctedDec = corrected.DEC;
	ASSERT_EQ(Calculations::debiasObservations(observations, path), 1u);
	EXPECT_EQ(corrected.RA, correctedRa);
	EXPECT_EQ(corrected.DEC, correctedDec);
}

TEST(AstrometryBias, ThrowsOnRowCountMismatch)
{
	namespace fs = std::filesystem;
	const fs::path tmp = fs::temp_directory_path() / "bias_tiny_rowcount.dat";
	{
		std::ofstream out(tmp);
		out << "! synthetic — only 5 rows, far short of NPIX\n";
		for (int r = 0; r < 5; ++r) {
			for (int k = 0; k < kRowColumns; ++k)
				out << "0.0 ";
			out << "\n";
		}
	}

	AstrometryBias bias;
	EXPECT_THROW(bias.load(tmp), std::runtime_error);
	fs::remove(tmp);
}

TEST(AstrometryBias, ThrowsOnTruncatedRow)
{
	namespace fs = std::filesystem;
	const fs::path tmp = fs::temp_directory_path() / "bias_truncated_row.dat";
	{
		std::ofstream out(tmp);
		for (int k = 0; k < kRowColumns / 2; ++k)
			out << "0.0 ";
		out << "\n";
	}

	AstrometryBias bias;
	EXPECT_THROW(bias.load(tmp), std::runtime_error);
	fs::remove(tmp);
}
