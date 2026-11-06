#include <gtest/gtest.h>

#include "fit/FitObservations.h"
#include "obs/ObservatoryCatalog.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace {
	constexpr double kTableSigmaArcsec = 0.6;

	orbital_core::obs::PreparedObservation observation(double jd, const char* station)
	{
		orbital_core::obs::PreparedObservation out;
		out.jdTdb = jd;
		out.observer.obsCode = station;
		out.note2 = 'C';
		out.sigmaRa = kTableSigmaArcsec * Constants::ARCSEC_RAD;
		out.sigmaDec = kTableSigmaArcsec * Constants::ARCSEC_RAD;
		return out;
	}

	std::vector<orbital_core::obs::PreparedObservation> sixInOneNight()
	{
		std::vector<orbital_core::obs::PreparedObservation> rows;
		for (int i = 0; i < 6; ++i)
			rows.push_back(observation(2460256.0 + i * 0.003, "C40"));
		return rows;
	}

	OrbitOptimizerConfig vfcc17()
	{
		OrbitOptimizerConfig config;
		config.useWeights = true;
		config.weightProfile = orbital_core::AstrometricWeightProfile::Vfcc17;
		config.logDiagnostics = false;
		return config;
	}

	double relaxedSigma(std::size_t batch)
	{
		return kTableSigmaArcsec * Constants::ARCSEC_RAD / std::sqrt(4.0 / static_cast<double>(batch));
	}
}  // namespace

TEST(FitObservationsTest, SubsetKeepsRelaxationCountedOverWholeSet)
{
	auto rows = sixInOneNight();
	orbital_core::fit::relaxSameStationBatches(rows, vfcc17());

	const std::vector<std::uint8_t> included{1, 1, 0, 0, 1, 1};
	const auto subset = orbital_core::fit::selectObservations(rows, included);

	ASSERT_EQ(subset.size(), 4u);
	for (const auto& row : subset) {
		EXPECT_DOUBLE_EQ(row.sigmaRa, relaxedSigma(6));
		EXPECT_DOUBLE_EQ(row.sigmaDec, relaxedSigma(6));
	}
}

TEST(FitObservationsTest, RelaxingASubsetRecountsItsBatches)
{
	auto rows = sixInOneNight();
	const std::vector<std::uint8_t> included{1, 1, 0, 0, 1, 1};
	auto subset = orbital_core::fit::selectObservations(rows, included);
	orbital_core::fit::relaxSameStationBatches(subset, vfcc17());

	for (const auto& row : subset)
		EXPECT_DOUBLE_EQ(row.sigmaRa, kTableSigmaArcsec * Constants::ARCSEC_RAD);
}

TEST(FitObservationsTest, RelaxationNeedsVfcc17AndObservationWeights)
{
	const double table = kTableSigmaArcsec * Constants::ARCSEC_RAD;

	auto unweighted = sixInOneNight();
	auto config = vfcc17();
	config.useWeights = false;
	orbital_core::fit::relaxSameStationBatches(unweighted, config);
	for (const auto& row : unweighted)
		EXPECT_DOUBLE_EQ(row.sigmaRa, table);

	auto uniform = sixInOneNight();
	config = vfcc17();
	config.weightProfile = orbital_core::AstrometricWeightProfile::Uniform;
	orbital_core::fit::relaxSameStationBatches(uniform, config);
	for (const auto& row : uniform)
		EXPECT_DOUBLE_EQ(row.sigmaRa, table);
}

TEST(FitObservationsTest, SelectionPreservesOrderAndRejectsMaskMismatch)
{
	std::vector<orbital_core::obs::PreparedObservation> rows;
	for (int i = 0; i < 4; ++i)
		rows.push_back(observation(2460256.0 + i, "F51"));

	const std::vector<std::uint8_t> included{0, 1, 0, 1};
	const auto subset = orbital_core::fit::selectObservations(rows, included);
	ASSERT_EQ(subset.size(), 2u);
	EXPECT_EQ(subset[0].jdTdb, rows[1].jdTdb);
	EXPECT_EQ(subset[1].jdTdb, rows[3].jdTdb);

	const std::vector<std::uint8_t> shortMask{1, 1};
	EXPECT_THROW((void)orbital_core::fit::selectObservations(rows, shortMask), std::invalid_argument);
}

TEST(ObservatoryCatalogTest, LoadsOnceFromNativeCatalogPath)
{
	ASSERT_TRUE(std::filesystem::is_regular_file(getExecutableDir() / "OBSCODE.dat"));
	const auto& catalog = orbital_core::obs::observatoryCatalog();
	EXPECT_EQ(&catalog, &orbital_core::obs::observatoryCatalog());
	ASSERT_TRUE(catalog.contains("500"));
	EXPECT_EQ(catalog.at("500").code, "500");
}
