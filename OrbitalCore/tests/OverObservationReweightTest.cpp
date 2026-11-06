#include <gtest/gtest.h>

#include "obs/OverObservationReweight.h"

#include <cmath>
#include <vector>

namespace {
	orbital_core::obs::PreparedObservation observation(double jd, const char* station)
	{
		orbital_core::obs::PreparedObservation out;
		out.jdTdb = jd;
		out.observer.obsCode = station;
		out.note2 = 'C';
		return out;
	}
}  // namespace

TEST(OverObservationReweightTest, VeresBatchUsesEightHourGapsAndCeilingFour)
{
	std::vector<orbital_core::obs::PreparedObservation> obs;
	for (int i = 0; i < 5; ++i)
		obs.push_back(observation(2460000.0 + i * 0.1, "I52"));

	const auto factors = orbital_core::obs::computeSameStationOverObservationFactors(
	    obs, orbital_core::AstrometricWeightProfile::Vfcc17);

	ASSERT_EQ(factors.size(), obs.size());
	const double expected = std::sqrt(4.0 / 5.0);
	for (const auto& factor : factors) {
		EXPECT_DOUBLE_EQ(factor.ra, expected);
		EXPECT_DOUBLE_EQ(factor.dec, expected);
	}
}

TEST(OverObservationReweightTest, ForcedCoordinatesCountInBatchButAreNotScaled)
{
	std::vector<orbital_core::obs::PreparedObservation> obs;
	for (int i = 0; i < 5; ++i)
		obs.push_back(observation(2460000.0 + i * 0.1, "I52"));
	obs[0].forceWeightRa = true;
	obs[0].forceWeightDec = true;
	obs[1].forceWeightRa = true;

	const auto factors = orbital_core::obs::computeSameStationOverObservationFactors(
	    obs, orbital_core::AstrometricWeightProfile::Vfcc17);

	ASSERT_EQ(factors.size(), obs.size());
	const double expected = std::sqrt(4.0 / 5.0);
	EXPECT_DOUBLE_EQ(factors[0].ra, 1.0);
	EXPECT_DOUBLE_EQ(factors[0].dec, 1.0);
	EXPECT_DOUBLE_EQ(factors[1].ra, 1.0);
	EXPECT_DOUBLE_EQ(factors[1].dec, expected);
	for (std::size_t i = 2; i < factors.size(); ++i) {
		EXPECT_DOUBLE_EQ(factors[i].ra, expected);
		EXPECT_DOUBLE_EQ(factors[i].dec, expected);
	}
}

TEST(OverObservationReweightTest, ObserverSigmasCountInBatchButAreNotScaled)
{
	std::vector<orbital_core::obs::PreparedObservation> obs;
	for (int i = 0; i < 6; ++i)
		obs.push_back(observation(2460000.0 + i * 0.1, "C40"));
	obs[0].observerSigma = true;
	obs[1].observerSigma = true;

	const auto factors = orbital_core::obs::computeSameStationOverObservationFactors(
	    obs, orbital_core::AstrometricWeightProfile::Vfcc17);

	ASSERT_EQ(factors.size(), obs.size());
	const double expected = std::sqrt(4.0 / 6.0);
	for (std::size_t i = 0; i < 2; ++i) {
		EXPECT_DOUBLE_EQ(factors[i].ra, 1.0);
		EXPECT_DOUBLE_EQ(factors[i].dec, 1.0);
	}
	for (std::size_t i = 2; i < factors.size(); ++i) {
		EXPECT_DOUBLE_EQ(factors[i].ra, expected);
		EXPECT_DOUBLE_EQ(factors[i].dec, expected);
	}
}

TEST(OverObservationReweightTest, BatchOfObserverSigmasOnlyIsInactive)
{
	std::vector<orbital_core::obs::PreparedObservation> obs;
	for (int i = 0; i < 5; ++i) {
		obs.push_back(observation(2460000.0 + i * 0.01, "V37"));
		obs.back().observerSigma = true;
	}

	EXPECT_TRUE(
	    orbital_core::obs::computeSameStationOverObservationFactors(obs, orbital_core::AstrometricWeightProfile::Vfcc17)
	        .empty());
}

TEST(OverObservationReweightTest, GapGreaterThanEightHoursStartsNewBatch)
{
	std::vector<orbital_core::obs::PreparedObservation> obs;
	for (int i = 0; i < 4; ++i)
		obs.push_back(observation(2460000.0 + i * 0.05, "I52"));
	for (int i = 0; i < 4; ++i)
		obs.push_back(observation(2460001.0 + i * 0.05, "I52"));

	EXPECT_TRUE(
	    orbital_core::obs::computeSameStationOverObservationFactors(obs, orbital_core::AstrometricWeightProfile::Vfcc17)
	        .empty());
}

TEST(OverObservationReweightTest, UnknownTechniqueAndLeadingSatelliteDoNotStartBatch)
{
	std::vector<orbital_core::obs::PreparedObservation> obs;
	for (int i = 0; i < 7; ++i)
		obs.push_back(observation(2460000.0 + i * 0.01, "I52"));
	obs[0].isSatellite = true;
	obs[1].note2 = '?';
	obs[6].isSatellite = true;

	const auto factors = orbital_core::obs::computeSameStationOverObservationFactors(
	    obs, orbital_core::AstrometricWeightProfile::Vfcc17);
	ASSERT_EQ(factors.size(), obs.size());
	EXPECT_DOUBLE_EQ(factors[0].ra, 1.0);
	EXPECT_DOUBLE_EQ(factors[1].ra, 1.0);
	for (std::size_t i = 2; i < obs.size(); ++i) {
		EXPECT_DOUBLE_EQ(factors[i].ra, std::sqrt(4.0 / 5.0));
		EXPECT_DOUBLE_EQ(factors[i].dec, std::sqrt(4.0 / 5.0));
	}
}

TEST(OverObservationReweightTest, OtherWeightProfilesDoNotApplyVfcc17Rule)
{
	std::vector<orbital_core::obs::PreparedObservation> obs;
	for (int i = 0; i < 6; ++i)
		obs.push_back(observation(2460000.0 + i * 0.01, "I52"));

	EXPECT_TRUE(orbital_core::obs::computeSameStationOverObservationFactors(
	                obs, orbital_core::AstrometricWeightProfile::Uniform)
	                .empty());
}
