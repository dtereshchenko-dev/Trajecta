#include <gtest/gtest.h>

#include "orbital_core/hazard/ImpactEnergy.h"
#include "orbital_core/hazard/PlanetData.h"

#include <cmath>

using namespace orbital_core::hazard;

namespace {
	constexpr double kDensityStony = 2.6;
	constexpr double kDensityIron = 7.8;
}  // namespace

TEST(ImpactEnergyTest, Diameter_FormulaClosedForm)
{
	EXPECT_NEAR(diameterFromMagnitudeKm(10.0, 1.00), 13.29, 1.0e-3);
	EXPECT_NEAR(diameterFromMagnitudeKm(15.0, 0.25), 2.658, 1.0e-3);
	EXPECT_NEAR(diameterFromMagnitudeKm(20.0, 0.04), 0.6645, 1.0e-4);
}

TEST(ImpactEnergyTest, Apophis_DiameterStonyAlbedo)
{
	const double D = diameterFromMagnitudeKm(19.7, 0.14);
	EXPECT_NEAR(D, 0.41, 0.02);
}

TEST(ImpactEnergyTest, Diameter_ZeroAlbedoIsZero)
{
	EXPECT_EQ(diameterFromMagnitudeKm(20.0, 0.0), 0.0);
	EXPECT_EQ(diameterFromMagnitudeKm(20.0, -0.1), 0.0);
}

TEST(ImpactEnergyTest, TunguskaScale)
{
	const double E = kineticEnergyMt(0.060, kDensityStony, 27.0);
	EXPECT_GT(E, 5.0);
	EXPECT_LT(E, 50.0);
	EXPECT_EQ(classify(E), EnergyClass::Tunguska);
}

TEST(ImpactEnergyTest, ChelyabinskScale)
{
	const double E = kineticEnergyMt(0.020, kDensityStony, 19.0);
	EXPECT_GT(E, 0.1);
	EXPECT_LT(E, 5.0);
	EXPECT_EQ(classify(E), EnergyClass::Chelyabinsk);
}

TEST(ImpactEnergyTest, KPgScale)
{
	const double E = kineticEnergyMt(10.0, kDensityStony, 20.0);
	EXPECT_GT(E, 1.0e7);
	EXPECT_LT(E, 1.0e9);
	EXPECT_EQ(classify(E), EnergyClass::KPg);
}

TEST(ImpactEnergyTest, Classify_BoundaryValues)
{
	EXPECT_EQ(classify(1.0e-4), EnergyClass::Trivial);
	EXPECT_EQ(classify(0.05), EnergyClass::Hiroshima);
	EXPECT_EQ(classify(2.0), EnergyClass::Chelyabinsk);
	EXPECT_EQ(classify(15.0), EnergyClass::Tunguska);
	EXPECT_EQ(classify(500.0), EnergyClass::Tsar);
	EXPECT_EQ(classify(1.0e4), EnergyClass::KPg);
	EXPECT_EQ(classify(1.0e8), EnergyClass::KPg);
	EXPECT_EQ(classify(1.0e10), EnergyClass::Sterilising);
}

TEST(ImpactEnergyTest, DensitySensitivity)
{
	const double Estony = kineticEnergyMt(0.5, kDensityStony, 20.0);
	const double Eiron = kineticEnergyMt(0.5, kDensityIron, 20.0);
	EXPECT_NEAR(Eiron / Estony, 3.0, 0.01);
}

TEST(TorinoScale, CharacterisationGrid)
{
	constexpr double energies[] = {1e-4, 1e-3, 1e-2, 0.5, 1.0, 5.0, 50.0, 100.0, 1e3, 1e5, 1e6, 1e8};

	struct Row {
		double probability;
		int level[12];
	};

	constexpr Row rows[] = {
	    {1e-12, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, {1e-10, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
	    {1e-08, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}, {1e-06, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}},
	    {1e-04, {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2}}, {1e-03, {0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 6}},
	    {1e-02, {0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 6, 6}}, {0.05, {0, 0, 0, 0, 0, 3, 3, 3, 4, 5, 7, 7}},
	    {0.5, {0, 0, 0, 0, 0, 3, 3, 3, 5, 5, 7, 7}},   {0.9, {0, 0, 0, 0, 0, 3, 3, 3, 5, 5, 7, 7}},
	    {0.98, {0, 0, 0, 0, 0, 3, 3, 3, 5, 5, 7, 7}},  {0.995, {0, 0, 0, 0, 0, 8, 8, 4, 9, 9, 10, 10}},
	    {1.0, {0, 0, 0, 0, 0, 8, 8, 4, 9, 9, 10, 10}},
	};

	for (const Row& row : rows) {
		for (std::size_t j = 0; j < std::size(energies); ++j) {
			EXPECT_EQ(trueTorinoLevel(row.probability, energies[j]), row.level[j])
			    << "p=" << row.probability << " E=" << energies[j] << " Mt";
		}
	}
}

TEST(TorinoScale, RejectsInputsOutsideTheScaleDomain)
{
	EXPECT_EQ(trueTorinoLevel(-0.1, 1.0), -1);
	EXPECT_EQ(trueTorinoLevel(1.1, 1.0), -1);
	EXPECT_EQ(trueTorinoLevel(0.5, -1.0), -1);
	EXPECT_EQ(trueTorinoLevel(0.5, 1e9), -1);

	EXPECT_EQ(trueTorinoLevel(0.0, 1.0), 0);
	EXPECT_EQ(trueTorinoLevel(0.5, 0.0), 0);
}

namespace {
	double escapeSpeedKmS(Planet p) { return std::sqrt(2.0 * planetMuKm3PerS2(p) / planetRadiusKm(p)); }

	constexpr double daysForYears(int years) { return (years + 0.5) * 365.25; }

	const double kTraceVInfKmS = std::sqrt(11.47 * 11.47 - 11.18 * 11.18);

	constexpr double kReferenceAlbedo = 0.154;
	constexpr double kReferenceDensityGcc = 2.6;

	PalermoInput traceInput(double probability, double days, double massKg)
	{
		PalermoInput in;
		in.vInfKmS = kTraceVInfKmS;
		in.probabilityImpact = probability;
		in.daysToImpact = days;
		in.massKg = massKg;
		return in;
	}

	PalermoInput magnitudeInput(double absoluteMag, double albedo, double densityGcc, double probability, double days)
	{
		PalermoInput in;
		in.vInfKmS = kTraceVInfKmS;
		in.probabilityImpact = probability;
		in.daysToImpact = days;
		in.absoluteMag = absoluteMag;
		in.albedo = albedo;
		in.densityGcc = densityGcc;
		return in;
	}

	double massForEnergyMt(double energyMt, double days)
	{
		const double perKg = palermoRisk(traceInput(1.0e-06, days, 1.0)).energyMt;
		return energyMt / perKg;
	}
}  // namespace

TEST(PalermoScale, YearDivisorCeilsThenDropsOne)
{
	const auto yearsFor = [](double days) { return palermoRisk(traceInput(1.0e-06, days, 2.33e6)).yearsUsed; };

	EXPECT_EQ(yearsFor(6607.4931326469959), 18);
	EXPECT_EQ(yearsFor(1.0), 1);
	EXPECT_EQ(yearsFor(300.0), 1);
	EXPECT_EQ(yearsFor(365.25), 1);
	EXPECT_EQ(yearsFor(365.30), 1);
	EXPECT_EQ(yearsFor(730.50), 1);
	EXPECT_EQ(yearsFor(730.60), 2);
	EXPECT_EQ(yearsFor(daysForYears(24)), 24);

	EXPECT_EQ(yearsFor(0.0), 0);
	EXPECT_EQ(yearsFor(-50.0), 0);
}

TEST(PalermoScale, ImpactInsideTheFirstYearStillScoresAgainstOneYear)
{
	for (const double days : {1.0, 300.0, 365.25}) {
		const PalermoRisk r = palermoRisk(traceInput(1.0e-04, days, 2.33e6));
		EXPECT_EQ(r.state, PalermoState::Computed) << "days=" << days;
		EXPECT_EQ(r.yearsUsed, 1) << "days=" << days;
		EXPECT_NEAR(r.relativeProbability, 1.0e-04 / r.backgroundPerYear, 1e-18) << "days=" << days;
	}

	EXPECT_NEAR(palermoRisk(traceInput(1.0e-04, 365.25, 2.33e6)).scale,
	            palermoRisk(traceInput(1.0e-04, 366.25, 2.33e6)).scale, 1e-15);
}

TEST(PalermoScale, EncounterNotAheadOfTheEpochHasNoInterval)
{
	for (const double days : {0.0, -50.0, -5000.0}) {
		const PalermoRisk r = palermoRisk(traceInput(1.0e-04, days, 2.33e6));
		EXPECT_EQ(r.state, PalermoState::NoInterval) << "days=" << days;
		EXPECT_GT(r.energyMt, 0.0) << "days=" << days;
	}
}

TEST(PalermoScale, CertainImpactScoresFarAboveBackground)
{
	constexpr double days = 3652.5;
	const double mass = massForEnergyMt(1000.0, days);

	const PalermoRisk r = palermoRisk(traceInput(1.0, days, mass));
	ASSERT_EQ(r.state, PalermoState::Computed);
	EXPECT_EQ(r.yearsUsed, 9);
	EXPECT_NEAR(r.energyMt, 1000.0, 1.0e-9);
	EXPECT_NEAR(r.expectedEnergyMt, r.energyMt, 1.0e-9);
	EXPECT_NEAR(r.relativeProbability, 1.0 / (r.backgroundPerYear * 9.0), 1e-12);
	EXPECT_GT(r.scale, 2.0);

	EXPECT_EQ(trueTorinoLevel(1.0, r.energyMt), 9);
}

TEST(PalermoScale, ReproducesReferenceTraceChain)
{
	const double prob = 3.53e-07 / 3.64e-02;
	const PalermoRisk r = palermoRisk(traceInput(prob, 6607.4931326469959, 2.33e6));

	ASSERT_EQ(r.state, PalermoState::Computed);
	EXPECT_EQ(r.yearsUsed, 18);
	EXPECT_NEAR(r.impactSpeedKmS, 11.47, 0.005);

	EXPECT_NEAR(r.energyMt / 3.64e-02, 1.0, 0.015);
	EXPECT_NEAR(r.backgroundPerYear / 4.25e-01, 1.0, 0.015);
	EXPECT_NEAR(r.relativeProbability / 1.27e-06, 1.0, 0.015);
	EXPECT_NEAR(r.expectedEnergyMt / 3.53e-07, 1.0, 0.015);

	EXPECT_NEAR(r.scale, -5.90, 0.01);
}

TEST(PalermoScale, ReproducesSecondReferenceTraceCase)
{
	const double prob = 1.30e-06 / 3.64e-02;
	const PalermoRisk r = palermoRisk(traceInput(prob, 6607.4905330949259, 2.33e6));

	ASSERT_EQ(r.state, PalermoState::Computed);
	EXPECT_EQ(r.yearsUsed, 18);
	EXPECT_NEAR(r.relativeProbability / 4.66e-06, 1.0, 0.015);
	EXPECT_NEAR(r.scale, -5.33, 0.01);
}

TEST(PalermoScale, MassFromMagnitudeMatchesReferenceAtItsImpliedAlbedo)
{
	const PalermoRisk r =
	    palermoRisk(magnitudeInput(27.26, kReferenceAlbedo, kReferenceDensityGcc, 1.0e-06, daysForYears(18)));

	ASSERT_EQ(r.state, PalermoState::Computed);
	EXPECT_NEAR(r.massKg, 2.33e6, 0.01e6);
	EXPECT_GT(r.diameterKm, 0.0);

	const double m20 =
	    palermoRisk(magnitudeInput(20.0, kReferenceAlbedo, kReferenceDensityGcc, 1.0e-06, daysForYears(18))).massKg;
	const double m25 =
	    palermoRisk(magnitudeInput(25.0, kReferenceAlbedo, kReferenceDensityGcc, 1.0e-06, daysForYears(18))).massKg;
	EXPECT_NEAR(m20 / m25, 1000.0, 1e-6);
}

TEST(PalermoScale, MatchesReferenceRiskTableRows)
{
	struct Row {
		double prob;
		double expectedEnergyMt;
		int years;
		double ps;
	};

	constexpr Row rows[] = {
	    {3.23e-06, 1.14e-07, 18, -6.38}, {9.45e-09, 3.34e-10, 23, -9.02}, {1.15e-08, 4.07e-10, 24, -8.96},
	    {3.14e-08, 1.11e-09, 24, -8.52}, {1.89e-04, 6.67e-06, 24, -4.74}, {5.80e-04, 2.05e-05, 24, -4.26},
	};

	for (const Row& row : rows) {
		const double energyMt = row.expectedEnergyMt / row.prob;
		const double days = daysForYears(row.years);
		const PalermoRisk r = palermoRisk(traceInput(row.prob, days, massForEnergyMt(energyMt, days)));

		ASSERT_EQ(r.state, PalermoState::Computed);
		EXPECT_EQ(r.yearsUsed, row.years);
		EXPECT_NEAR(r.scale, row.ps, 0.015) << "p=" << row.prob << " E=" << energyMt << " years=" << row.years;
	}
}

TEST(PalermoScale, UndefinedInputs)
{
	const double days = daysForYears(10);

	const PalermoRisk zeroProb = palermoRisk(traceInput(0.0, days, 2.33e6));
	EXPECT_EQ(zeroProb.state, PalermoState::Undefined);
	EXPECT_TRUE(std::isinf(zeroProb.scale));
	EXPECT_LT(zeroProb.scale, 0.0);

	PalermoInput negativeSpeed = magnitudeInput(22.0, 0.15, 2.6, 1.0e-04, days);
	negativeSpeed.vInfKmS = -1.0;
	EXPECT_EQ(palermoRisk(negativeSpeed).state, PalermoState::Undefined);

	EXPECT_EQ(palermoRisk(magnitudeInput(22.0, 0.0, 2.6, 1.0e-04, days)).state, PalermoState::Undefined);
	EXPECT_EQ(palermoRisk(magnitudeInput(22.0, 0.15, 0.0, 1.0e-04, days)).state, PalermoState::Undefined);
	EXPECT_EQ(palermoRisk(PalermoInput{}).state, PalermoState::Undefined);
}

TEST(PalermoScale, UnitRelativeRiskIsZeroAndScaleIsLogarithmic)
{
	constexpr int years = 10;
	const double days = daysForYears(years);
	const double mass = massForEnergyMt(5.0, days);

	const double fb = palermoRisk(traceInput(1.0e-06, days, mass)).backgroundPerYear;
	ASSERT_GT(fb, 0.0);

	EXPECT_NEAR(palermoRisk(traceInput(fb * years, days, mass)).scale, 0.0, 1e-12);
	EXPECT_NEAR(palermoRisk(traceInput(10.0 * fb * years, days, mass)).scale, 1.0, 1e-12);
}

TEST(PalermoScale, ImpactSpeedFoldsInEscapeSpeed)
{
	PalermoInput in = traceInput(1.0e-06, daysForYears(10), 2.33e6);

	const double vEscEarth = escapeSpeedKmS(Planet::Earth);
	EXPECT_NEAR(vEscEarth, 11.18, 5.0e-4);

	in.vInfKmS = 0.0;
	EXPECT_NEAR(palermoRisk(in).impactSpeedKmS, vEscEarth, 1e-9);

	in.vInfKmS = 30.0;
	EXPECT_NEAR(palermoRisk(in).impactSpeedKmS, std::sqrt(30.0 * 30.0 + vEscEarth * vEscEarth), 1e-9);
}

TEST(PalermoScale, UsesTheTargetBodyEscapeSpeed)
{
	for (const Planet target : {Planet::Earth, Planet::Moon, Planet::Mars, Planet::Jupiter}) {
		PalermoInput in = traceInput(1.0e-06, daysForYears(10), 2.33e6);
		in.target = target;
		in.vInfKmS = 17.0;

		const double vEsc = escapeSpeedKmS(target);
		EXPECT_NEAR(palermoRisk(in).impactSpeedKmS, std::sqrt(17.0 * 17.0 + vEsc * vEsc), 1e-9) << longLabel(target);
	}
}

TEST(PalermoScale, ZeroApproachSpeedStillImpactsAtEscapeSpeed)
{
	PalermoInput in = traceInput(1.0e-06, daysForYears(10), 2.33e6);
	in.vInfKmS = 0.0;

	in.target = Planet::Moon;
	EXPECT_NEAR(palermoRisk(in).impactSpeedKmS, escapeSpeedKmS(Planet::Moon), 1e-9);
	EXPECT_NEAR(escapeSpeedKmS(Planet::Moon), 2.376, 1e-3);

	in.target = Planet::Mars;
	EXPECT_NEAR(palermoRisk(in).impactSpeedKmS, escapeSpeedKmS(Planet::Mars), 1e-9);
	EXPECT_NEAR(escapeSpeedKmS(Planet::Mars), 5.022, 1e-3);
}

TEST(PalermoScale, LunarImpactCarriesLessEnergyThanTerrestrial)
{
	PalermoInput in = traceInput(1.0e-06, daysForYears(10), 2.33e6);
	in.vInfKmS = 17.0;

	in.target = Planet::Earth;
	const PalermoRisk earth = palermoRisk(in);

	in.target = Planet::Moon;
	const PalermoRisk moon = palermoRisk(in);

	EXPECT_EQ(moon.massKg, earth.massKg);
	EXPECT_LT(moon.impactSpeedKmS, earth.impactSpeedKmS);
	EXPECT_LT(moon.energyMt, earth.energyMt);

	const double ratio = earth.energyMt / moon.energyMt;
	EXPECT_NEAR(ratio, 1.405, 0.01);
}

TEST(PalermoScale, ReportsNoScaleForTargetsOtherThanEarth)
{
	for (const Planet target : {Planet::Mercury, Planet::Venus, Planet::Mars, Planet::Jupiter, Planet::Saturn,
	                            Planet::Uranus, Planet::Neptune, Planet::Moon}) {
		PalermoInput in = traceInput(1.0e-06, daysForYears(10), 2.33e6);
		in.target = target;
		in.vInfKmS = 17.0;

		const PalermoRisk r = palermoRisk(in);

		EXPECT_EQ(r.state, PalermoState::NotApplicable) << longLabel(target);

		EXPECT_GT(r.massKg, 0.0) << longLabel(target);
		EXPECT_GT(r.impactSpeedKmS, 0.0) << longLabel(target);
		EXPECT_GT(r.energyMt, 0.0) << longLabel(target);
		EXPECT_GT(r.expectedEnergyMt, 0.0) << longLabel(target);
		EXPECT_EQ(r.yearsUsed, 10) << longLabel(target);

		EXPECT_EQ(r.backgroundPerYear, 0.0) << longLabel(target);
		EXPECT_EQ(r.relativeProbability, 0.0) << longLabel(target);
		EXPECT_EQ(r.scale, 0.0) << longLabel(target);
	}
}

TEST(PalermoScale, TakesMassFromTheMostDirectInputAvailable)
{
	constexpr double days = 3652.5;
	constexpr double albedo = 0.154;
	constexpr double density = 2.6;
	constexpr double absMag = 22.0;

	const double diameter = diameterFromMagnitudeKm(absMag, albedo);
	ASSERT_GT(diameter, 0.0);

	PalermoInput fromMag;
	fromMag.vInfKmS = 17.0;
	fromMag.probabilityImpact = 1.0e-06;
	fromMag.daysToImpact = days;
	fromMag.absoluteMag = absMag;
	fromMag.albedo = albedo;
	fromMag.densityGcc = density;

	PalermoInput fromDiameter = fromMag;
	fromDiameter.absoluteMag = 0.0;
	fromDiameter.albedo = 0.0;
	fromDiameter.diameterKm = diameter;

	const PalermoRisk viaMag = palermoRisk(fromMag);
	const PalermoRisk viaDia = palermoRisk(fromDiameter);

	ASSERT_EQ(viaMag.state, PalermoState::Computed);
	ASSERT_EQ(viaDia.state, PalermoState::Computed);
	EXPECT_NEAR(viaDia.diameterKm, viaMag.diameterKm, 1e-15);
	EXPECT_NEAR(viaDia.massKg, viaMag.massKg, 1e-6);
	EXPECT_NEAR(viaDia.scale, viaMag.scale, 1e-12);

	PalermoInput fromMass = fromDiameter;
	fromMass.massKg = viaMag.massKg * 8.0;
	const PalermoRisk viaMass = palermoRisk(fromMass);
	ASSERT_EQ(viaMass.state, PalermoState::Computed);
	EXPECT_NEAR(viaMass.massKg, viaMag.massKg * 8.0, 1e-6);
	EXPECT_NEAR(viaMass.energyMt, viaMag.energyMt * 8.0, 1e-9);

	PalermoInput noDensity = fromDiameter;
	noDensity.densityGcc = 0.0;
	EXPECT_EQ(palermoRisk(noDensity).state, PalermoState::Undefined);
}

TEST(PalermoScale, DefaultsToEarth)
{
	EXPECT_EQ(PalermoInput{}.target, Planet::Earth);

	PalermoInput in = traceInput(1.0e-06, daysForYears(10), 2.33e6);
	EXPECT_EQ(palermoRisk(in).state, PalermoState::Computed);
}
