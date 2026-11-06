#include "orbital_core/hazard/ImpactEnergy.h"
#include "orbital_core/Constants.h"
#include "orbital_core/hazard/PlanetData.h"

#include <cmath>
#include <limits>

namespace orbital_core::hazard {

	namespace {


		constexpr double kMegatonJoules = 4.184e15;

		constexpr double kBackgroundRate = 0.03;
		constexpr double kBackgroundExp = -0.8;

		constexpr double kDaysPerJulianYear = 365.25;

		[[nodiscard]] double sphereMassKg(double diameterKm, double densityGcc) noexcept
		{
			const double diameterM = diameterKm * 1000.0;
			const double densityKgM3 = densityGcc * 1000.0;
			const double volumeM3 = (Constants::PI / 6.0) * diameterM * diameterM * diameterM;
			return densityKgM3 * volumeM3;
		}

		[[nodiscard]] double energyFromMassMt(double massKg, double speedKmS) noexcept
		{
			const double speedMs = speedKmS * 1000.0;
			return 0.5 * massKg * speedMs * speedMs / kMegatonJoules;
		}

		[[nodiscard]] double backgroundPerYear(double energyMt) noexcept
		{
			if (!std::isfinite(energyMt) || energyMt <= 0.0)
				return 0.0;

			return kBackgroundRate * std::pow(energyMt, kBackgroundExp);
		}

		[[nodiscard]] int wholeYearsToImpact(double daysToImpact) noexcept
		{
			if (!std::isfinite(daysToImpact))
				return 0;

			const double years = std::ceil(daysToImpact / kDaysPerJulianYear);
			if (years <= 0.0)
				return 0;

			return years > 1.0 ? static_cast<int>(years) - 1 : 1;
		}

		struct TorinoDiagonals {
			double d1{0.0};
			double d2{0.0};
			double d3{0.0};
		};

		[[nodiscard]] TorinoDiagonals torinoDiagonals(double logP, double logE) noexcept
		{
			const double p = (logP + 2.0) / 2.0;
			return {(logE + 1.0) / 3.0 + p, (logE - 2.0) / 3.0 + p, (logE - 5.0) / 3.0 + p};
		}

		[[nodiscard]] bool insideTorinoDomain(double probabilityImpact, double energyMt) noexcept
		{
			return probabilityImpact >= 0.0 && probabilityImpact <= 1.0 && energyMt >= 0.0 && energyMt <= 1.0e8;
		}

	}  // namespace

	double impactSpeedKmS(double vInfKmS, Planet target) noexcept
	{
		const double mu = planetMuKm3PerS2(target);
		const double r = planetRadiusKm(target);
		const double escapeSq = (mu > 0.0 && r > 0.0) ? 2.0 * mu / r : 0.0;
		return std::sqrt(vInfKmS * vInfKmS + escapeSq);
	}

	double diameterFromMagnitudeKm(double absoluteMag, double albedo) noexcept
	{
		if (albedo <= 0.0)
			return 0.0;

		return 1329.0 / std::sqrt(albedo) * std::pow(10.0, -absoluteMag / 5.0);
	}

	double kineticEnergyMt(double diameterKm, double densityGcc, double speedKmS) noexcept
	{
		return energyFromMassMt(sphereMassKg(diameterKm, densityGcc), speedKmS);
	}

	EnergyClass classify(double energyMt) noexcept
	{
		if (energyMt < 1.0e-3)
			return EnergyClass::Trivial;
		if (energyMt < 0.1)
			return EnergyClass::Hiroshima;
		if (energyMt < 5.0)
			return EnergyClass::Chelyabinsk;
		if (energyMt < 50.0)
			return EnergyClass::Tunguska;
		if (energyMt < 1.0e3)
			return EnergyClass::Tsar;
		if (energyMt < 1.0e9)
			return EnergyClass::KPg;
		return EnergyClass::Sterilising;
	}

	int trueTorinoLevel(double probabilityImpact, double energyMt) noexcept
	{
		if (!insideTorinoDomain(probabilityImpact, energyMt))
			return -1;

		if (probabilityImpact == 0.0 || energyMt == 0.0)
			return 0;

		const double logP = std::log10(probabilityImpact);
		const double logE = std::log10(energyMt);
		const TorinoDiagonals d = torinoDiagonals(logP, logE);

		if (logP <= -2.0) {
			if (logE > 0.0 && d.d1 > 0.0 && d.d2 <= 0.0)
				return 1;
			if (d.d2 > 0.0 && d.d3 <= 0.0)
				return 2;
		}
		if (logP > -2.0 && probabilityImpact <= 0.99 && logE > 0.0 && logE <= 2.0)
			return 3;

		if (logP > -2.0 && logE >= 2.0 && d.d3 <= 0.0)
			return 4;

		if (probabilityImpact <= 0.99 && logE <= 5.0 && d.d3 > 0.0)
			return 5;

		if (logP <= -2.0 && d.d3 > 0.0)
			return 6;

		if (logP > -2.0 && probabilityImpact <= 0.99 && logE > 5.0)
			return 7;

		if (probabilityImpact > 0.99 && logE > 0.0 && logE <= 2.0)
			return 8;

		if (probabilityImpact > 0.99 && logE > 2.0 && logE <= 5.0)
			return 9;

		if (probabilityImpact > 0.99 && logE > 5.0)
			return 10;

		return 0;
	}

	PalermoRisk palermoRisk(const PalermoInput& in) noexcept
	{
		PalermoRisk out;

		const bool haveDensity = std::isfinite(in.densityGcc) && in.densityGcc > 0.0;
		const bool haveMass = in.massKg > 0.0;
		const bool haveDiameter = std::isfinite(in.diameterKm) && in.diameterKm > 0.0 && haveDensity;
		const bool haveMagnitude =
		    std::isfinite(in.absoluteMag) && std::isfinite(in.albedo) && in.albedo > 0.0 && haveDensity;

		const bool usable = std::isfinite(in.vInfKmS) && in.vInfKmS >= 0.0 && std::isfinite(in.probabilityImpact) &&
		                    in.probabilityImpact >= 0.0 && (haveMass || haveDiameter || haveMagnitude);

		if (!usable)
			return out;

		out.impactSpeedKmS = impactSpeedKmS(in.vInfKmS, in.target);

		if (haveMass) {
			out.massKg = in.massKg;
			out.diameterKm = haveDiameter ? in.diameterKm : 0.0;
		}

		else if (haveDiameter) {
			out.diameterKm = in.diameterKm;
			out.massKg = sphereMassKg(out.diameterKm, in.densityGcc);
		}

		else {
			out.diameterKm = diameterFromMagnitudeKm(in.absoluteMag, in.albedo);
			out.massKg = sphereMassKg(out.diameterKm, in.densityGcc);
		}

		out.energyMt = energyFromMassMt(out.massKg, out.impactSpeedKmS);
		out.expectedEnergyMt = out.energyMt * in.probabilityImpact;
		out.yearsUsed = wholeYearsToImpact(in.daysToImpact);

		if (in.target != Planet::Earth) {
			out.state = PalermoState::NotApplicable;
			return out;
		}

		out.backgroundPerYear = backgroundPerYear(out.energyMt);

		if (out.yearsUsed <= 0) {
			out.state = PalermoState::NoInterval;
			return out;
		}

		if (out.backgroundPerYear <= 0.0 || in.probabilityImpact <= 0.0) {
			out.scale = -std::numeric_limits<double>::infinity();
			return out;
		}

		out.relativeProbability = in.probabilityImpact / (out.backgroundPerYear * out.yearsUsed);
		out.scale = std::log10(out.relativeProbability);
		out.state = PalermoState::Computed;
		return out;
	}
}  // namespace orbital_core::hazard
