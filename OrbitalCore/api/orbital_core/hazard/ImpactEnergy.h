#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/hazard/Planet.h"

#include <cstdint>
#include <string_view>

// Impact severity for a close-approach solution: the energy the body would
// deliver, and where that puts it on the two published hazard scales.
// Energy-only — does NOT model entry angle, ablation, or target geology.
//
// This header is the interface. The working parts of the formulas — sphere
// mass, the megaton constant, escape speed, the background impact frequency,
// the whole-year divisor of the Palermo scale — live in ImpactEnergy.cpp,
// together with the notes on how each was checked against the reference.
namespace orbital_core::hazard {

	enum class EnergyClass : std::uint8_t {
		// Below 1 kt — disintegrates harmlessly in the upper atmosphere
		// (most fireballs reported by NEO survey).
		Trivial = 0,

		// 1 kt – 100 kt.  Small airburst class.  Hiroshima reference
		// (15 kt fission yield).
		Hiroshima,

		// 100 kt – 5 Mt.  Chelyabinsk (2013) was ~500 kt.  Class covers
		// regional damage from airbursts.
		Chelyabinsk,

		// 5 Mt – 50 Mt.  Tunguska (1908) ~12 Mt — destroyed 2000 km² of
		// forest.  Region-scale damage on impact, transcontinental
		// climate impact possible.
		Tunguska,

		// 50 Mt – 1 Gt (1000 Mt).  Above any historical weapon test
		// (Tsar Bomba 50 Mt); civilisation-scale regional event.
		Tsar,

		// 1 Gt – 1 Pt (10⁹ Mt).  Continental-scale damage up through
		// the K-Pg / Chicxulub event (~10⁸ Mt = 100 Tt — solidly inside
		// this band).  Mass-extinction territory at the upper end.
		KPg,

		// > 1 Pt (10⁹ Mt).  Beyond K-Pg.  Sterilising impact; preserved
		// here as an explicit ceiling so the UI doesn't clamp silently.
		Sterilising,
		Count_
	};

	[[nodiscard]] constexpr std::string_view label(EnergyClass c) noexcept
	{
		switch (c) {
			case EnergyClass::Trivial:
				return "Trivial";
			case EnergyClass::Hiroshima:
				return "Hiroshima";
			case EnergyClass::Chelyabinsk:
				return "Chelyabinsk";
			case EnergyClass::Tunguska:
				return "Tunguska";
			case EnergyClass::Tsar:
				return "Tsar";
			case EnergyClass::KPg:
				return "K-Pg";
			case EnergyClass::Sterilising:
				return "Sterilising";
			case EnergyClass::Count_:
				break;
		}
		return "?";
	}

	[[nodiscard]] ORBITAL_CORE_API EnergyClass classify(double energyMt) noexcept;

	[[nodiscard]] ORBITAL_CORE_API double diameterFromMagnitudeKm(double absoluteMag, double albedo) noexcept;

	// Speed at the surface for a body that approached `target` with this
	// asymptotic speed: v_imp = sqrt(v_inf^2 + v_esc^2), v_esc = sqrt(2mu/R)
	// from the shared planet table.
	//
	// Same identity as v_inf * gravitationalFocusing(v_inf, target) from
	// BPlane, but that factor is clamped to 1 below v_inf = 1e-6 km/s, which
	// is right for a capture cross-section and wrong here: a body arriving
	// from rest still lands at v_esc.  Energy callers want this one.
	[[nodiscard]] ORBITAL_CORE_API double impactSpeedKmS(double vInfKmS, Planet target) noexcept;

	// Megatons of TNT for a sphere of this diameter and density arriving at
	// `speedKmS`.  Pass the IMPACT speed, not the encounter v_inf - see
	// impactSpeedKmS above.
	[[nodiscard]] ORBITAL_CORE_API double kineticEnergyMt(double diameterKm, double densityGcc,
	                                                      double speedKmS) noexcept;

	// Torino Impact Hazard Scale, 0..10.  Returns -1 outside the scale's
	// domain: probability not in [0,1], or energy negative or above 1e8 Mt.
	//
	// EARTH ONLY.  The energy boundaries are Earth damage regimes - the
	// "no consequence" floor exists because our atmosphere breaks small
	// bodies up, and levels 3, 7 and 10 are worded in terms of a threatened
	// population.  Neither survives a move to the Moon or Mars, so callers
	// must not display a level for another target body.  There is no planet
	// argument on purpose: the scale has no reading to give for one.
	[[nodiscard]] ORBITAL_CORE_API int trueTorinoLevel(double probabilityImpact, double energyMt) noexcept;

	enum class PalermoState : std::uint8_t {
		Computed = 0,
		// The encounter is at or before the epoch the interval is measured
		// from, so there is no span to divide by.  NOT the same as "soon":
		// an impact inside the first year divides by one year and does get
		// a reading.
		NoInterval,
		// Target is not Earth.  The physical fields are filled in and valid;
		// the scale is not, because its denominator is the Earth impact flux.
		NotApplicable,
		// No scale possible — zero probability or energy, or bad arguments.
		Undefined
	};

	// What a caller already has in hand for one risk row.  One argument
	// instead of seven positional doubles that are easy to transpose.
	struct PalermoInput {
		double vInfKmS{0.0};  // asymptotic relative speed, km/s
		double probabilityImpact{0.0};
		double daysToImpact{0.0};  // now → encounter
		double absoluteMag{0.0};   // H
		double albedo{0.0};        // impact-monitor setting
		double densityGcc{0.0};    // impact-monitor setting

		// Three ways to reach the mass, most direct first: a measured mass
		// wins, then a known diameter with the density, then H / albedo /
		// density.  Leave the shorter routes at zero to fall through.
		double massKg{0.0};
		double diameterKm{0.0};
		// Body that would be struck.  Sets the escape speed added to v_inf,
		// so mass, impact speed and energy are right for any target; the
		// scale itself is only defined for Earth.
		Planet target{Planet::Earth};
	};

	// The whole chain, so the UI can show how the value was reached and tests
	// can check it step by step.  `scale` means something only when
	// `state == PalermoState::Computed`.
	struct PalermoRisk {
		double diameterKm{0.0};
		double massKg{0.0};
		double impactSpeedKmS{0.0};
		double energyMt{0.0};
		double expectedEnergyMt{0.0};  // energy · probability
		double backgroundPerYear{0.0};
		double relativeProbability{0.0};
		double scale{0.0};
		int yearsUsed{0};
		PalermoState state{PalermoState::Undefined};
	};

	// Palermo Technical Impact Hazard Scale:
	//   PS = log10( P / (f_B(E) · ΔT) ),  f_B(E) = 0.03 · E^(-0.8) per year.
	// PS = 0 is "as likely as the background"; -2 is the usual threshold
	// below which a case is not worth publicising.
	//
	// f_B is the frequency of EARTH impacts, fitted to the flux Earth
	// collects with its own capture cross-section.  The Moon's is ~19x
	// smaller at v_inf = 17 km/s and ~66x smaller at 5 km/s, so the error
	// is not even a fixed offset; Mars draws on a different population
	// altogether.  For a target other than Earth the physical fields are
	// still computed and `state` comes back NotApplicable - rescaling f_B
	// would produce a number this scale does not define.
	[[nodiscard]] ORBITAL_CORE_API PalermoRisk palermoRisk(const PalermoInput& in) noexcept;
}  // namespace orbital_core::hazard
