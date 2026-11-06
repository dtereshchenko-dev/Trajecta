#pragma once

#include "orbital_core/OrbitalCoreExport.h"

class Orbit;

// PHA gates (IAU/MPC): q < 1.3 AU, Earth MOID < 0.05 AU, H < 22.
namespace orbital_core::hazard {
	struct PlanetMoidTable;
	inline constexpr double kPhaMoidThresholdAu{0.05};       // AU
	inline constexpr double kPhaPerihelionThresholdAu{1.3};  // AU
	inline constexpr double kPhaAbsoluteMagThreshold{22.0};  // mag

	struct ORBITAL_CORE_API PhaCriteria {
		bool isPha{false};
		bool perihelionOk{false};   // q < 1.3 AU (NEA orbit class)
		bool moidOk{false};         // |MOID(Earth)| < 0.05 AU
		bool absoluteMagOk{false};  // H < 22 mag
		double moidEarthAu{};       // computed Earth MOID, NaN on failure
		double perihelionAu{};      // input q (echo for UI)
		double semiMajorAxisAu{};   // input a (echo for UI / Amor reason)
		double absoluteMag{};       // input H (echo for UI)
	};

	// Uses the existing MOID table for this orbit and epoch; no ephemeris queries
	// or MOID calculations. A missing Earth MOID leaves the MOID gate false.
	[[nodiscard]] ORBITAL_CORE_API PhaCriteria classifyPha(const Orbit& orbit, double absoluteMag,
	                                                       const PlanetMoidTable& moids) noexcept;
}  // namespace orbital_core::hazard
