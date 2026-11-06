#pragma once

#include "orbital_core/Orbit.h"
#include <cmath>
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/hazard/Planet.h"

#include <algorithm>
#include <array>
#include <limits>

class IEphemerisSource;

// MOID via 256x256 coarse grid → 2D Newton-Raphson (Wisniowski & Rickman
// 2013).  Accurate to ~1e-9 AU for e ≤ 0.99; hyperbolic orbits return NaN.
namespace orbital_core::hazard {
	struct ORBITAL_CORE_API MoidResult {
		double moidAu{std::numeric_limits<double>::quiet_NaN()};
		double trueAnomaly1Rad{};
		double trueAnomaly2Rad{};
		int iterations{};
		bool converged{false};

		[[nodiscard]] bool isValid() const noexcept { return converged && std::isfinite(moidAu) && moidAu >= 0.0; }
	};

	[[nodiscard]] ORBITAL_CORE_API MoidResult computeMoid(const Orbit& a, const Orbit& b) noexcept;

	struct ORBITAL_CORE_API PlanetMoidTable {
		std::array<double, kPlanetCount> values = []() {
			std::array<double, kPlanetCount> a{};
			a.fill(std::numeric_limits<double>::quiet_NaN());
			return a;
		}();

		[[nodiscard]] double operator[](Planet p) const noexcept { return values[static_cast<std::size_t>(p)]; }

		[[nodiscard]] double& operator[](Planet p) noexcept { return values[static_cast<std::size_t>(p)]; }

		[[nodiscard]] bool any() const noexcept
		{
			return std::any_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); });
		}
	};

	// Planet orbits built fresh from ephemeris (R,V) at jdEpochTdb — no mean
	// elements.  `eph = nullptr` uses the global default; missing ephemeris
	// leaves all entries NaN.
	[[nodiscard]] ORBITAL_CORE_API PlanetMoidTable computeMoidsToAllPlanets(const Orbit& target, double jdEpochTdb,
	                                                                        IEphemerisSource* eph = nullptr) noexcept;
}  // namespace orbital_core::hazard
