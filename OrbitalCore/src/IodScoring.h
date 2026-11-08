#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/Vec3D.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Calculations {
	// Sum of squared unit-vector residuals for ranking initial orbit candidates.
	// Uses observer positions, but not light-time or astrometric corrections.
	// Approximate RMS in radians is sqrt(cost / N); returns 1e30 on failure.
	[[nodiscard]] ORBITAL_CORE_API double keplerResidualCost(const Vec3D& R1, const Vec3D& V1, double t1JD,
	                                                         const std::vector<ObsNode>& obs);
}  // namespace Calculations

namespace Calculations::iod_scoring {

	inline constexpr double kBadScore = 1.0e30;
	inline constexpr double kFailureScore = 10000.0;
	inline constexpr double kAcceptableExtensionScore = 10.0;

	struct ClassicalElementsForScore {
		double eccentricity{};
		double inclinationRad{};
		double semiMajorAxisAu{};
		double perihelionAu{};
	};

	[[nodiscard]] double scoreFromElements(double weightedRms, std::size_t observationCount,
	                                       const ClassicalElementsForScore& elements, bool heliocentric = true,
	                                       bool interstellar = false) noexcept;

	[[nodiscard]] double weightedRmsKepler(
	    const Vec3D& R, const Vec3D& V, double epochJD, const std::vector<ObsNode>& obs,
	    orbital_core::AstrometricWeightProfile weightProfile = orbital_core::AstrometricWeightProfile::Vfcc17) noexcept;

	[[nodiscard]] double scoreStateKepler(
	    const Vec3D& R, const Vec3D& V, double epochJD, const std::vector<ObsNode>& obs, bool interstellar = false,
	    orbital_core::AstrometricWeightProfile weightProfile = orbital_core::AstrometricWeightProfile::Vfcc17) noexcept;

	[[nodiscard]] double weightedRmsFromTotalSigma(const std::vector<double>& residualsSigma,
	                                               const std::vector<std::uint8_t>& included) noexcept;

	[[nodiscard]] double scoreOrbit(const Orbit_p& orbit, double weightedRms, std::size_t observationCount,
	                                bool interstellar = false) noexcept;

}  // namespace Calculations::iod_scoring
