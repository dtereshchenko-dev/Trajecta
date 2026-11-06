#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <Eigen/Core>
#include <optional>

namespace orbital_core::cometary {
	using Vector6 = Eigen::Matrix<double, 6, 1>;
	using Matrix6 = Eigen::Matrix<double, 6, 6>;

	// COT = [q (AU), e, i, ascending node, argument of perihelion, true anomaly]
	// with all angles in radians. CAR = [position (AU), velocity (AU/day)] in
	// the same inertial frame; mu = Constants::GAUSS_K squared. These conversions
	// preserve the epoch and support elliptic, exactly parabolic and hyperbolic
	// conics without clamping eccentricity or introducing a semimajor axis.
	// Coordinate operations return nullopt for nonfinite or unrepresentable geometry,
	// q <= 0, e < 0, or 1 + e*cos(nu) <= 0.

	// Requires nonradial, noncircular motion. Returns i in [0, pi] and
	// node/peri/nu in [0, 2*pi). Uses node = 0 for equatorial orbits, including
	// retrograde ones. Circularity (e <= 64*double epsilon) is rejected;
	// abs(sin(i)) <= 64*double epsilon selects the equatorial convention.
	[[nodiscard]] ORBITAL_CORE_API std::optional<Vector6> fromCartesian(const Vector6& state);

	// Without a Jacobian, circular/equatorial elements are allowed with the
	// supplied angular convention. With a Jacobian, rejects those singular
	// charts at the same 64*epsilon thresholds and returns d(CAR)/d(COT);
	// the output matrix is unchanged on failure.
	[[nodiscard]] ORBITAL_CORE_API std::optional<Vector6> toCartesian(const Vector6& elements,
	                                                                  Matrix6* jacobian = nullptr);

	// Local d(COT)/d(CAR), not a derivative across a wrapped-angle discontinuity.
	// Inverts the analytic forward Jacobian using row/column-equilibrated full
	// pivoting; numerically rank-deficient charts return nullopt.
	[[nodiscard]] ORBITAL_CORE_API std::optional<Matrix6> cartesianToElementsJacobian(const Vector6& elements);

	// Conditioning policy shared by Orbit and coordinate charts: near-parabolic
	// semimajor axis / mean anomaly are ill-conditioned. This threshold selects
	// coordinates only; it neither classifies a conic nor clamps eccentricity.
	[[nodiscard]] constexpr bool usesCometaryElements(double e) noexcept { return e >= 0.99; }

	// Signed days from the nearest perihelion, reducing nu with remainder to
	// [-pi, pi]. Exactly parabolic motion uses Barker's equation; elliptic and
	// hyperbolic Kepler expressions use cancellation-resistant small-anomaly
	// series. q is AU; invalid/nonfinite geometry or an unrepresentable result
	// returns NaN. The e >= 0.99 selection policy does not affect this function.
	[[nodiscard]] ORBITAL_CORE_API double timeSincePerihelion(double qAu, double e, double trueAnomalyRad);
}  // namespace orbital_core::cometary
