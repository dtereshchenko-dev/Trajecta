#pragma once

#include <cmath>

namespace Calculations::detail {

	inline constexpr double kStumpffTaylorPsiThreshold = 1.0e-6;

	[[nodiscard]] inline double stumpffC(double psi) noexcept
	{
		if (psi > kStumpffTaylorPsiThreshold)
			return (1.0 - std::cos(std::sqrt(psi))) / psi;
		if (psi < -kStumpffTaylorPsiThreshold)
			return (std::cosh(std::sqrt(-psi)) - 1.0) / (-psi);
		return 0.5 - psi / 24.0 + psi * psi / 720.0;
	}

	[[nodiscard]] inline double stumpffS(double psi) noexcept
	{
		if (psi > kStumpffTaylorPsiThreshold) {
			const double s = std::sqrt(psi);
			return (s - std::sin(s)) / (s * s * s);
		}
		if (psi < -kStumpffTaylorPsiThreshold) {
			const double s = std::sqrt(-psi);
			return (std::sinh(s) - s) / (s * s * s);
		}
		return 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0;
	}

}  // namespace Calculations::detail
