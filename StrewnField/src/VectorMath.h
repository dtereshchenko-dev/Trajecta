#pragma once

#include "strewn_core/Types.h"
#include <cmath>
#include <numbers>

namespace strewn::detail {

	constexpr double pi = std::numbers::pi;
	constexpr double deg = pi / 180.0;
	constexpr double earthA = 6378137.0;
	constexpr double earthF = 1.0 / 298.257223563;
	constexpr double earthE2 = earthF * (2.0 - earthF);
	constexpr double earthMu = 3.986004418e14;
	constexpr double earthJ2 = 1.0826267e-3;
	constexpr double gravityRadius = 6378136.3;
	constexpr double omega = 7.292115146706979e-5;

	inline double norm(const Vector3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

	inline Vector3 cross(const Vector3& a, const Vector3& b)
	{
		return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
	}

	inline bool finite(const Vector3& v) { return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]); }

}  // namespace strewn::detail
