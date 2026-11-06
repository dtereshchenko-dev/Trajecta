#pragma once

#include "orbital_core/Constants.h"
#include "orbital_core/Vec3D.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace orbital_core {

	// Reference close_app.f90::stepcon. The caller retains the step at encounter
	// entry as an upper bound, and restores adaptive control on encounter exit.
	inline double encounterStepDays(const Vec3D& r, const Vec3D& v, int npoint, double direction)
	{
		if (npoint <= 0)
			return std::numeric_limits<double>::infinity();
		const Vec3D h(r.y * v.z - r.z * v.y, r.z * v.x - r.x * v.z, r.x * v.y - r.y * v.x);
		if (!(h.magnitude() > 0.0))
			return std::copysign(std::numeric_limits<double>::infinity(), direction);
		const double angle = Constants::PI / npoint;
		double dt = r.dot(r) * angle / h.magnitude() * direction;
		const Vec3D next = r + v * dt;
		const double denominator = r.magnitude() * next.magnitude();
		if (denominator > 0.0) {
			const double actual = std::acos(std::clamp(r.dot(next) / denominator, -1.0, 1.0));
			if (actual > angle)
				dt *= angle / actual;
		}
		return std::copysign(std::max(std::abs(dt), 1.0e-5), direction);
	}

}  // namespace orbital_core
