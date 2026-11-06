#pragma once

#include "orbital_core/Vec3D.h"

#include <Eigen/Core>
#include <cmath>

namespace orbital_core::forces {
	inline void accumulatePointMass_dR(const Vec3D& q, double gm, Eigen::Matrix3d& M) noexcept
	{
		const double q2 = q.dot(q);
		if (q2 < 1e-24)
			return;

		const double qMag = std::sqrt(q2);
		const double q3 = q2 * qMag;
		const double k = gm / q3;
		const double s = 3.0 * k / q2;

		M(0, 0) += s * q.x * q.x - k;
		M(1, 1) += s * q.y * q.y - k;
		M(2, 2) += s * q.z * q.z - k;

		const double xy = s * q.x * q.y;
		const double xz = s * q.x * q.z;
		const double yz = s * q.y * q.z;

		M(0, 1) += xy;
		M(1, 0) += xy;
		M(0, 2) += xz;
		M(2, 0) += xz;
		M(1, 2) += yz;
		M(2, 1) += yz;
	}
}  // namespace orbital_core::forces
