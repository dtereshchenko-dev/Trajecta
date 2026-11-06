#pragma once

#include "orbital_core/Vec3D.h"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>

namespace orbital_core::obs {
	struct PhotocenterCorrection {
		Vec3D offset;
		Eigen::Matrix3d derivative{Eigen::Matrix3d::Zero()};
	};

	// The phase law estimates the photocenter relative to the centre of mass.
	// Source and observer are heliocentric, in AU; radius is a physical radius.
	[[nodiscard]] inline PhotocenterCorrection gaiaPhotocenterCorrection(const Vec3D& source, const Vec3D& observer,
	                                                                     double radiusAu) noexcept
	{
		PhotocenterCorrection result;
		const Vec3D relative = source - observer;
		const double rho = relative.magnitude(), r = source.magnitude();
		if (!(radiusAu > 0.0) || !(rho > 0.0) || !(r > 0.0) || !std::isfinite(radiusAu + rho + r))
			return result;
		const Eigen::Vector3d p(relative.x / rho, relative.y / rho, relative.z / rho);
		const Eigen::Vector3d q(source.x / r, source.y / r, source.z / r);
		const double cosine = std::clamp(p.dot(q), -1.0, 1.0);
		const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
		if (!(sine > 1.0e-10))
			return result;
		const double phase = std::acos(cosine);
		constexpr double a = -0.02352667223191772, b = 0.05403589930067203;
		constexpr double c = 0.3318343597581827, d = -0.001233060221632019;
		const double amplitude = radiusAu * (((a * phase + b) * phase + c) * phase + d);
		const double slope = radiusAu * ((3.0 * a * phase + 2.0 * b) * phase + c);
		const Eigen::Vector3d direction = cosine * p - q;
		const Eigen::Vector3d offset = amplitude / sine * direction;
		result.offset = Vec3D(offset.x(), offset.y(), offset.z());
		const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
		const Eigen::Matrix3d dp = (identity - p * p.transpose()) / rho;
		const Eigen::Matrix3d dq = (identity - q * q.transpose()) / r;
		const Eigen::RowVector3d du = q.transpose() * dp + p.transpose() * dq;
		const Eigen::Matrix3d directionJac = p * du + cosine * dp - dq;
		const Eigen::RowVector3d factorJac = (amplitude * cosine / (sine * sine * sine) - slope / (sine * sine)) * du;
		result.derivative = (amplitude / sine) * directionJac + direction * factorJac;
		return result;
	}
}  // namespace orbital_core::obs
