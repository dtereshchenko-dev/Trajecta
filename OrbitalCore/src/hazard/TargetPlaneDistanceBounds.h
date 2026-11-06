#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <optional>

namespace orbital_core::hazard::detail {

	struct DistanceBounds {
		double minimum;
		double maximum;
	};

	// Distances from the origin to a filled ellipse, in its principal-axis frame.
	// Semi-axes are ordered major/minor and already include the desired sigma scale.
	// Solve the two trust-region secular equations rather than sampling the boundary:
	// the nearest point of a thin, offset ellipse can fall between angular samples.
	inline std::optional<DistanceBounds> ellipseDistanceBounds(const Eigen::Vector2d& center,
	                                                           const Eigen::Vector2d& semiAxes)
	{
		if (!center.allFinite() || !semiAxes.allFinite() || semiAxes.y() < 0.0 || semiAxes.x() < semiAxes.y())
			return std::nullopt;
		const double scale = std::max(center.cwiseAbs().maxCoeff(), semiAxes.x());
		if (scale == 0.0)
			return DistanceBounds{0.0, 0.0};
		const Eigen::Vector2d q = center.cwiseAbs() / scale;
		const double a = semiAxes.x() / scale, b = semiAxes.y() / scale;
		const double a2 = a * a, b2 = b * b;
		if (a == b)
			return DistanceBounds{scale * std::max(0.0, q.norm() - a), scale * (q.norm() + a)};
		if (b2 == 0.0)
			return DistanceBounds{scale * std::hypot(std::max(0.0, q.x() - a), q.y()),
			                      scale * std::hypot(q.x() + a, q.y())};

		const double upper = std::hypot(a * q.x(), b * q.y());
		const auto solve = [upper](auto norm) {
			double low = 0.0, high = upper;
			for (int iteration = 0; iteration < 256; ++iteration) {
				const double mid = low + (high - low) * 0.5;
				if (mid == low || mid == high)
					break;
				if (norm(mid) > 1.0)
					low = mid;
				else
					high = mid;
			}
			return high;
		};

		double minimum = 0.0;
		if (std::hypot(q.x() / a, q.y() / b) > 1.0) {
			const double lambda =
			    solve([&](double value) { return std::hypot(a * q.x() / (a2 + value), b * q.y() / (b2 + value)); });
			minimum = std::hypot(q.x() * lambda / (a2 + lambda), q.y() * lambda / (b2 + lambda));
		}

		const double gap = (a - b) * (a + b);
		double maximum;
		if (q.x() == 0.0 && b * q.y() <= gap) {
			// Hard case: the maximum lies off the minor axis even though the center
			// lies on it. The Lagrange multiplier equals the largest squared axis.
			const double u = b * q.y() / gap;
			maximum = std::hypot(a * std::sqrt(std::max(0.0, 1.0 - u * u)), q.y() + b * u);
		} else {
			const double delta =
			    solve([&](double value) { return std::hypot(a * q.x() / value, b * q.y() / (gap + value)); });
			maximum = std::hypot(q.x() + a2 * q.x() / delta, q.y() + b2 * q.y() / (gap + delta));
		}
		DistanceBounds result{scale * minimum, scale * maximum};
		if (!std::isfinite(result.minimum) || !std::isfinite(result.maximum))
			return std::nullopt;
		return result;
	}

	// Two-body inverse focusing at the nominal v-infinity, alpha = mu / vInf^2.
	// This is a local, linearized B-plane estimate, not nonlinear clone propagation.
	inline double pericentreFromImpactParameter(double b, double alpha)
	{
		return b * (b / (std::hypot(b, alpha) + alpha));
	}

}  // namespace orbital_core::hazard::detail
