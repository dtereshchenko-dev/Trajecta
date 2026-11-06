#include "orbital_core/hazard/EquinoctialCoordinates.h"
#include "orbital_core/Constants.h"
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <cmath>
#include <limits>
#include <numbers>

namespace orbital_core::hazard::equinoctial {
	namespace {
		using Basis = Eigen::Matrix<double, 3, 2>;
		constexpr double twoPi = 2.0 * std::numbers::pi;

		double principal(double angle)
		{
			const double result = std::fmod(angle, twoPi);
			return result < 0.0 ? result + twoPi : result;
		}

		Basis planeBasis(double p, double q)
		{
			const double denominator = 1.0 + p * p + q * q;
			Basis basis;
			basis.col(0) << 1.0 - p * p + q * q, 2.0 * p * q, -2.0 * p;
			basis.col(1) << 2.0 * p * q, 1.0 + p * p - q * q, 2.0 * q;
			return basis / denominator;
		}

		std::optional<double> eccentricLongitude(double lambda, double h, double k)
		{
			lambda = std::remainder(lambda, twoPi);
			const double e = std::hypot(h, k);
			double lower = lambda - e, upper = lambda + e;
			double f = lambda;
			for (int iteration = 0; iteration < 100; ++iteration) {
				const double residual = (f - lambda) - k * std::sin(f) + h * std::cos(f);
				const double derivative = 1.0 - k * std::cos(f) - h * std::sin(f);
				const double step = residual / derivative;
				if (std::abs(step) <= 4.0 * std::numeric_limits<double>::epsilon() * (1.0 + std::abs(f)))
					return f - step;
				if (residual > 0.0)
					upper = f;
				else
					lower = f;
				const double next = f - step;
				f = next > lower && next < upper ? next : 0.5 * (lower + upper);
			}
			return std::nullopt;
		}
	}  // namespace

	std::optional<Vector6> fromCartesian(const Vector6& state)
	{
		if (!state.allFinite())
			return std::nullopt;
		const Eigen::Vector3d r = state.head<3>(), v = state.tail<3>();
		const double radius = r.norm();
		const Eigen::Vector3d angular = r.cross(v);
		const double angularNorm = angular.norm();
		if (!(radius > 0.0) || !(angularNorm > 0.0))
			return std::nullopt;
		const Eigen::Vector3d normal = angular / angularNorm;
		const double denominator = 1.0 + normal.z();
		if (!(denominator > 0.0))
			return std::nullopt;
		constexpr double mu = Constants::GAUSS_K * Constants::GAUSS_K;
		const double energyTwice = v.squaredNorm() - 2.0 * mu / radius;
		if (!(energyTwice < 0.0))
			return std::nullopt;
		const Eigen::Vector3d eccentricity = v.cross(angular) / mu - r / radius;
		const double p = normal.x() / denominator, q = -normal.y() / denominator;
		const Basis basis = planeBasis(p, q);
		const double h = eccentricity.dot(basis.col(1));
		const double k = eccentricity.dot(basis.col(0));
		const double e2 = h * h + k * k;
		if (!(e2 < 1.0))
			return std::nullopt;
		const double a = -mu / energyTwice;
		const double root = std::sqrt(1.0 - e2), beta = 1.0 / (1.0 + root);
		const Eigen::Vector2d xy = basis.transpose() * r;
		const double cross = h * k * beta;
		const double cosF = k + ((1.0 - k * k * beta) * xy.x() - cross * xy.y()) / (a * root);
		const double sinF = h + ((1.0 - h * h * beta) * xy.y() - cross * xy.x()) / (a * root);
		Vector6 elements;
		elements << a, h, k, p, q, principal(std::atan2(sinF, cosF) + h * cosF - k * sinF);
		if (!elements.allFinite())
			return std::nullopt;
		return elements;
	}

	std::optional<Vector6> toCartesian(const Vector6& elements, Matrix6* jacobian)
	{
		if (!elements.allFinite())
			return std::nullopt;
		const double a = elements(0), h = elements(1), k = elements(2);
		const double p = elements(3), q = elements(4);
		const double e2 = h * h + k * k;
		if (!(a > 0.0) || !(e2 < 1.0))
			return std::nullopt;
		const auto f = eccentricLongitude(elements(5), h, k);
		if (!f)
			return std::nullopt;
		const double root = std::sqrt(1.0 - e2), beta = 1.0 / (1.0 + root);
		const Eigen::Vector2d u(std::cos(*f), std::sin(*f));
		const Eigen::Vector2d du(-u.y(), u.x());
		Eigen::Matrix2d ellipse;
		ellipse << 1.0 - beta * h * h, beta * h * k, beta * h * k, 1.0 - beta * k * k;
		const double radiusOverA = 1.0 - k * u.x() - h * u.y();
		const double speed = Constants::GAUSS_K / (std::sqrt(a) * radiusOverA);
		const Eigen::Vector2d position = a * (ellipse * u - Eigen::Vector2d(k, h));
		const Eigen::Vector2d velocity = speed * (ellipse * du);
		const Basis basis = planeBasis(p, q);
		Vector6 state;
		state.head<3>() = basis * position;
		state.tail<3>() = basis * velocity;
		if (!state.allFinite())
			return std::nullopt;
		if (!jacobian)
			return state;

		Matrix6& j = *jacobian;
		j.col(0).head<3>() = state.head<3>() / a;
		j.col(0).tail<3>() = -state.tail<3>() / (2.0 * a);
		for (const int column : {1, 2, 5}) {
			const double dh = column == 1 ? 1.0 : 0.0;
			const double dk = column == 2 ? 1.0 : 0.0;
			const double df = (column == 5 ? 1.0 : dk * u.y() - dh * u.x()) / radiusOverA;
			const double dbeta = beta * beta * (h * dh + k * dk) / root;
			Eigen::Matrix2d dEllipse;
			dEllipse << -dbeta * h * h - 2.0 * beta * h * dh, dbeta * h * k + beta * (dh * k + h * dk),
			    dbeta * h * k + beta * (dh * k + h * dk), -dbeta * k * k - 2.0 * beta * k * dk;
			const double dRadius = -dk * u.x() - dh * u.y() + (k * u.y() - h * u.x()) * df;
			j.col(column).head<3>() = basis * (a * (dEllipse * u + ellipse * du * df - Eigen::Vector2d(dk, dh)));
			j.col(column).tail<3>() =
			    basis * (speed * (dEllipse * du - ellipse * u * df) - velocity * (dRadius / radiusOverA));
		}
		const double denominator = 1.0 + p * p + q * q;
		Basis dp, dq;
		dp.col(0) << -2.0 * p, 2.0 * q, -2.0;
		dp.col(1) << 2.0 * q, 2.0 * p, 0.0;
		dq.col(0) << 2.0 * q, 2.0 * p, 0.0;
		dq.col(1) << 2.0 * p, -2.0 * q, 2.0;
		dp = (dp - 2.0 * p * basis) / denominator;
		dq = (dq - 2.0 * q * basis) / denominator;
		j.col(3).head<3>() = dp * position;
		j.col(3).tail<3>() = dp * velocity;
		j.col(4).head<3>() = dq * position;
		j.col(4).tail<3>() = dq * velocity;
		if (!j.allFinite())
			return std::nullopt;
		return state;
	}

	std::optional<Matrix6> cartesianToElementsJacobian(const Vector6& elements)
	{
		Matrix6 reverse;
		if (!toCartesian(elements, &reverse))
			return std::nullopt;
		const Eigen::FullPivLU<Matrix6> factorization(reverse);
		if (!factorization.isInvertible())
			return std::nullopt;
		const Matrix6 result = factorization.inverse();
		if (!result.allFinite())
			return std::nullopt;
		return result;
	}
}  // namespace orbital_core::hazard::equinoctial
