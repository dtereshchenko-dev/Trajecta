#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/Constants.h"
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <cmath>
#include <limits>
#include <numbers>

namespace orbital_core::cometary {
	namespace {
		constexpr double mu = Constants::GAUSS_K * Constants::GAUSS_K;
		constexpr double twoPi = 2.0 * std::numbers::pi;
		constexpr double angularTolerance = 64.0 * std::numeric_limits<double>::epsilon();

		double principal(double angle)
		{
			const double result = std::fmod(angle, twoPi);
			const double positive = result < 0.0 ? result + twoPi : result;
			return positive < twoPi ? positive : 0.0;
		}

		double anomalyDifference(double anomaly, bool hyperbolic)
		{
			if (std::abs(anomaly) > 0.5)
				return hyperbolic ? std::sinh(anomaly) - anomaly : anomaly - std::sin(anomaly);
			const double square = anomaly * anomaly;
			double term = anomaly * square / 6.0;
			double result = term;
			for (int order = 2; order <= 9; ++order) {
				term *= (hyperbolic ? square : -square) / ((2 * order) * (2 * order + 1));
				result += term;
			}
			return result;
		}
	}  // namespace

	std::optional<Vector6> fromCartesian(const Vector6& state)
	{
		if (!state.allFinite())
			return std::nullopt;
		const Eigen::Vector3d r = state.head<3>(), v = state.tail<3>();
		const double radius = r.stableNorm();
		const Eigen::Vector3d angular = r.cross(v);
		const double angularNorm = angular.stableNorm();
		if (!(radius > 0.0) || !std::isfinite(radius) || !(angularNorm > 0.0) || !std::isfinite(angularNorm))
			return std::nullopt;

		const Eigen::Vector3d radial = r / radius;
		const Eigen::Vector3d normal = angular / angularNorm;
		const Eigen::Vector3d eccentricity = v.cross(angular) / mu - radial;
		const double e = eccentricity.stableNorm();
		const double sinInclination = std::hypot(normal.x(), normal.y());
		if (!eccentricity.allFinite() || !std::isfinite(e) || !(e > angularTolerance))
			return std::nullopt;

		Eigen::Vector3d node = Eigen::Vector3d::UnitX();
		if (sinInclination > angularTolerance)
			node = Eigen::Vector3d(-normal.y(), normal.x(), 0.0) / sinInclination;
		const Eigen::Vector3d peri = eccentricity / e;
		const double angularOverK = angularNorm / Constants::GAUSS_K;
		const double q = (angularOverK / (1.0 + e)) * angularOverK;
		Vector6 elements;
		elements << q, e, std::atan2(sinInclination, normal.z()), principal(std::atan2(node.y(), node.x())),
		    principal(std::atan2(normal.dot(node.cross(peri)), node.dot(peri))),
		    principal(std::atan2(normal.dot(peri.cross(radial)), peri.dot(radial)));
		const double denominator = 1.0 + e * std::cos(elements(5));
		if (!elements.allFinite() || !(q > 0.0) || !(denominator > 0.0) || !std::isfinite(denominator) ||
		    !std::isfinite(q * (1.0 + e)))
			return std::nullopt;
		return elements;
	}

	std::optional<Vector6> toCartesian(const Vector6& elements, Matrix6* jacobian)
	{
		if (!elements.allFinite())
			return std::nullopt;
		const double q = elements(0), e = elements(1);
		if (!(q > 0.0) || !(e >= 0.0))
			return std::nullopt;
		const double ci = std::cos(elements(2)), si = std::sin(elements(2));
		const double cn = std::cos(elements(3)), sn = std::sin(elements(3));
		const double cw = std::cos(elements(4)), sw = std::sin(elements(4));
		const double cv = std::cos(elements(5)), sv = std::sin(elements(5));
		const double p = q * (1.0 + e);
		const double denominator = 1.0 + e * cv;
		if (!(p > 0.0) || !std::isfinite(p) || !(denominator > 0.0) || !std::isfinite(denominator))
			return std::nullopt;
		const double radius = p / denominator;
		const double speed = Constants::GAUSS_K / std::sqrt(p);
		if (!(radius > 0.0) || !std::isfinite(radius) || !(speed > 0.0) || !std::isfinite(speed))
			return std::nullopt;

		const Eigen::Vector3d peri(cn * cw - sn * sw * ci, sn * cw + cn * sw * ci, sw * si);
		const Eigen::Vector3d transverse(-cn * sw - sn * cw * ci, -sn * sw + cn * cw * ci, cw * si);
		Vector6 state;
		state.head<3>() = radius * (cv * peri + sv * transverse);
		state.tail<3>() = speed * (-sv * peri + (e + cv) * transverse);
		if (!state.allFinite())
			return std::nullopt;
		if (!jacobian)
			return state;
		if (!(e > angularTolerance) || !(std::abs(si) > angularTolerance))
			return std::nullopt;

		Matrix6 j;
		j.col(0).head<3>() = state.head<3>() / q;
		j.col(0).tail<3>() = -0.5 * (state.tail<3>() / q);
		j.col(1).head<3>() = state.head<3>() * ((1.0 - cv) / ((1.0 + e) * denominator));
		j.col(1).tail<3>() = speed * transverse - 0.5 * (state.tail<3>() / (1.0 + e));
		const Eigen::Vector3d dPeri(sn * sw * si, -cn * sw * si, sw * ci);
		const Eigen::Vector3d dTransverse(sn * cw * si, -cn * cw * si, cw * ci);
		j.col(2).head<3>() = radius * (cv * dPeri + sv * dTransverse);
		j.col(2).tail<3>() = speed * (-sv * dPeri + (e + cv) * dTransverse);
		j.col(3).head<3>() = Eigen::Vector3d::UnitZ().cross(state.head<3>());
		j.col(3).tail<3>() = Eigen::Vector3d::UnitZ().cross(state.tail<3>());
		j.col(4).head<3>() = radius * (cv * transverse - sv * peri);
		j.col(4).tail<3>() = speed * (-sv * transverse - (e + cv) * peri);
		j.col(5).head<3>() = state.head<3>() * (e * sv / denominator) + radius * (-sv * peri + cv * transverse);
		j.col(5).tail<3>() = -speed * (cv * peri + sv * transverse);
		if (!j.allFinite())
			return std::nullopt;
		*jacobian = j;
		return state;
	}

	std::optional<Matrix6> cartesianToElementsJacobian(const Vector6& elements)
	{
		Matrix6 forward;
		if (!toCartesian(elements, &forward))
			return std::nullopt;
		const Vector6 rowScale = forward.cwiseAbs().rowwise().maxCoeff();
		if (!(rowScale.array() > 0.0).all())
			return std::nullopt;
		Matrix6 scaled = forward;
		for (int row = 0; row < 6; ++row)
			scaled.row(row) /= rowScale(row);
		const Vector6 columnScale = scaled.cwiseAbs().colwise().maxCoeff().transpose();
		if (!(columnScale.array() > 0.0).all())
			return std::nullopt;
		for (int column = 0; column < 6; ++column)
			scaled.col(column) /= columnScale(column);
		const Eigen::FullPivLU<Matrix6> factorization(scaled);
		if (!factorization.isInvertible())
			return std::nullopt;
		Matrix6 result = factorization.solve(Matrix6::Identity());
		for (int row = 0; row < 6; ++row)
			result.row(row) /= columnScale(row);
		for (int column = 0; column < 6; ++column)
			result.col(column) /= rowScale(column);
		if (!result.allFinite())
			return std::nullopt;
		return result;
	}

	double timeSincePerihelion(double qAu, double e, double trueAnomalyRad)
	{
		const double invalid = std::numeric_limits<double>::quiet_NaN();
		if (!(qAu > 0.0) || !(e >= 0.0) || !std::isfinite(qAu) || !std::isfinite(e) || !std::isfinite(trueAnomalyRad))
			return invalid;
		const double nu = std::remainder(trueAnomalyRad, twoPi);
		const double denominator = 1.0 + e * std::cos(nu);
		if (!(denominator > 0.0) || !std::isfinite(denominator))
			return invalid;
		if (nu == 0.0)
			return nu;
		const double scale = qAu * std::sqrt(qAu) / Constants::GAUSS_K;
		double time;
		if (e == 1.0) {
			const double d = std::tan(0.5 * nu);
			time = std::sqrt(2.0) * scale * d * (1.0 + d * d / 3.0);
		} else if (e < 1.0) {
			const double delta = 1.0 - e;
			const double anomaly =
			    2.0 * std::atan2(std::sqrt(delta) * std::sin(0.5 * nu), std::sqrt(1.0 + e) * std::cos(0.5 * nu));
			const double mean = anomalyDifference(anomaly, false) + delta * std::sin(anomaly);
			time = scale * (mean / (delta * std::sqrt(delta)));
		} else {
			const double delta = e - 1.0;
			const double sinhAnomaly = std::sqrt(delta) * (std::sqrt(1.0 + e) * std::sin(nu) / denominator);
			const double anomaly = std::asinh(sinhAnomaly);
			const double mean = anomalyDifference(anomaly, true) + delta * sinhAnomaly;
			time = scale * (mean / (delta * std::sqrt(delta)));
		}
		return std::isfinite(time) ? time : invalid;
	}
}  // namespace orbital_core::cometary
