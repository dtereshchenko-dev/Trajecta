#include "MovCoordinates.h"
#include "orbital_core/Constants.h"
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <unsupported/Eigen/Polynomials>
#include <unsupported/Eigen/AutoDiff>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core::hazard::mov {
	namespace {
		constexpr double lightSpeed = 299792.458 * 86400 / (Constants::AU_IN_METERS / 1000);
		constexpr double energyLimit = -Constants::GM_SUN / 200;

		Eigen::Vector3d eigen(const Vec3D& v) { return {v.x, v.y, v.z}; }

		Eigen::Matrix3d rotation()
		{
			const double c = std::cos(Constants::EARTH_INCL), s = std::sin(Constants::EARTH_INCL);
			Eigen::Matrix3d r;
			r << 1, 0, 0, 0, c, s, 0, -s, c;
			return r;
		}

		struct Directions {
			Eigen::Vector3d e, a, d, aa, ad;
		};

		Directions directions(double a, double d)
		{
			const double ca = std::cos(a), sa = std::sin(a), cd = std::cos(d), sd = std::sin(d);
			const auto r = rotation();
			return {r * Eigen::Vector3d(ca * cd, sa * cd, sd), r * Eigen::Vector3d(-sa * cd, ca * cd, 0),
			        r * Eigen::Vector3d(-ca * sd, -sa * sd, cd), r * Eigen::Vector3d(-ca * cd, -sa * cd, 0),
			        r * Eigen::Vector3d(sa * sd, -ca * sd, 0)};
		}

		std::vector<double> roots(const Eigen::VectorXd& coefficients)
		{
			Eigen::PolynomialSolver<double, Eigen::Dynamic> solver(coefficients);
			std::vector<double> values;
			for (const auto& root : solver.roots())
				if (root.real() > 0 && std::abs(root.imag()) < 1e-8 * std::max(1.0, root.real()))
					values.push_back(root.real());
			std::sort(values.begin(), values.end());
			values.erase(
			    std::unique(values.begin(), values.end(),
			                [](double a, double b) { return std::abs(a - b) < 1e-10 * std::max({1.0, a, b}); }),
			    values.end());
			return values;
		}

		Eigen::VectorXd product(const Eigen::VectorXd& a, const Eigen::VectorXd& b)
		{
			Eigen::VectorXd c = Eigen::VectorXd::Zero(a.size() + b.size() - 1);
			for (int i = 0; i < a.size(); ++i)
				for (int j = 0; j < b.size(); ++j)
					c(i + j) += a(i) * b(j);
			return c;
		}
	}  // namespace

	Attributable fitAttributable(std::span<const obs::PreparedObservation> obs)
	{
		if (obs.size() < 3)
			throw std::invalid_argument("MOV needs at least three optical observations");
		const double first = obs.front().mjdTdb;
		double sum = 0, time = 0;
		for (const auto& o : obs) {
			if (o.radar || o.isSatellite)
				throw std::invalid_argument("MOV currently requires ground-based optical observations");
			const double w = 1 / (o.sigmaRa * o.sigmaDec);
			if (!std::isfinite(w) || w <= 0)
				throw std::invalid_argument("Invalid MOV astrometric uncertainty");
			sum += w;
			time += (o.mjdTdb - first) * w;
		}
		Attributable result;
		result.epoch = first + time / sum;
		result.observer = obs.front().observer;
		for (const auto& o : obs)
			if (o.observer.obsCode != result.observer.obsCode) {
				result.observer = {};
				result.observer.obsCode = "500";
				break;
			}
		const double duration = obs.back().mjdTdb - first;
		if (!(duration > 0) || duration > 30)
			throw std::invalid_argument("MOV requires a short, time-ordered optical arc");
		Eigen::Matrix<double, 6, 1> polynomial;
		Eigen::Matrix<double, 6, 6> covariance = Matrix6::Zero();
		for (int axis = 0; axis < 2; ++axis) {
			Eigen::MatrixXd a(obs.size(), 3);
			Eigen::VectorXd y(obs.size());
			double previous = obs.front().raObs, unwrapped = previous;
			for (int i = 0; i < static_cast<int>(obs.size()); ++i) {
				const auto& o = obs[i];
				const double t = (o.mjdTdb - result.epoch) / duration;
				const double sigma = axis ? o.sigmaDec : o.sigmaRa / std::abs(o.raMetricCosDec);
				a.row(i) = Eigen::Vector3d(1, t, t * t / 2).transpose() / sigma;
				unwrapped += std::remainder(o.raObs - previous, 2 * Constants::PI);
				previous = o.raObs;
				y(i) = (axis ? o.decObs : unwrapped) / sigma;
			}
			const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);
			if (qr.rank() != 3)
				throw std::invalid_argument("MOV needs three distinct observing epochs");
			const Eigen::Vector3d scale(1, 1 / duration, 1 / (duration * duration));
			polynomial.segment<3>(3 * axis) = scale.asDiagonal() * qr.solve(y);
			const Eigen::Matrix3d inverse = (a.transpose() * a).ldlt().solve(Eigen::Matrix3d::Identity());
			covariance.block<3, 3>(3 * axis, 3 * axis) = scale.asDiagonal() * inverse * scale.asDiagonal();
		}
		result.angles << polynomial(0), polynomial(3), polynomial(1), polynomial(4);
		using AD = Eigen::AutoDiffScalar<Eigen::Matrix<double, 6, 1>>;
		std::array<AD, 6> p;
		for (int i = 0; i < 6; ++i) {
			p[i].value() = polynomial(i);
			p[i].derivatives() = Vector6::Unit(i);
		}
		const AD cd = cos(p[3]), sd = sin(p[3]);
		const AD eta2 = p[1] * p[1] * cd * cd + p[4] * p[4];
		if (!(eta2.value() > 1e-20))
			throw std::invalid_argument("MOV angular motion is unresolved");
		const AD curvature =
		    ((p[5] * p[1] - p[2] * p[4]) * cd + p[1] * sd * (eta2 + p[4] * p[4])) / (eta2 * sqrt(eta2));
		const double variance = curvature.derivatives().dot(covariance * curvature.derivatives());
		result.curvatureSnr = variance > 0 ? std::abs(curvature.value()) / std::sqrt(variance) : 0;
		return result;
	}

	Vector6 Chart::state(const Vector6& att) const
	{
		const auto q = observer(epoch + att(4) / lightSpeed);
		const auto u = directions(att(0), att(1));
		Vector6 x;
		x.head<3>() = eigen(q.R) + att(4) * u.e;
		x.tail<3>() = eigen(q.V) + att(5) * u.e + att(4) * (att(2) * u.a + att(3) * u.d);
		return x;
	}

	Matrix6 Chart::derivative(const Vector6& att) const
	{
		const auto u = directions(att(0), att(1));
		Matrix6 j = Matrix6::Zero();
		j.block<3, 1>(0, 0) = att(4) * u.a;
		j.block<3, 1>(0, 1) = att(4) * u.d;
		j.block<3, 1>(0, 4) = u.e;
		j.block<3, 1>(3, 0) = att(5) * u.a + att(4) * (att(2) * u.aa + att(3) * u.ad);
		j.block<3, 1>(3, 1) = att(5) * u.d + att(4) * (att(2) * u.ad - att(3) * u.e);
		j.block<3, 2>(3, 2) = j.block<3, 2>(0, 0);
		j.block<3, 1>(3, 4) = att(2) * u.a + att(3) * u.d;
		j.block<3, 1>(3, 5) = u.e;
		const double t = epoch + att(4) / lightSpeed;
		const auto q = observer(t);
		const double t0 = t - .5 / 86400, t1 = t + .5 / 86400;
		j.block<3, 1>(0, 4) += eigen(q.V) / lightSpeed;
		j.block<3, 1>(3, 4) += (eigen(observer(t1).V) - eigen(observer(t0).V)) / ((t1 - t0) * lightSpeed);
		return j;
	}

	Vector6 Chart::inverse(const Vector6& x) const
	{
		double rho = (x.head<3>() - eigen(observer(epoch).R)).norm();
		for (int i = 0; i < 12; ++i)
			rho = (x.head<3>() - eigen(observer(epoch + rho / lightSpeed).R)).norm();
		const auto q = observer(epoch + rho / lightSpeed);
		const Eigen::Vector3d r = rotation().transpose() * (x.head<3>() - eigen(q.R));
		const Eigen::Vector3d v = rotation().transpose() * (x.tail<3>() - eigen(q.V));
		const double xy = std::hypot(r.x(), r.y());
		if (!(xy > 1e-12) || !(rho > 0))
			throw std::runtime_error("Singular MOV angular chart");
		Vector6 att;
		att << std::atan2(r.y(), r.x()), std::atan2(r.z(), xy), (r.x() * v.y() - r.y() * v.x()) / (xy * xy),
		    (v.z() * xy - r.z() * (r.x() * v.x() + r.y() * v.y()) / xy) / (rho * rho), rho, r.dot(v) / rho;
		return att;
	}

	double AdmissibleRegion::discriminant(double range) const
	{
		const Eigen::Vector3d v = velocity + range * angularVelocity;
		return 2 * Constants::GM_SUN / (position + range * direction).norm() + 2 * energyLimit - v.squaredNorm() +
		       std::pow(direction.dot(velocity), 2);
	}

	bool AdmissibleRegion::contains(double range, double rate) const
	{
		if (!(range > 0))
			return false;
		if (std::none_of(ranges.begin(), ranges.end(), [&](auto r) { return range >= r.first && range <= r.second; }))
			return false;
		return std::pow(rate + direction.dot(velocity), 2) <= discriminant(range);
	}

	AdmissibleRegion admissibleRegion(const Chart& chart, const Eigen::Vector4d& att)
	{
		const auto q = chart.observer(chart.epoch);
		const auto u = directions(att(0), att(1));
		AdmissibleRegion region{eigen(q.R), eigen(q.V), u.e, att(2) * u.a + att(3) * u.d};
		const double a = region.angularVelocity.squaredNorm(), b = 2 * region.velocity.dot(region.angularVelocity);
		const double c =
		    region.velocity.squaredNorm() - std::pow(region.velocity.dot(region.direction), 2) - 2 * energyLimit;
		const double d = 2 * region.position.dot(region.direction), e = region.position.squaredNorm();
		const Eigen::Vector3d kinetic(c, b, a), radius(e, d, 1);
		Eigen::VectorXd polynomial = product(product(kinetic, kinetic), radius);
		polynomial(0) -= 4 * Constants::GM_SUN * Constants::GM_SUN;
		auto boundaries = roots(polynomial);
		boundaries.erase(
		    std::remove_if(boundaries.begin(), boundaries.end(),
		                   [&](double r) { return std::abs(region.discriminant(r)) > 1e-8 * (1 + a * r * r); }),
		    boundaries.end());
		double left = 1e-4;
		for (double right : boundaries) {
			if (right > left && region.discriminant((left + right) / 2) > 0)
				region.ranges.emplace_back(left, right);
			left = right;
		}
		if (region.ranges.empty())
			throw std::runtime_error("Empty admissible region");
		if (region.ranges.size() > 1) {
			region.ranges.front().first = std::max(region.ranges.front().first, 1e-3);
			if (region.ranges.front().first >= region.ranges.front().second)
				region.ranges.erase(region.ranges.begin());
		}
		Eigen::Vector2d slope(b, 2 * a);
		Eigen::VectorXd extrema = product(product(slope, slope), product(product(radius, radius), radius));
		extrema.head<3>() -= 4 * Constants::GM_SUN * Constants::GM_SUN * Eigen::Vector3d(d * d / 4, d, 1);
		auto critical = roots(extrema);
		for (auto r : region.ranges) {
			critical.push_back(r.first);
			critical.push_back(r.second);
		}
		double maximum = 0;
		for (double r : critical)
			if (std::any_of(region.ranges.begin(), region.ranges.end(),
			                [&](auto v) { return r >= v.first && r <= v.second; }))
				maximum = std::max(maximum, region.discriminant(r));
		region.rateMin = -region.direction.dot(region.velocity) - std::sqrt(maximum);
		region.rateMax = -region.direction.dot(region.velocity) + std::sqrt(maximum);
		return region;
	}

	double manifoldJacobian(const Eigen::MatrixXd& design)
	{
		if (design.cols() != 6 || !design.allFinite())
			throw std::runtime_error("Invalid MOV design matrix");
		Eigen::MatrixXd a = design.leftCols(4);
		const Eigen::Vector4d scale = a.colwise().norm();
		if ((scale.array() <= 0).any())
			throw std::runtime_error("Singular constrained normal matrix");
		a *= scale.cwiseInverse().asDiagonal();
		const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);
		if (qr.rank() != 4)
			throw std::runtime_error("Rank-deficient constrained fit");
		const Eigen::Matrix<double, 4, 2> tangent =
		    (-scale.cwiseInverse()).asDiagonal() * qr.solve(design.rightCols(2));
		const Eigen::JacobiSVD<Eigen::Matrix<double, 4, 2>> svd(tangent);
		const auto s = svd.singularValues();
		const double volume = std::hypot(1.0, s(0)) * std::hypot(1.0, s(1));
		if (!(volume > 0) || !std::isfinite(volume))
			throw std::runtime_error("Invalid MOV volume element");
		return volume;
	}

}  // namespace orbital_core::hazard::mov
