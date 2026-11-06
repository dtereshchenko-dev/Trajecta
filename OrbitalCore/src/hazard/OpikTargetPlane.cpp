#include "hazard/OpikTargetPlane.h"

#include "orbital_core/hazard/PlanetData.h"

#include <Eigen/Geometry>

#include <array>
#include <cmath>

namespace orbital_core::hazard {
	namespace {
		using Matrix3 = Eigen::Matrix3d;
		using Matrix6 = Eigen::Matrix<double, 6, 6>;
		using Vector3 = Eigen::Vector3d;

		[[nodiscard]] double impactCrossSectionBodyRadii(double vInfKmS, Planet planet) noexcept
		{
			if (planet != Planet::Earth)
				return gravitationalFocusing(vInfKmS, planet);

			constexpr double impactRadiusKm = 6420.0;
			const double focusing =
			    std::sqrt(1.0 + 2.0 * planetMuKm3PerS2(planet) / (impactRadiusKm * vInfKmS * vInfKmS));
			return focusing * impactRadiusKm / planetRadiusKm(planet);
		}

		struct OpikEvaluation {
			Vector3 rTp{Vector3::Zero()};
			Vector3 vTp{Vector3::Zero()};
			Vector3 s{Vector3::Zero()};
			Vector3 xi{Vector3::Zero()};
			Vector3 zeta{Vector3::Zero()};
			Eigen::Vector2d coords{Eigen::Vector2d::Zero()};
			double vInf{0.0};
			double d{0.0};
			double speed{0.0};
			double bOverD{0.0};
			double sinHalf{0.0};
			double cosHalf{0.0};
			Vector3 angularMomentumUnit{Vector3::Zero()};
			Matrix3 rotation{Matrix3::Zero()};
			bool valid{false};
		};

		[[nodiscard]] bool finite(double value) noexcept { return std::isfinite(value); }

		[[nodiscard]] bool finite(const Vec3D& value) noexcept
		{
			return finite(value.x) && finite(value.y) && finite(value.z);
		}

		[[nodiscard]] Vector3 toEigen(const Vec3D& value) noexcept { return {value.x, value.y, value.z}; }

		[[nodiscard]] Vec3D toVec3D(const Vector3& value) noexcept { return {value.x(), value.y(), value.z()}; }

		[[nodiscard]] Matrix3 crossMatrix(const Vector3& axis) noexcept
		{
			Matrix3 result;
			result << 0.0, -axis.z(), axis.y(), axis.z(), 0.0, -axis.x(), -axis.y(), axis.x(), 0.0;
			return result;
		}

		[[nodiscard]] Matrix3 rotationMinusAngle(const Vector3& axis, double sinAngle, double cosAngle) noexcept
		{
			return cosAngle * Matrix3::Identity() - sinAngle * crossMatrix(axis) +
			       (1.0 - cosAngle) * (axis * axis.transpose());
		}

		[[nodiscard]] OpikEvaluation evaluateOpikTargetPlane(const Vector3& r, const Vector3& v, double mu) noexcept
		{
			OpikEvaluation out;
			const double v2 = v.squaredNorm();
			if (!(v2 > 0.0))
				return out;

			const Vector3 angularMomentum = r.cross(v);
			const double angularMomentumMagnitude = angularMomentum.norm();
			if (!(angularMomentumMagnitude > 0.0))
				return out;

			const Vector3 rMtp = r - v * (r.dot(v) / v2);
			out.d = rMtp.norm();
			out.speed = std::sqrt(v2);
			if (!(out.d > 0.0))
				return out;

			const double denominator = out.d * v2 - 2.0 * mu;
			if (!(denominator > 0.0))
				return out;
			const double bOverDSquared = out.d * v2 / denominator;
			if (!(bOverDSquared > 0.0))
				return out;

			out.bOverD = std::sqrt(bOverDSquared);
			out.sinHalf = mu / (out.d * v2 - mu);
			if (!(std::abs(out.sinHalf) < 1.0))
				return out;
			out.cosHalf = std::sqrt(std::max(0.0, 1.0 - out.sinHalf * out.sinHalf));
			if (!(out.cosHalf > 0.0))
				return out;

			out.angularMomentumUnit = angularMomentum / angularMomentumMagnitude;
			out.rotation = rotationMinusAngle(out.angularMomentumUnit, out.sinHalf, out.cosHalf);
			const Vector3 rotatedR = out.rotation * r;
			const Vector3 rotatedV = out.rotation * v;
			out.rTp = rotatedR * out.bOverD;
			out.vTp = rotatedV / out.bOverD;

			out.vInf = out.vTp.norm();
			const double rho = std::hypot(out.vTp.x(), out.vTp.y());
			if (!(out.vInf > 0.0) || !(rho > 1.0e-12))
				return out;

			out.s = out.vTp / out.vInf;
			out.xi = {out.vTp.y() / rho, -out.vTp.x() / rho, 0.0};
			out.zeta = {-out.vTp.x() * out.vTp.z() / (out.vInf * rho), -out.vTp.y() * out.vTp.z() / (out.vInf * rho),
			            rho / out.vInf};
			out.coords = {out.rTp.dot(out.xi), out.rTp.dot(out.zeta)};
			out.valid = out.coords.array().isFinite().all();
			return out;
		}

		[[nodiscard]] int leviCivita(int i, int j, int k) noexcept
		{
			if (i == j || j == k || i == k)
				return 0;
			return ((i - j) * (j - k) * (k - i) / 2) > 0 ? 1 : -1;
		}

		[[nodiscard]] Matrix6 tpcarMtpcarJacobian(const Vector3& r, const Vector3& v, double mu,
		                                          const OpikEvaluation& evaluation) noexcept
		{
			const double d = evaluation.d;
			const double speed = evaluation.speed;
			const double v2 = speed * speed;
			const double bsd = evaluation.bOverD;
			const double sinHalf = evaluation.sinHalf;
			const double cosHalf = evaluation.cosHalf;
			const Vector3& axis = evaluation.angularMomentumUnit;
			const Matrix3& rotation = evaluation.rotation;
			const Vector3 xp = rotation * r;
			const Vector3 yp = rotation * v;
			const Vector3 xhat = r / d;
			const Vector3 yhat = v / speed;
			const double normJ = r.cross(v).norm();

			const Matrix3 dRotationDSin =
			    (-sinHalf / cosHalf) * (Matrix3::Identity() - axis * axis.transpose()) - crossMatrix(axis);

			std::array<Matrix3, 3> dRotationDAxis;
			for (int h = 0; h < 3; ++h) {
				Vector3 unit = Vector3::Zero();
				unit(h) = 1.0;
				dRotationDAxis[h] = -sinHalf * crossMatrix(unit) +
				                    (1.0 - cosHalf) * (unit * axis.transpose() + axis * unit.transpose());
			}

			Matrix3 dAxisDR = Matrix3::Zero();
			Matrix3 dAxisDV = Matrix3::Zero();
			for (int h = 0; h < 3; ++h) {
				for (int j = 0; j < 3; ++j) {
					double crossR = 0.0;
					double crossV = 0.0;
					for (int k = 0; k < 3; ++k) {
						crossR += static_cast<double>(leviCivita(j, k, h)) * v(k);
						crossV += static_cast<double>(leviCivita(k, j, h)) * r(k);
					}
					dAxisDR(h, j) = (crossR - speed / d * r(j) * axis(h)) / normJ;
					dAxisDV(h, j) = (crossV - d / speed * v(j) * axis(h)) / normJ;
				}
			}

			const double dgdd = -v2 / mu * sinHalf * sinHalf;
			const double dgdv = -2.0 * speed * d / mu * sinHalf * sinHalf;
			Matrix6 dRotatedDState = Matrix6::Zero();
			dRotatedDState.block<3, 3>(0, 0) = rotation;
			dRotatedDState.block<3, 3>(3, 3) = rotation;
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < 3; ++k) {
						dRotatedDState(i, j) += dRotationDSin(i, k) * r(k) * dgdd * xhat(j);
						dRotatedDState(i + 3, j + 3) += dRotationDSin(i, k) * v(k) * dgdv * yhat(j);
						dRotatedDState(i, j + 3) += dRotationDSin(i, k) * r(k) * dgdv * yhat(j);
						dRotatedDState(i + 3, j) += dRotationDSin(i, k) * v(k) * dgdd * xhat(j);
						for (int h = 0; h < 3; ++h) {
							dRotatedDState(i, j) += dRotationDAxis[h](i, k) * dAxisDR(h, j) * r(k);
							dRotatedDState(i, j + 3) += dRotationDAxis[h](i, k) * dAxisDV(h, j) * r(k);
							dRotatedDState(i + 3, j) += dRotationDAxis[h](i, k) * dAxisDR(h, j) * v(k);
							dRotatedDState(i + 3, j + 3) += dRotationDAxis[h](i, k) * dAxisDV(h, j) * v(k);
						}
					}
				}
			}

			const Vector3 xphat = xp / d;
			const Vector3 yphat = yp / speed;
			const double denominator = v2 * d - 2.0 * mu;
			const double dbsddd = -mu * v2 / (bsd * denominator * denominator);
			const double dbsddv = -2.0 * mu * speed * d / (bsd * denominator * denominator);
			Matrix6 scaleJacobian = Matrix6::Zero();
			for (int i = 0; i < 3; ++i) {
				scaleJacobian(i, i) = bsd;
				scaleJacobian(i + 3, i + 3) = 1.0 / bsd;
				for (int j = 0; j < 3; ++j) {
					scaleJacobian(i, j) += dbsddd * xphat(j) * xp(i);
					scaleJacobian(i, j + 3) = dbsddv * yphat(j) * xp(i);
					scaleJacobian(i + 3, j) = -dbsddd * xphat(j) * yp(i) / (bsd * bsd);
					scaleJacobian(i + 3, j + 3) += -dbsddv * yphat(j) * yp(i) / (bsd * bsd);
				}
			}
			return scaleJacobian * dRotatedDState;
		}

		[[nodiscard]] Eigen::Matrix<double, 2, 6> opikJacobian(const OpikEvaluation& evaluation) noexcept
		{
			Eigen::Matrix<double, 2, 6> result = Eigen::Matrix<double, 2, 6>::Zero();
			result.block<1, 3>(0, 0) = evaluation.xi.transpose();
			result.block<1, 3>(1, 0) = evaluation.zeta.transpose();

			const double x = evaluation.vTp.x();
			const double y = evaluation.vTp.y();
			const double z = evaluation.vTp.z();
			const double u = evaluation.vInf;
			const double rho2 = x * x + y * y;
			const double rho = std::sqrt(rho2);
			const double rho3 = rho * rho2;
			const double d1 = u * u * u * rho3;
			const double d2 = u * u * u * rho;

			const Vector3 dXiDx(-x * y / rho3, -y * y / rho3, 0.0);
			const Vector3 dXiDy(x * x / rho3, x * y / rho3, 0.0);
			const Vector3 dZetaDx(z * ((x * x - u * u) * rho2 + u * u * x * x) / d1, x * y * z * (rho2 + u * u) / d1,
			                      x * z * z / d2);
			const Vector3 dZetaDy(x * y * z * (rho2 + u * u) / d1, z * ((y * y - u * u) * rho2 + u * u * y * y) / d1,
			                      y * z * z / d2);

			result(0, 3) = dXiDx.dot(evaluation.rTp);
			result(0, 4) = dXiDy.dot(evaluation.rTp);
			result(1, 3) = dZetaDx.dot(evaluation.rTp);
			result(1, 4) = dZetaDy.dot(evaluation.rTp);
			return result;
		}
	}  // namespace

	OpikTargetPlanePoint computeOpikTargetPlanePoint(const Vec3D& rGeoKm, const Vec3D& vGeoKmS, Planet planet) noexcept
	{
		OpikTargetPlanePoint out;
		const double mu = getPlanetData(planet).muKm3PerS2;
		if (!(mu > 0.0) || !finite(rGeoKm) || !finite(vGeoKmS))
			return out;

		const OpikEvaluation evaluation = evaluateOpikTargetPlane(toEigen(rGeoKm), toEigen(vGeoKmS), mu);
		if (!evaluation.valid)
			return out;

		const double gamma = impactCrossSectionBodyRadii(evaluation.vInf, planet);
		if (!(gamma > 0.0) || !finite(gamma))
			return out;

		out.coordsKm = evaluation.coords;
		out.basis.s = toVec3D(evaluation.s);
		out.basis.xi = toVec3D(evaluation.xi);
		out.basis.zeta = toVec3D(evaluation.zeta);
		out.basis.valid = finite(out.basis.s) && finite(out.basis.xi) && finite(out.basis.zeta);
		out.gamma = gamma;
		out.vInfKmS = evaluation.vInf;
		out.valid = out.basis.valid;
		return out;
	}

	OpikTargetPlaneDifferential computeOpikTargetPlaneDifferential(const Vec3D& rGeoKm, const Vec3D& vGeoKmS,
	                                                               Planet planet) noexcept
	{
		OpikTargetPlaneDifferential out;
		const double mu = getPlanetData(planet).muKm3PerS2;
		if (!(mu > 0.0) || !finite(rGeoKm) || !finite(vGeoKmS))
			return out;

		const Vector3 r = toEigen(rGeoKm);
		const Vector3 v = toEigen(vGeoKmS);
		const OpikEvaluation evaluation = evaluateOpikTargetPlane(r, v, mu);
		if (!evaluation.valid)
			return out;

		out.point = computeOpikTargetPlanePoint(rGeoKm, vGeoKmS, planet);
		if (!out.point.valid)
			return out;

		const Matrix6 del1 = tpcarMtpcarJacobian(r, v, mu, evaluation);
		const Eigen::Matrix<double, 2, 6> del2 = opikJacobian(evaluation);
		out.dCoordsDGeoState = del2 * del1;
		out.valid = out.dCoordsDGeoState.array().isFinite().all();
		return out;
	}

}  // namespace orbital_core::hazard
