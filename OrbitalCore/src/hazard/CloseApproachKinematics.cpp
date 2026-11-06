#include "hazard/CloseApproachKinematics.h"

#include "orbital_core/Constants.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/hazard/PlanetData.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbital_core::hazard {
	namespace {
		constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
		constexpr double kAuDayToKmS = kAuKm / Constants::SECONDS_PER_DAY;
		constexpr double kKm3S2ToAuDay2 =
		    (Constants::SECONDS_PER_DAY * Constants::SECONDS_PER_DAY) / (kAuKm * kAuKm * kAuKm);

		[[nodiscard]] double planetMuAuDay2(Planet planet) noexcept
		{
			return planetMuKm3PerS2(planet) * kKm3S2ToAuDay2;
		}

		struct StumpffValues {
			double s0{};
			double s1{};
			double s2{};
			double s3{};
			bool valid{false};
		};

		[[nodiscard]] StumpffValues stumpffSeries(double psi, double alpha) noexcept
		{
			constexpr int kMaxTerms = 70;
			constexpr int kMaxHalvings = 30;
			constexpr double kBetaLimit = 100.0;
			constexpr double kTolerance = 100.0 * std::numeric_limits<double>::epsilon();

			double reducedPsi = psi;
			double beta = alpha * reducedPsi * reducedPsi;
			int halvings = 0;
			while (std::abs(beta) >= kBetaLimit && halvings < kMaxHalvings) {
				reducedPsi *= 0.5;
				beta = alpha * reducedPsi * reducedPsi;
				++halvings;
			}

			double term2 = reducedPsi * reducedPsi * 0.5;
			double term3 = term2 * reducedPsi / 3.0;
			double s2 = term2;
			double s3 = term3;
			for (int j = 1; j <= kMaxTerms; ++j) {
				term2 *= beta / ((2.0 * j + 1.0) * (2.0 * j + 2.0));
				s2 += term2;
				if (std::abs(term2) < kTolerance)
					break;
			}
			for (int j = 1; j <= kMaxTerms; ++j) {
				term3 *= beta / ((2.0 * j + 2.0) * (2.0 * j + 3.0));
				s3 += term3;
				if (std::abs(term3) < kTolerance)
					break;
			}

			double s1 = reducedPsi + alpha * s3;
			double s0 = 1.0 + alpha * s2;
			for (int i = 0; i < halvings; ++i) {
				const double doubledS0 = 2.0 * s0 * s0 - 1.0;
				const double doubledS1 = 2.0 * s0 * s1;
				s0 = doubledS0;
				s1 = doubledS1;
			}
			if (halvings > 0) {
				if (alpha == 0.0)
					return {};
				s2 = (s0 - 1.0) / alpha;
				s3 = (s1 - psi) / alpha;
			}

			const bool valid = std::isfinite(s0) && std::isfinite(s1) && std::isfinite(s2) && std::isfinite(s3);
			return {s0, s1, s2, s3, valid};
		}

		struct PericenterAnomaly {
			double psi{};
			double deltaDays{};
			bool valid{false};
		};

		[[nodiscard]] PericenterAnomaly solvePericenterAnomaly(double r0, double sigma0, double pericenter, double mu,
		                                                       double alpha) noexcept
		{
			constexpr int kMaxIterations = 60;
			constexpr double kTolerance = 10.0 * std::numeric_limits<double>::epsilon();
			if (!(r0 > 0.0) || !(mu > 0.0) || !std::isfinite(alpha) || std::abs(alpha) < 1.0e-18) {
				return {};
			}

			double psi = 0.0;
			double deltaDays = 0.0;
			const double a = -mu / alpha;
			const double eccentricity = 1.0 - pericenter / a;
			if (alpha <= 0.0) {
				if (!(a > 0.0) || !(eccentricity > 0.0))
					return {};
				const double cosE = std::clamp((1.0 - r0 / a) / eccentricity, -1.0, 1.0);
				const double eccentricAnomaly = std::acos(cosE);
				psi = -std::copysign(1.0, sigma0) * eccentricAnomaly / std::sqrt(-alpha);
				double meanAnomaly = eccentricAnomaly - eccentricity * std::sin(eccentricAnomaly);
				meanAnomaly = std::remainder(meanAnomaly, 2.0 * Constants::PI);
				deltaDays = meanAnomaly / (std::sqrt(-alpha * -alpha * -alpha) / mu);
			} else {
				if (!(eccentricity > 1.0))
					return {};
				const double coshF = (1.0 - r0 / a) / eccentricity;
				const double hyperbolicAnomaly = coshF > 1.0 ? std::acosh(coshF) : 0.0;
				psi = -std::copysign(1.0, sigma0) * hyperbolicAnomaly / std::sqrt(alpha);
				deltaDays = (eccentricity * std::sinh(hyperbolicAnomaly) - hyperbolicAnomaly) /
				            (std::sqrt(alpha * alpha * alpha) / mu);
			}

			for (int iteration = 0; iteration < kMaxIterations; ++iteration) {
				const StumpffValues s = stumpffSeries(psi, alpha);
				if (!s.valid)
					return {};
				const double radius = r0 * s.s0 + sigma0 * s.s1 + mu * s.s2;
				const double first = (r0 * alpha + mu) * s.s1 + sigma0 * s.s0;
				const double second = (r0 * alpha + mu) * s.s0 + sigma0 * alpha * s.s1;
				if (std::abs(first) < kTolerance || radius < pericenter)
					break;
				if (std::abs(second) < kTolerance)
					return {};
				const double correction = -first / second;
				psi += correction;
				if (std::abs(correction) < kTolerance)
					break;
			}

			const StumpffValues s = stumpffSeries(psi, alpha);
			if (!s.valid)
				return {};
			deltaDays = r0 * s.s1 + sigma0 * s.s2 + mu * s.s3;
			return {psi, deltaDays, std::isfinite(psi) && std::isfinite(deltaDays)};
		}
	}  // namespace

	double circularEarthVInfinityKmS(const Orbit& orbit) noexcept
	{
		const double a = orbit.getSemiMajorAxis();
		const double e = orbit.getExcentricity();
		const double i = orbit.getInclination();
		if (!(a > 0.0) || e < 0.0 || e >= 1.0 || !std::isfinite(i))
			return std::numeric_limits<double>::quiet_NaN();

		const double oneMinusE2 = std::max(0.0, 1.0 - e * e);
		const double v2 = 3.0 - 1.0 / a - 2.0 * std::sqrt(a * oneMinusE2) * std::cos(i);
		const double auPerDay = (v2 > 0.0 ? std::sqrt(v2) : 0.0) * Constants::GAUSS_K;
		return auPerDay * kAuDayToKmS;
	}

	double closeApproachVInfinityKmS(const Vec3D& vRelAuPerDay, double rAu, Planet planet) noexcept
	{
		const double rForEnergyAu = std::max(rAu, 1.0e-12);
		const double v2 = vRelAuPerDay.dot(vRelAuPerDay);
		const double escape2 = 2.0 * planetMuAuDay2(planet) / rForEnergyAu;
		const double vInf2 = v2 - escape2;
		if (!(vInf2 > 0.0))
			return 0.0;
		return std::sqrt(vInf2) * kAuDayToKmS;
	}

	void refreshCloseApproachKinematics(CloseApproachEvent& event) noexcept
	{
		const Vec3D rRel = event.rTarget - event.rPlanet;
		const Vec3D vRel = event.vTarget - event.vPlanet;
		const double d2 = rRel.dot(rRel);
		event.distAu = d2 > 0.0 ? std::sqrt(d2) : 0.0;
		event.vInfKmS = closeApproachVInfinityKmS(vRel, event.distAu, event.planet);
	}

	ORBITAL_CORE_API InteriorPericenterMapping mapInteriorEncounterToPericenter(const Vec3D& relativePosition,
	                                                                            const Vec3D& relativeVelocity,
	                                                                            Planet planet) noexcept
	{
		InteriorPericenterMapping out;
		const double mu = planetMuAuDay2(planet);
		const double r0 = relativePosition.magnitude();
		const double velocity2 = relativeVelocity.dot(relativeVelocity);
		const Vec3D angularMomentum = relativePosition.cross(relativeVelocity);
		const double angularMomentum2 = angularMomentum.dot(angularMomentum);
		if (!(mu > 0.0) || !(r0 > 0.0) || !(angularMomentum2 > 0.0))
			return out;

		const Vec3D eccentricityVector = relativeVelocity.cross(angularMomentum) / mu - relativePosition / r0;
		const double eccentricity = eccentricityVector.magnitude();
		const double pericenter = angularMomentum2 / (mu * (1.0 + eccentricity));
		const double sigma0 = relativePosition.dot(relativeVelocity);
		const double alpha = velocity2 - 2.0 * mu / r0;
		const PericenterAnomaly anomaly = solvePericenterAnomaly(r0, sigma0, pericenter, mu, alpha);
		if (!anomaly.valid)
			return out;

		const StumpffValues s = stumpffSeries(anomaly.psi, alpha);
		if (!s.valid)
			return out;
		const double radius = r0 * s.s0 + sigma0 * s.s1 + mu * s.s2;
		if (!(radius > 0.0))
			return out;

		const double f = 1.0 - mu * s.s2 / r0;
		const double g = anomaly.deltaDays - mu * s.s3;
		const double fDot = -mu * s.s1 / (r0 * radius);
		const double gDot = 1.0 - mu * s.s2 / radius;
		out.relativePosition = relativePosition * f + relativeVelocity * g;
		out.relativeVelocity = relativePosition * fDot + relativeVelocity * gDot;

		const double s4 = (s.s2 - anomaly.psi * anomaly.psi * 0.5) / alpha;
		const double s5 = (s.s3 - anomaly.psi * anomaly.psi * anomaly.psi / 6.0) / alpha;
		const double u = s.s2 * anomaly.deltaDays + mu * (anomaly.psi * s4 - 3.0 * s5);
		const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
		const Eigen::Vector3d x0(relativePosition.x, relativePosition.y, relativePosition.z);
		const Eigen::Vector3d y0(relativeVelocity.x, relativeVelocity.y, relativeVelocity.z);
		const Eigen::Vector3d x(out.relativePosition.x, out.relativePosition.y, out.relativePosition.z);
		const Eigen::Vector3d y(out.relativeVelocity.x, out.relativeVelocity.y, out.relativeVelocity.z);
		const Eigen::Vector3d xdd = -mu * x / (radius * radius * radius);
		const Eigen::Vector3d xdd0 = -mu * x0 / (r0 * r0 * r0);
		Eigen::Matrix<double, 3, 2> xAndY;
		xAndY << x, y;
		Eigen::Matrix<double, 2, 3> x0AndY0;
		x0AndY0.row(0) = x0.transpose();
		x0AndY0.row(1) = y0.transpose();

		Eigen::Matrix2d coefficients;
		coefficients << -(fDot * s.s1 + (f - 1.0) / r0) / r0, -fDot * s.s2, (f - 1.0) * s.s1 / r0, (f - 1.0) * s.s2;
		const Eigen::Matrix3d dXdX0 = identity * f + y * xdd0.transpose() * u + xAndY * coefficients * x0AndY0;

		coefficients << -fDot * s.s2, -(gDot - 1.0) * s.s2, (f - 1.0) * s.s2, g * s.s2;
		const Eigen::Matrix3d dXdY0 = identity * g - y * y0.transpose() * u + xAndY * coefficients * x0AndY0;

		coefficients << -fDot * (s.s0 / (radius * r0) + 1.0 / (radius * radius) + 1.0 / (r0 * r0)),
		    -(fDot * s.s1 + (gDot - 1.0) / radius) / radius, (fDot * s.s1 + (f - 1.0) / r0) / r0, fDot * s.s2;
		const Eigen::Matrix3d dYdX0 = identity * fDot + xdd * xdd0.transpose() * u + xAndY * coefficients * x0AndY0;

		coefficients << -(fDot * s.s1 + (gDot - 1.0) / radius) / radius, -(gDot - 1.0) * s.s1 / radius, fDot * s.s2,
		    (gDot - 1.0) * s.s2;
		const Eigen::Matrix3d dYdY0 = identity * gDot - xdd * y0.transpose() * u + xAndY * coefficients * x0AndY0;

		out.stateTransition.block<3, 3>(0, 0) = dXdX0;
		out.stateTransition.block<3, 3>(0, 3) = dXdY0;
		out.stateTransition.block<3, 3>(3, 0) = dYdX0;
		out.stateTransition.block<3, 3>(3, 3) = dYdY0;
		out.deltaDays = anomaly.deltaDays;
		out.valid = out.relativePosition.magnitude() > 0.0 && out.stateTransition.array().isFinite().all();
		return out;
	}

	std::optional<Eigen::Matrix<double, 6, 6>> InteriorPericenterMapping::eventJacobian(Planet planet) const noexcept
	{
		if (!valid)
			return std::nullopt;
		const Eigen::Vector3d r(relativePosition.x, relativePosition.y, relativePosition.z);
		const Eigen::Vector3d v(relativeVelocity.x, relativeVelocity.y, relativeVelocity.z);
		const double radius = r.norm();
		const double mu = planetMuAuDay2(planet);
		if (!(radius > 0.0) || !(mu > 0.0))
			return std::nullopt;
		Eigen::Matrix<double, 6, 1> flow;
		flow << v, -mu * r / (radius * radius * radius);
		Eigen::Matrix<double, 1, 6> normal;
		normal << v.transpose(), r.transpose();
		const double denominator = (normal * flow)(0);
		const double scale = v.squaredNorm() + mu / radius;
		if (!(denominator > 64.0 * std::numeric_limits<double>::epsilon() * scale))
			return std::nullopt;
		const Eigen::Matrix<double, 6, 6> result =
		    (Eigen::Matrix<double, 6, 6>::Identity() - flow * normal / denominator) * stateTransition;
		if (!result.allFinite())
			return std::nullopt;
		return result;
	}

}  // namespace orbital_core::hazard
