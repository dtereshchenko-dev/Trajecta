#include "obs/AstrometricResidual.h"

#include "obs/ObserverGeometry.h"
#include "orbital_core/Constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbital_core::obs {
	namespace {

		void errorEllipseToQuadratic(double a, double b, double angle, double& A, double& B, double& C) noexcept
		{
			const double cosAng = std::cos(angle);
			const double sinAng = std::sin(angle);
			const double a2 = a * a;
			const double b2 = b * b;

			A = -cosAng * cosAng / a2 - sinAng * sinAng / b2;
			B = -sinAng * cosAng * (1.0 / a2 - 1.0 / b2);
			C = -sinAng * sinAng / a2 - cosAng * cosAng / b2;
		}

		[[nodiscard]] bool quadraticToErrorEllipse(double A, double B, double C, double& a, double& b,
		                                           double& angle) noexcept
		{
			const double tval2 = (A - C) * (A - C) + 4.0 * B * B;
			const double tval = std::sqrt(tval2);
			const double eigen2 = (A + C - tval) * 0.5;
			if (eigen2 == 0.0)
				return false;

			const double eigen1 = (A * C - B * B) / eigen2;
			if (eigen1 >= 0.0 || eigen2 >= 0.0 || !std::isfinite(tval2))
				return false;

			a = 1.0 / std::sqrt(-eigen1);
			b = 1.0 / std::sqrt(-eigen2);
			angle = std::atan2(eigen1 - A, B);
			return std::isfinite(a) && std::isfinite(b) && std::isfinite(angle);
		}

		[[nodiscard]] bool adjustErrorEllipseForTiming(double& sigma1, double& sigma2, double& angle, double vx,
		                                               double vy) noexcept
		{
			if (sigma1 <= 0.0 || sigma2 <= 0.0)
				return false;

			double A = 0.0;
			double B = 0.0;
			double C = 0.0;
			errorEllipseToQuadratic(sigma1, sigma2, angle, A, B, C);

			const double E = A * vx * vx + 2.0 * B * vx * vy + C * vy * vy - 1.0;
			if (E == 0.0 || !std::isfinite(E))
				return false;

			const double Fx = A * vx + B * vy;
			const double Fy = C * vy + B * vx;
			const double A1 = A - Fx * Fx / E;
			const double B1 = B - Fx * Fy / E;
			const double C1 = C - Fy * Fy / E;
			return quadraticToErrorEllipse(A1, B1, C1, sigma1, sigma2, angle);
		}

	}  // namespace

	Eigen::Matrix<double, 2, 3> astrometricAngularJacobianEcl(const Vec3D& d, double observedCosDec) noexcept
	{
		const double xy2 = d.x * d.x + d.y * d.y;
		const double r2 = xy2 + d.z * d.z;
		if (!(xy2 > 0.0) || !std::isfinite(r2) || !std::isfinite(observedCosDec)) {
			return Eigen::Matrix<double, 2, 3>::Constant(std::numeric_limits<double>::quiet_NaN());
		}

		const double xy = std::sqrt(xy2);
		Eigen::Matrix<double, 2, 3> equatorial;
		equatorial << observedCosDec * d.y / xy2, -observedCosDec * d.x / xy2, 0.0, d.z * d.x / (xy * r2),
		    d.z * d.y / (xy * r2), -xy / r2;
		const double c = std::cos(Constants::EARTH_INCL);
		const double s = std::sin(Constants::EARTH_INCL);
		Eigen::Matrix3d eclToEq;
		eclToEq << 1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c;
		return equatorial * eclToEq;
	}

	AstrometricResidualTransform astrometricResidualTransform(const PreparedObservation& observation,
	                                                          const Vec3D& objectHelioEclAu,
	                                                          const Vec3D* objectVelocityEclAuPerDay) noexcept
	{
		double sigmaRa = observation.sigmaRa;
		double sigmaDec = observation.sigmaDec;
		double correlation = observation.sigmaCorr;
		if (observation.astrometricConvention == AstrometricConvention::GaiaFpr &&
		    observation.photocenterRadiusAu > 0.0) {
			const double distance =
			    (objectHelioEclAu - observation.earthR - observation.observerOffsetEclAu).magnitude();
			if (distance > 0.0 && std::isfinite(distance)) {
				const double extra = 0.3 * observation.photocenterRadiusAu / distance;
				const double covariance = correlation * sigmaRa * sigmaDec;
				sigmaRa = std::hypot(sigmaRa, extra);
				sigmaDec = std::hypot(sigmaDec, extra);
				correlation = covariance / (sigmaRa * sigmaDec);
			}
		}
		double raMotion = 0.0;
		double decMotion = 0.0;
		if (objectVelocityEclAuPerDay && observation.timeSigmaDays > 0.0) {
			const auto motion =
			    apparentMotionArcsecPerMinute(objectHelioEclAu, *objectVelocityEclAuPerDay, observation);
			if (motion.valid) {
				raMotion = motion.raArcsecPerMinute;
				decMotion = motion.decArcsecPerMinute;
			}
		}

		return astrometricResidualTransformFromMotion(sigmaRa, sigmaDec, correlation, observation.timeSigmaDays,
		                                              raMotion, decMotion);
	}

	AstrometricResidualTransform astrometricResidualTransformFromMotion(double sigmaRaRad, double sigmaDecRad,
	                                                                    double sigmaCorr, double timeSigmaDays,
	                                                                    double raArcsecPerMinute,
	                                                                    double decArcsecPerMinute) noexcept
	{
		double sigma1 = ((sigmaDecRad > 0.0) ? sigmaDecRad : Constants::ARCSEC_RAD) / Constants::ARCSEC_RAD;
		double sigma2 = ((sigmaRaRad > 0.0) ? sigmaRaRad : Constants::ARCSEC_RAD) / Constants::ARCSEC_RAD;
		double tilt = 0.0;

		if (sigmaCorr != 0.0 && std::isfinite(sigmaCorr)) {
			const double rho = std::clamp(sigmaCorr, -0.999999, 0.999999);
			const double denom = 1.0 - rho * rho;
			const double A = -1.0 / (denom * sigma1 * sigma1);
			const double B = rho / (denom * sigma1 * sigma2);
			const double C = -1.0 / (denom * sigma2 * sigma2);
			if (!quadraticToErrorEllipse(A, B, C, sigma1, sigma2, tilt))
				return {};
		}

		if (timeSigmaDays > 0.0 && (raArcsecPerMinute != 0.0 || decArcsecPerMinute != 0.0)) {
			const double timingMinutes = timeSigmaDays * 1440.0;
			(void)adjustErrorEllipseForTiming(sigma1, sigma2, tilt, timingMinutes * decArcsecPerMinute,
			                                  timingMinutes * raArcsecPerMinute);
		}

		if (sigma1 <= 0.0 || sigma2 <= 0.0 || !std::isfinite(sigma1) || !std::isfinite(sigma2) ||
		    !std::isfinite(tilt)) {
			return {};
		}

		const double invArcsec = 1.0 / Constants::ARCSEC_RAD;
		const double cosTilt = std::cos(tilt);
		const double sinTilt = std::sin(tilt);

		return AstrometricResidualTransform{sinTilt * invArcsec / sigma1, cosTilt * invArcsec / sigma1,
		                                    cosTilt * invArcsec / sigma2, -sinTilt * invArcsec / sigma2, true};
	}

	AstrometricResidualPair applyAstrometricResidualTransform(double raCosDecResidualRad, double decResidualRad,
	                                                          const AstrometricResidualTransform& transform) noexcept
	{
		if (!transform.valid)
			return {};

		const double first = transform.m00 * raCosDecResidualRad + transform.m01 * decResidualRad;
		const double second = transform.m10 * raCosDecResidualRad + transform.m11 * decResidualRad;
		if (!std::isfinite(first) || !std::isfinite(second))
			return {};

		return AstrometricResidualPair{first, second, std::hypot(first, second), true};
	}

	ResidualWeight residualWeightForReport(const PreparedObservation& observation, double raCosDecResidualRad,
	                                       double decResidualRad, bool useWeights,
	                                       const AstrometricResidualTransform& transform) noexcept
	{
		ResidualWeight result;
		result.weighted = useWeights;
		if (useWeights) {
			result.sigmaRaArcsec = observation.sigmaRa / Constants::ARCSEC_RAD;
			result.sigmaDecArcsec = observation.sigmaDec / Constants::ARCSEC_RAD;
			result.forceWeightRa = observation.forceWeightRa;
			result.forceWeightDec = observation.forceWeightDec;
			const auto pair = applyAstrometricResidualTransform(raCosDecResidualRad, decResidualRad, transform);
			if (pair.valid) {
				const double determinant = transform.m00 * transform.m11 - transform.m01 * transform.m10;
				const double scale = std::abs(determinant) * Constants::ARCSEC_RAD;
				if (scale > 0.0 && std::isfinite(scale)) {
					result.sigmaRaArcsec = std::hypot(transform.m11, transform.m01) / scale;
					result.sigmaDecArcsec = std::hypot(transform.m10, transform.m00) / scale;
					result.chiSquared = pair.first * pair.first + pair.second * pair.second;
					return result;
				}
			}
		}

		const double ra = raCosDecResidualRad / Constants::ARCSEC_RAD / result.sigmaRaArcsec;
		const double dec = decResidualRad / Constants::ARCSEC_RAD / result.sigmaDecArcsec;
		result.chiSquared = ra * ra + dec * dec;
		return result;
	}

}  // namespace orbital_core::obs
