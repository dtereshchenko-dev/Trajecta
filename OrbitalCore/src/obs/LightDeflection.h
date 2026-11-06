#pragma once

#include "orbital_core/Vec3D.h"
#include "obs/PreparedObservation.h"
#include "obs/GaiaPhotocenter.h"
#include <Eigen/Core>

#include <algorithm>
#include <cmath>

namespace orbital_core::obs {

	// Schwarzschild radius of the Sun, 2 GM/c^2, in AU (SOFA SRS).
	inline constexpr double kSunSchwarzschildRadiusAu = 1.97412574336e-8;

	// Light deflection by a body, as SOFA iauLd. bodyMassRatio is the body mass in
	// solar masses; p is the unit direction from the observer to the source, q the
	// unit direction from the body to the source, e the unit direction from the
	// body to the observer, em the body-observer distance in AU and dlim the
	// deflection limiter. Returns the deflected direction, not renormalised.
	[[nodiscard]] inline Vec3D deflectLight(double bodyMassRatio, const Vec3D& p, const Vec3D& q, const Vec3D& e,
	                                        double em, double dlim) noexcept
	{
		const double w = bodyMassRatio * kSunSchwarzschildRadiusAu / em / std::max(q.dot(q + e), dlim);
		return p + p.cross(e.cross(q)) * w;
	}

	// Catalogue-referenced astrometry combines an undeflected star direction
	// with the measured object-star offset, including occultation astrometry.
	// The reported direction is therefore the light-time (astrometric) direction minus
	// the solar deflection of a star in that direction plus the smaller deflection
	// of the source itself (SOFA iauLd, limiter as iauLdsun).
	// topocentric is the observer-to-source vector at emission; source and observer
	// are heliocentric positions in the same frame, AU. The returned vector keeps
	// the topocentric length. Degenerate geometry returns the input unchanged.
	[[nodiscard]] inline Vec3D catalogReferencedDirection(const Vec3D& topocentric, const Vec3D& sourceHeliocentric,
	                                                      const Vec3D& observerHeliocentric) noexcept
	{
		const double rho = topocentric.magnitude();
		const double em = observerHeliocentric.magnitude();
		const double r = sourceHeliocentric.magnitude();
		if (!(rho > 0.0) || !(em > 0.0) || !(r > 0.0) || !std::isfinite(rho + em + r))
			return topocentric;
		const Vec3D p = topocentric / rho;
		const Vec3D e = observerHeliocentric / em;
		const Vec3D q = sourceHeliocentric / r;
		const double dlim = 1.0e-6 / std::max(em * em, 1.0);
		const Vec3D source = deflectLight(1.0, p, q, e, em, dlim) - p;
		const Vec3D star = deflectLight(1.0, p, p, e, em, dlim) - p;
		return (p + source - star) * rho;
	}

	[[nodiscard]] inline Vec3D absoluteDeflectedDirection(const Vec3D& topocentric, const Vec3D& sourceHeliocentric,
	                                                      const Vec3D& observerHeliocentric) noexcept
	{
		const double rho = topocentric.magnitude(), em = observerHeliocentric.magnitude(),
		             r = sourceHeliocentric.magnitude();
		if (!(rho > 0.0) || !(em > 0.0) || !(r > 0.0) || !std::isfinite(rho + em + r))
			return topocentric;
		return deflectLight(1.0, topocentric / rho, sourceHeliocentric / r, observerHeliocentric / em, em,
		                    1.0e-6 / std::max(em * em, 1.0)) *
		       rho;
	}

	// Optical catalogue positions are referenced to stars at effectively infinite
	// distance. Gaia spacecraft astrometry requires an explicitly identified
	// absolute convention and must not be inferred from the star catalogue.
	[[nodiscard]] inline Vec3D opticalObservationDirection(const Vec3D& sourceHeliocentric,
	                                                       const PreparedObservation& observation) noexcept
	{
		const Vec3D observer = observation.earthR + observation.observerOffsetEclAu;
		const Vec3D topocentric = sourceHeliocentric - observer;
		if (!observation.radar) {
			if (observation.astrometricConvention != AstrometricConvention::CatalogReferenced) {
				Vec3D direction = absoluteDeflectedDirection(topocentric, sourceHeliocentric, observer);
				if (observation.astrometricConvention == AstrometricConvention::GaiaFpr)
					direction +=
					    gaiaPhotocenterCorrection(sourceHeliocentric, observer, observation.photocenterRadiusAu).offset;
				return direction;
			}
			return catalogReferencedDirection(topocentric, sourceHeliocentric, observer);
		}
		return topocentric;
	}

	// Derivative of the same finite-distance, catalogue-referenced direction
	// with respect to the heliocentric source position, at fixed observer.
	[[nodiscard]] inline Eigen::Matrix3d opticalObservationDirectionJacobian(
	    const Vec3D& source, const PreparedObservation& observation) noexcept
	{
		if (observation.radar)
			return Eigen::Matrix3d::Identity();
		const Vec3D observer = observation.earthR + observation.observerOffsetEclAu;
		const Vec3D rhoVec = source - observer;
		const double rho = rhoVec.magnitude(), em = observer.magnitude(), r = source.magnitude();
		if (!(rho > 0.0) || !(em > 0.0) || !(r > 0.0) || !std::isfinite(rho + em + r))
			return Eigen::Matrix3d::Identity();
		const Eigen::Vector3d p(rhoVec.x / rho, rhoVec.y / rho, rhoVec.z / rho);
		const Eigen::Vector3d q(source.x / r, source.y / r, source.z / r);
		const Eigen::Vector3d e(observer.x / em, observer.y / em, observer.z / em);
		const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
		const Eigen::Matrix3d dp = (identity - p * p.transpose()) / rho;
		const Eigen::Matrix3d dq = (identity - q * q.transpose()) / r;
		const double limit = 1.0e-6 / std::max(em * em, 1.0);
		const double factor = kSunSchwarzschildRadiusAu / em;
		const double pe = p.dot(e), pq = p.dot(q);
		const double sourceDenomRaw = q.dot(q + e);
		const double sourceDenom = std::max(sourceDenomRaw, limit);
		const Eigen::Vector3d sourceNumerator = e * pq - q * pe;
		const Eigen::Matrix3d sourceNumeratorJac =
		    e * (q.transpose() * dp + p.transpose() * dq) - q * (e.transpose() * dp) - pe * dq;
		const Eigen::RowVector3d sourceDenomJac =
		    sourceDenomRaw > limit ? (2.0 * q + e).transpose() * dq : Eigen::RowVector3d::Zero().eval();
		const Eigen::Vector3d sourceShift = factor * sourceNumerator / sourceDenom;
		const Eigen::Matrix3d sourceShiftJac =
		    factor *
		    (sourceNumeratorJac / sourceDenom - sourceNumerator * sourceDenomJac / (sourceDenom * sourceDenom));
		if (observation.astrometricConvention != AstrometricConvention::CatalogReferenced) {
			Eigen::Matrix3d jacobian = (p + sourceShift) * p.transpose() + rho * (dp + sourceShiftJac);
			if (observation.astrometricConvention == AstrometricConvention::GaiaFpr)
				jacobian += gaiaPhotocenterCorrection(source, observer, observation.photocenterRadiusAu).derivative;
			return jacobian;
		}
		const double starDenomRaw = p.dot(p + e);
		const double starDenom = std::max(starDenomRaw, limit);
		const Eigen::Vector3d starNumerator = e * p.squaredNorm() - p * pe;
		const Eigen::Matrix3d starNumeratorJac = 2.0 * e * (p.transpose() * dp) - pe * dp - p * (e.transpose() * dp);
		const Eigen::RowVector3d starDenomJac =
		    starDenomRaw > limit ? (2.0 * p + e).transpose() * dp : Eigen::RowVector3d::Zero().eval();
		const Eigen::Vector3d starShift = factor * starNumerator / starDenom;
		const Eigen::Matrix3d starShiftJac =
		    factor * (starNumeratorJac / starDenom - starNumerator * starDenomJac / (starDenom * starDenom));
		return (p + sourceShift - starShift) * p.transpose() + rho * (dp + sourceShiftJac - starShiftJac);
	}

}  // namespace orbital_core::obs
