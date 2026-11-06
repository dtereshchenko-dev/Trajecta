#pragma once

#include "obs/PreparedObservation.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/Residuals.h"

#include <array>
#include <Eigen/Core>

namespace orbital_core::obs {

	struct AstrometricResidualTransform {
		double m00{0.0};
		double m01{0.0};
		double m10{0.0};
		double m11{0.0};
		bool valid{false};
	};

	struct AstrometricResidualPair {
		double first{0.0};
		double second{0.0};
		double norm{0.0};
		bool valid{false};
	};

	// d(observed - computed RA*cosDec, DEC)/d(apparent ecliptic position), rad/AU.
	// Returns non-finite entries at a singular direction; callers must fall back.
	[[nodiscard]] ORBITAL_CORE_API Eigen::Matrix<double, 2, 3> astrometricAngularJacobianEcl(
	    const Vec3D& topocentricEquatorialAu, double observedCosDec) noexcept;

	[[nodiscard]] ORBITAL_CORE_API AstrometricResidualTransform
	astrometricResidualTransform(const PreparedObservation& observation, const Vec3D& objectHelioEclAu,
	                             const Vec3D* objectVelocityEclAuPerDay = nullptr) noexcept;

	// `sigmaCorr` is the RA/DEC correlation of the observation (ADES rmsCorr),
	// in [-1, 1]; 0 leaves the ellipse axis-aligned.  It enters before the timing
	// term, seeding the covariance with
	//   [[sRA^2, rho*sRA*sDec], [rho*sRA*sDec, sDec^2]]
	// and only then adds sigma_t^2 * v v^T.
	[[nodiscard]] AstrometricResidualTransform astrometricResidualTransformFromMotion(
	    double sigmaRaRad, double sigmaDecRad, double sigmaCorr, double timeSigmaDays, double raArcsecPerMinute,
	    double decArcsecPerMinute) noexcept;

	[[nodiscard]] ORBITAL_CORE_API AstrometricResidualPair applyAstrometricResidualTransform(
	    double raCosDecResidualRad, double decResidualRad, const AstrometricResidualTransform& transform) noexcept;

	// Capture the metric actually used by the optimizer, not a new rule lookup.
	[[nodiscard]] ORBITAL_CORE_API ResidualWeight
	residualWeightForReport(const PreparedObservation& observation, double raCosDecResidualRad, double decResidualRad,
	                        bool useWeights, const AstrometricResidualTransform& transform = {}) noexcept;

}  // namespace orbital_core::obs
