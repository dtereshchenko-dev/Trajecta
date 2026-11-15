#pragma once

#include "orbital_core/Ephemeris.h"
#include "orbital_core/IEphemerisSource.h"

#include <optional>

namespace orbital_core {
	// Geometric states at one epoch; no observational light-time correction.
	// Rotating axes: x away from the Sun, z along the planet's orbital momentum.
	struct RelativeMotionSample {
		double mjdTdb{};
		Vec3D positionAu;
		Vec3D velocityAuPerDay;
		Vec3D rotatingPositionAu;
		Vec3D moonPositionAu;
		Vec3D rotatingMoonPositionAu;
		Vec3D sunDirection;
		double distanceAu{};
		double hillRadiusAu{};
		double specificEnergyAu2PerDay2{};
		bool hasMoon{};
	};

	struct RelativeMotionRequest {
		EphemerisSolution solution;
		EphemerisPlanet planet{EphemerisPlanet::Earth};
		std::vector<double> epochsMjdTdb;
	};

	struct TemporaryCaptureInterval {
		double startMjdTdb{};
		double endMjdTdb{};
		double signedRevolutions{};
		double inertialSignedRevolutions{};
		double angularTravelRevolutions{};
		double inertialAngularTravelRevolutions{};
		bool startClipped{};
		bool endClipped{};
		bool revolutionsValid{true};
		bool inertialRevolutionsValid{true};
	};

	struct CaptureTimeUncertainty {
		std::optional<double> start3SigmaDays;
		std::optional<double> end3SigmaDays;
	};

	// The Hill sphere and two-body energy are diagnostics, not capture decisions.
	[[nodiscard]] ORBITAL_CORE_API std::vector<RelativeMotionSample> computeRelativeMotion(
	    const RelativeMotionRequest& request, const FitRunControl& control = {});

	// Adds propagated samples where the body travels rapidly relative to the planet.
	// The requested epochs are retained; this only refines the displayed trajectory.
	[[nodiscard]] ORBITAL_CORE_API std::vector<RelativeMotionSample> computeRefinedRelativeMotion(
	    const RelativeMotionRequest& request, const FitRunControl& control = {});

	// Sampled intervals where planetocentric two-body energy is negative and distance < 3 Hill radii.
	// Boundaries and positions are linearly interpolated between samples.
	[[nodiscard]] ORBITAL_CORE_API std::vector<TemporaryCaptureInterval> findTemporaryCaptures(
	    const std::vector<RelativeMotionSample>& samples);

	// Linearized uncertainty of the interpolated event boundaries, using the full fit covariance.
	[[nodiscard]] ORBITAL_CORE_API std::vector<CaptureTimeUncertainty> captureTimeUncertainty(
	    const RelativeMotionRequest& request, const std::vector<RelativeMotionSample>& samples,
	    const std::vector<TemporaryCaptureInterval>& captures, const FitRunControl& control = {});

}  // namespace orbital_core
