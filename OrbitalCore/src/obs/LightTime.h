#pragma once

#include "TwoBodyPropagation.h"
#include "orbital_core/Constants.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Types.h"
#include "orbital_core/Vec3D.h"

#include <cmath>
#include <Eigen/Core>

namespace orbital_core::obs {

	struct LightTimeSolution {
		double emissionJd{0.0};
		Vec3D position;
		Vec3D velocityAuPerDay;
		bool valid{false};
	};

	// Barycentric displacement of the Sun during the light delay, S(t - tau) - S(t),
	// in J2000 ecliptic AU. Light travels straight in the barycentric frame, so a
	// heliocentric light-time solution must add it to the emission point; without it
	// the predicted direction is off by v_sun/c (up to ~10 mas) for any distance.
	// Sources without a barycentric Sun (test doubles) return zero, which keeps the
	// purely heliocentric two-pass solution.
	[[nodiscard]] inline Vec3D solarDisplacementDuringDelay(double receptionJd, double delayDays,
	                                                        Vec3D* emissionSunVelocity = nullptr)
	{
		if (emissionSunVelocity)
			*emissionSunVelocity = {};
		IEphemerisSource* ephemeris = getDefaultEphemerisSource();
		if (!ephemeris)
			return {};
		const double receptionMjd = receptionJd - 2400000.5;
		const auto reception = ephemeris->barycentricSunMjd(receptionMjd);
		const auto emission = ephemeris->barycentricSunMjd(receptionMjd - delayDays);
		if (!reception || !emission)
			return {};
		if (emissionSunVelocity)
			*emissionSunVelocity = emission->V;
		return emission->R - reception->R;
	}

	// Two-pass light time (estimate rho/c at reception, propagate the
	// asteroid two-body over that delay, recompute rho/c, propagate again from the
	// reception state), but measured in the barycentric frame: the returned position
	// is the heliocentric emission point plus the solar displacement during the
	// delay, so position - observer is the barycentric emission-to-reception ray.
	// The velocity stays the heliocentric emission velocity.
	[[nodiscard]] inline LightTimeSolution solveLightTimeTwoPass(
	    double receptionJd, const Vec3D& receptionPosition, const Vec3D& receptionVelocityAuPerDay,
	    const Vec3D& observerHeliocentricPosition, Eigen::Matrix<double, 3, 6>* positionJacobian = nullptr)
	{
		LightTimeSolution out;
		const double firstDelayDays =
		    (receptionPosition - observerHeliocentricPosition).magnitude() / Constants::AU_PER_DAY;
		if (!std::isfinite(firstDelayDays))
			return out;

		const StateVectors receptionState(receptionPosition, receptionVelocityAuPerDay / Constants::GAUSS_K);
		const StateVectors firstEmission = Calculations::propagateKeplerState(receptionState, -firstDelayDays);
		if (!firstEmission.isValid())
			return out;

		Vec3D firstSunVelocity;
		const Vec3D firstShift = solarDisplacementDuringDelay(receptionJd, firstDelayDays, &firstSunVelocity);
		const double refinedDelayDays =
		    (firstEmission.first + firstShift - observerHeliocentricPosition).magnitude() / Constants::AU_PER_DAY;
		if (!std::isfinite(refinedDelayDays))
			return out;

		const StateVectors finalEmission = Calculations::propagateKeplerState(receptionState, -refinedDelayDays);
		if (!finalEmission.isValid())
			return out;

		out.emissionJd = receptionJd - refinedDelayDays;
		Vec3D finalSunVelocity;
		out.position =
		    finalEmission.first + solarDisplacementDuringDelay(receptionJd, refinedDelayDays, &finalSunVelocity);
		out.velocityAuPerDay = finalEmission.second * Constants::GAUSS_K;
		if (positionJacobian) {
			if (!(firstDelayDays > 0) || !(refinedDelayDays > 0))
				return out;
			const Vec3D firstDirection =
			    (receptionPosition - observerHeliocentricPosition) / (firstDelayDays * Constants::AU_PER_DAY);
			const Vec3D refinedDirection = (firstEmission.first + firstShift - observerHeliocentricPosition) /
			                               (refinedDelayDays * Constants::AU_PER_DAY);
			// Differentiate both delay estimates, the Kepler maps and the moving
			// heliocentric origin. Columns use reception [AU, AU/day].
			for (int column = 0; column < 6; ++column) {
				const Vec3D unit(column % 3 == 0 ? 1.0 : 0.0, column % 3 == 1 ? 1.0 : 0.0, column % 3 == 2 ? 1.0 : 0.0);
				const StateVectors tangent(column < 3 ? unit : Vec3D{},
				                           column >= 3 ? unit / Constants::GAUSS_K : Vec3D{});
				const double firstDelayRate = firstDirection.dot(tangent.first) / Constants::AU_PER_DAY;
				const Vec3D firstRate = Calculations::propagateKeplerPositionRate(receptionState, -firstDelayDays,
				                                                                  tangent, -firstDelayRate) -
				                        firstSunVelocity * firstDelayRate;
				const double refinedDelayRate = refinedDirection.dot(firstRate) / Constants::AU_PER_DAY;
				const Vec3D finalRate = Calculations::propagateKeplerPositionRate(receptionState, -refinedDelayDays,
				                                                                  tangent, -refinedDelayRate) -
				                        finalSunVelocity * refinedDelayRate;
				positionJacobian->col(column) << finalRate.x, finalRate.y, finalRate.z;
			}
		}
		out.valid = true;
		return out;
	}

}  // namespace orbital_core::obs
