#pragma once

#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Types.h"
#include "orbital_core/Vec3D.h"

namespace orbital_core::obs {

	// Barycentric displacement of the Sun during the light delay, S(t - tau) - S(t),
	// in J2000 ecliptic AU. Light travels straight in the barycentric frame, so a
	// heliocentric light-time solution must add it to the emission point; without it
	// the predicted direction is off by v_sun/c (up to ~10 mas) for any distance.
	// Sources without a barycentric Sun return zero.
	[[nodiscard]] inline Vec3D solarDisplacementDuringDelay(double receptionJd, double delayDays)
	{
		IEphemerisSource* ephemeris = getDefaultEphemerisSource();
		if (!ephemeris)
			return {};
		const double receptionMjd = receptionJd - 2400000.5;
		const auto reception = ephemeris->barycentricSunMjd(receptionMjd);
		const auto emission = ephemeris->barycentricSunMjd(receptionMjd - delayDays);
		if (!reception || !emission)
			return {};
		return emission->R - reception->R;
	}

}  // namespace orbital_core::obs
