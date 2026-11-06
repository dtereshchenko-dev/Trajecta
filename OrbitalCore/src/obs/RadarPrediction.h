#pragma once

#include "obs/PreparedObservation.h"
#include "orbital_core/IEphemerisSource.h"
#include <Eigen/Core>

namespace orbital_core::obs {

	// Two-way terrestrial radar: prepared epoch and reception state use TDB,
	// in heliocentric J2000 ecliptic AU
	// and AU/day. Delay is seconds; Doppler is received
	// minus transmitted Hz. Rows are delay then Doppler, omitting absent channels.
	struct RadarPrediction {
		double delaySeconds{};
		double dopplerHz{};
	};

	RadarPrediction predictRadar(const PreparedObservation& observation, const Vec3D& receptionPosition,
	                             const Vec3D& receptionVelocityAuPerDay, IEphemerisSource& ephemeris);

	Eigen::MatrixXd radarResidualJacobian(const PreparedObservation& observation, const Vec3D& receptionPosition,
	                                      const Vec3D& receptionVelocityAuPerDay, IEphemerisSource& ephemeris);

}  // namespace orbital_core::obs
