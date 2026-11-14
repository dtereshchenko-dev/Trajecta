#include "obs/RadarPrediction.h"
#include "obs/ObserverGeometry.h"
#include "TwoBodyPropagation.h"
#include "earth_orientation/EarthRotation.h"
#include "earth_orientation/EarthOrientation.h"
#include "orbital_core/Types.h"
#include "orbital_core/TimeScale.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core::obs {
	namespace {
		constexpr double c = Constants::AU_PER_DAY;
		constexpr double daySeconds = Constants::SECONDS_PER_DAY;

		struct StationState {
			Vec3D r, v, ground, groundVelocity;
		};

		using Calculations::tdbMinusTtSeconds;

		StationState stationState(double mjdTdb, const ObserverData& station, IEphemerisSource& eph)
		{
			const auto earth = eph.queryMjd(EphemerisPlanet::Earth, mjdTdb);
			const auto ground = radarStationStateEclAuMjd(mjdTdb, time::fromMjdTdb(mjdTdb, time::Scale::UTC), station);
			return {earth.R + ground.offsetEclAu, earth.V + ground.velocityEclAuPerDay, ground.offsetEclAu,
			        ground.velocityEclAuPerDay};
		}

		double shapiro(const Vec3D& object, const Vec3D& station)
		{
			const double sum = object.magnitude() + station.magnitude();
			const double distance = (object - station).magnitude();
			if (!(sum > distance))
				throw std::runtime_error("Invalid radar Shapiro geometry");
			return 2.0 * Constants::GM_SUN / (c * c * c) * std::log1p(2 * distance / (sum - distance));
		}

		double tropo(const Vec3D& object, const StationState& station)
		{
			const auto ray = object - station.r;
			const double cosz =
			    std::clamp(ray.dot(station.ground) / (ray.magnitude() * station.ground.magnitude()), -1.0, 1.0);
			const double sinz = std::sqrt(std::max(0.0, 1 - cosz * cosz));
			const double mapping = cosz + 0.0014 * sinz / (0.045 * sinz + cosz);
			return 7e-9 / daySeconds / mapping;
		}

		double propagationDelay(const Vec3D& object, const StationState& station, double radiusAu)
		{
			return ((object - station.r).magnitude() - radiusAu) / c + shapiro(object, station.r) +
			       tropo(object, station);
		}

		double delayRateCorrection(const Vec3D& object, const Vec3D& velocity, const StationState& station)
		{
			constexpr double h = 0.25 / daySeconds;
			const auto shifted = [&](double dt) {
				const StationState s{station.r + station.v * dt, station.v,
				                     station.ground + station.groundVelocity * dt, station.groundVelocity};
				return shapiro(object + velocity * dt, s.r) + tropo(object + velocity * dt, s);
			};
			return c * (shifted(h) - shifted(-h)) / (2 * h);
		}
	}  // namespace

	RadarPrediction predictRadar(const PreparedObservation& observation, const Vec3D& receptionPosition,
	                             const Vec3D& receptionVelocityAuPerDay, IEphemerisSource& ephemeris)
	{
		if (!observation.radar || !observation.radar->valid() || !ephemeris.isAvailable())
			throw std::invalid_argument("Radar prediction requires a valid measurement and ephemeris");
		const auto& radar = *observation.radar;
		const double radius = radar.returnType == RadarReturn::Surface
		                          ? radar.surfaceRadiusKm.value_or(-1.0) * 1000.0 / Constants::AU_IN_METERS
		                          : 0.0;
		if (radius < 0.0)
			throw std::invalid_argument("Surface radar prediction requires a body radius");
		const double reception = observation.mjdTdb;
		const auto receiver = stationState(reception, observation.observer, ephemeris);
		const StateVectors initial(receptionPosition, receptionVelocityAuPerDay / Constants::GAUSS_K);
		double down = (receptionPosition - receiver.r).magnitude() / c;
		StateVectors bounce;
		bool converged = false;
		for (int iteration = 0; iteration < 16; ++iteration) {
			bounce = Calculations::propagateKeplerState(initial, -down);
			if (!bounce.isValid())
				throw std::runtime_error("Radar bounce propagation failed");
			const double next = propagationDelay(bounce.first, receiver, radius);
			converged = std::abs(next - down) < 1e-15;
			down = next;
			if (converged)
				break;
		}
		if (!converged)
			throw std::runtime_error("Radar down-leg light time did not converge");
		bounce = Calculations::propagateKeplerState(initial, -down);
		double up = down;
		StationState transmitter;
		converged = false;
		for (int iteration = 0; iteration < 16; ++iteration) {
			const double transmitTdb = reception - down - up;
			transmitter = stationState(transmitTdb, observation.transmitter, ephemeris);
			const double next = propagationDelay(bounce.first, transmitter, radius);
			converged = std::abs(next - up) < 1e-15;
			up = next;
			if (converged)
				break;
		}
		if (!converged)
			throw std::runtime_error("Radar up-leg light time did not converge");
		const Vec3D velocity = bounce.second * Constants::GAUSS_K;
		const auto sunR = ephemeris.barycentricSunMjd(reception);
		const auto sunB = ephemeris.barycentricSunMjd(reception - down);
		const auto sunT = ephemeris.barycentricSunMjd(reception - down - up);
		if (!sunR || !sunB || !sunT)
			throw std::runtime_error("Radar Doppler requires barycentric Sun ephemerides");
		const auto rayD = bounce.first - receiver.r;
		const auto rayU = bounce.first - transmitter.r;
		const auto unitD = rayD / rayD.magnitude(), unitU = rayU / rayU.magnitude();
		const double rd = unitD.dot(velocity - receiver.v);
		const double ru = unitU.dot(velocity - transmitter.v);
		const auto txV = transmitter.v + sunT->V, rxV = receiver.v + sunR->V;
		const double s1 = unitU.dot(txV) / c, s2 = unitD.dot(velocity + sunB->V) / c;
		const double s3 = 0.5 * (txV.dot(txV) - rxV.dot(rxV)) +
		                  Constants::GM_SUN * (1.0 / transmitter.r.magnitude() - 1.0 / receiver.r.magnitude());
		const double level = ru + rd;
		double twoWayRate =
		    level + ru * s1 * (1 + s1) - rd * s2 * (1 - s2) - (ru * rd * (1 + s1 - s2) - s3 * (1 - level / c)) / c;
		twoWayRate += delayRateCorrection(bounce.first, velocity, receiver) +
		              delayRateCorrection(bounce.first, velocity, transmitter);
		const double clockCorrection = tdbMinusTtSeconds(time::fromMjdTdb(reception - down - up, time::Scale::TT)) -
		                               tdbMinusTtSeconds(time::fromMjdTdb(reception, time::Scale::TT));
		const RadarPrediction out{(down + up) * daySeconds + clockCorrection, -radar.frequencyHz * twoWayRate / c};
		if (!std::isfinite(out.delaySeconds) || !std::isfinite(out.dopplerHz))
			throw std::runtime_error("Non-finite radar prediction");
		return out;
	}

	Eigen::MatrixXd radarResidualJacobian(const PreparedObservation& observation, const Vec3D& position,
	                                      const Vec3D& velocity, IEphemerisSource& ephemeris)
	{
		Eigen::MatrixXd jac(observation.measurementDimensions(), 6);
		for (int column = 0; column < 6; ++column) {
			Vec3D rp = position, rm = position, vp = velocity, vm = velocity;
			const double h = column < 3 ? 1e-7 : 1e-8;
			const Vec3D delta(column % 3 == 0 ? h : 0, column % 3 == 1 ? h : 0, column % 3 == 2 ? h : 0);
			if (column < 3) {
				rp = rp + delta;
				rm = rm - delta;
			} else {
				vp = vp + delta;
				vm = vm - delta;
			}
			const auto plus = predictRadar(observation, rp, vp, ephemeris);
			const auto minus = predictRadar(observation, rm, vm, ephemeris);
			int row = 0;
			if (observation.radar->delaySeconds)
				jac(row++, column) = -(plus.delaySeconds - minus.delaySeconds) / (2 * h);
			if (observation.radar->dopplerHz)
				jac(row, column) = -(plus.dopplerHz - minus.dopplerHz) / (2 * h);
		}
		return jac;
	}
}  // namespace orbital_core::obs
