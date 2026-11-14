#include "ObserverGeometry.h"
#include <span>

#include "earth_orientation/EarthRotation.h"
#include "earth_orientation/EarthOrientation.h"
#include "earth_orientation/EarthPck.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "orbital_core/TimeScale.h"

#include <cmath>
#include <stdexcept>

namespace orbital_core::obs {

	namespace {
		Vec3D groundObserverPositionAtTimes(double mjdTt, double mjdUtc, const ObserverData& observer,
		                                    const EarthOrientation& orientation);

		[[nodiscard]] Vec3D eqToEcl(const Vec3D& eq) noexcept
		{
			const double cosE = std::cos(Constants::EARTH_INCL);
			const double sinE = std::sin(Constants::EARTH_INCL);
			return Vec3D(eq.x, eq.y * cosE + eq.z * sinE, -eq.y * sinE + eq.z * cosE);
		}

		[[nodiscard]] Vec3D eclToEq(const Vec3D& ecl) noexcept
		{
			const double cosE = std::cos(Constants::EARTH_INCL);
			const double sinE = std::sin(Constants::EARTH_INCL);
			return Vec3D(ecl.x, ecl.y * cosE - ecl.z * sinE, ecl.y * sinE + ecl.z * cosE);
		}

		[[nodiscard]] Vec3D terrestrialObserver(const ObserverData& observer)
		{
			const double lon = Convertors::degreesToRadians(observer.longitude);
			return {Constants::EARTH_RADIUS_IN_AU * observer.cosPhi * std::cos(lon),
			        Constants::EARTH_RADIUS_IN_AU * observer.cosPhi * std::sin(lon),
			        Constants::EARTH_RADIUS_IN_AU * observer.sinPhi};
		}
	}  // namespace

	ObserverData rovingObserverData(const RovingObserver& position, const std::string& code)
	{
		if (!position.valid())
			throw std::invalid_argument("Invalid terrestrial observer position");
		constexpr double e2 = Constants::EARTH_FLATTENING * (2.0 - Constants::EARTH_FLATTENING);
		const double latitude = position.latitudeDeg * Constants::DEG_TO_RAD;
		const double sinLat = std::sin(latitude);
		const double radiusKm = Constants::EARTH_EQ_RADIUS_KM / std::sqrt(1.0 - e2 * sinLat * sinLat);
		const double heightKm = position.heightMeters / 1000.0;
		return {position.longitudeDeg, (radiusKm + heightKm) * std::cos(latitude) / Constants::EARTH_PARALLAX_RADIUS_KM,
		        (radiusKm * (1.0 - e2) + heightKm) * sinLat / Constants::EARTH_PARALLAX_RADIUS_KM, code};
	}

	ObserverState groundObserverStateEclAu(double jdTdb, double mjdUtc, const ObserverData& observer)
	{
		return groundObserverStateEclAuMjd(time::mjdFromJd(jdTdb), mjdUtc, observer);
	}

	ObserverState groundObserverStateEclAuMjd(double mjdTdb, double mjdUtc, const ObserverData& observer)
	{
		if (observer.cosPhi == 0.0 && observer.sinPhi == 0.0)
			return {};
		Vec3D position, velocity;
		if (earthPckState(mjdTdb, terrestrialObserver(observer), position, velocity))
			return {eqToEcl(position), eqToEcl(velocity)};
		double ut1;
		if (earthOrientation().ttToUt1Mjd(time::fromMjdTdb(mjdTdb, time::Scale::TT), ut1))
			return groundObserverStateEclAuMjd(mjdTdb, mjdUtc, observer, earthOrientation());
		throw std::runtime_error(
		    "Observer prediction requires Earth PCK or observed Earth orientation data at observation epoch");
	}

	ObserverState radarStationStateEclAuMjd(double mjdTdb, double mjdUtc, const ObserverData& observer)
	{
		Vec3D position, velocity;
		if (earthPckState(mjdTdb, terrestrialObserver(observer), position, velocity))
			return {eqToEcl(position), eqToEcl(velocity)};
		constexpr double h = 1.0 / Constants::SECONDS_PER_DAY;
		const auto covered = [&](const EarthOrientation& orientation) {
			for (const double t : {mjdTdb - h, mjdTdb, mjdTdb + h})
				if (double ut1; !orientation.ttToUt1Mjd(time::fromMjdTdb(t, time::Scale::TT), ut1))
					return false;
			return true;
		};
		if (covered(earthOrientation()))
			return groundObserverStateEclAuMjd(mjdTdb, mjdUtc, observer, earthOrientation());
		throw std::runtime_error("Radar prediction requires Earth PCK or Earth orientation data at observation epoch");
	}

	ObserverState groundObserverStateEclAuMjd(double mjdTdb, double mjdUtc, const ObserverData& observer,
	                                          const EarthOrientation& orientation)
	{
		const double mjdTt = time::fromMjdTdb(mjdTdb, time::Scale::TT);
		constexpr double h = 1.0 / Constants::SECONDS_PER_DAY;
		const Vec3D position = groundObserverPositionAtTimes(mjdTt, mjdUtc, observer, orientation);
		const double beforeTdb = mjdTdb - h, afterTdb = mjdTdb + h;
		const Vec3D before = groundObserverPositionAtTimes(time::fromMjdTdb(beforeTdb, time::Scale::TT), mjdUtc - h,
		                                                   observer, orientation);
		const Vec3D after = groundObserverPositionAtTimes(time::fromMjdTdb(afterTdb, time::Scale::TT), mjdUtc + h,
		                                                  observer, orientation);
		return {position, (after - before) / (afterTdb - beforeTdb)};
	}

	namespace {
		Vec3D groundObserverPositionAtTimes(double mjdTt, double mjdUtc, const ObserverData& observer,
		                                    const EarthOrientation& orientation)
		{
			auto eop = orientation.lookupMjdTt(mjdTt);
			if (!orientation.coversObservedMjdTt(mjdTt)) {
				double mjdUt1;
				if (!orientation.ttToUt1Mjd(mjdTt, mjdUt1))
					throw std::runtime_error("Earth orientation and UT1 data unavailable at observation epoch");
				eop = {};
				eop.ttMinusUt1Seconds = (mjdTt - mjdUt1) * Constants::SECONDS_PER_DAY;
			}
			const Vec3D body = terrestrialObserver(observer);
			return eqToEcl(EarthRotation::WithEop::terrestrialToEquatorialMjd(mjdTt, body, eop));
		}
	}  // namespace

	Vec3D groundObserverOffsetEclAu(double jdTdb, double mjdUtc, const ObserverData& observer)
	{
		return groundObserverOffsetEclAuMjd(time::mjdFromJd(jdTdb), mjdUtc, observer);
	}

	Vec3D groundObserverOffsetEclAuMjd(double mjdTdb, double mjdUtc, const ObserverData& observer)
	{
		return groundObserverStateEclAuMjd(mjdTdb, mjdUtc, observer).offsetEclAu;
	}

	Vec3D observerOffsetAt(double jdTdb, const PreparedObservation& observation)
	{
		if (observation.isSatellite)
			return observation.observerOffsetEclAu;

		const double mjdUtc = observation.mjdUtc + (jdTdb - observation.jdTdb);
		return groundObserverOffsetEclAu(jdTdb, mjdUtc, observation.observer);
	}

	ApparentMotion apparentMotionArcsecPerMinute(const Vec3D& objectHelioEclAu, const Vec3D& objectVelocityEclAuPerDay,
	                                             const PreparedObservation& observation)
	{
		const Vec3D topoEcl = objectHelioEclAu - observation.earthR - observation.observerOffsetEclAu;
		const Vec3D relVelEcl =
		    objectVelocityEclAuPerDay - observation.earthV - observation.observerVelocityEclAuPerDay;

		const Vec3D topoEq = eclToEq(topoEcl);
		const Vec3D relVelEq = eclToEq(relVelEcl);
		const double rho = topoEq.magnitude();
		const double xy = std::hypot(topoEq.x, topoEq.y);
		if (rho <= 0.0 || xy <= 0.0 || !std::isfinite(rho) || !std::isfinite(xy)) {
			return {};
		}

		const double invRho = 1.0 / rho;
		const double invXy = 1.0 / xy;
		const Vec3D east(-topoEq.y * invXy, topoEq.x * invXy, 0.0);
		const Vec3D north(-topoEq.x * topoEq.z * invXy * invRho, -topoEq.y * topoEq.z * invXy * invRho, xy * invRho);

		constexpr double kRadPerDayToArcsecPerMinute = 1.0 / (Constants::ARCSEC_RAD * 1440.0);
		return ApparentMotion{east.dot(relVelEq) * invRho * kRadPerDayToArcsecPerMinute,
		                      north.dot(relVelEq) * invRho * kRadPerDayToArcsecPerMinute, true};
	}

}  // namespace orbital_core::obs
