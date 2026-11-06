#include "ObserverGeometry.h"
#include <span>

#include "earth_orientation/EarthRotation.h"
#include "earth_orientation/EarthOrientation.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "orbital_core/TimeScale.h"

#include <cmath>
#include <stdexcept>

namespace orbital_core::obs {

	namespace {
		ObserverState groundObserverStateAtTimes(double mjdTt, double mjdUtc, const ObserverData& observer);

		constexpr double kEarthRotationRadPerDay = Constants::TWOPI * 1.00273790934;

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
		return groundObserverStateAtTimes(time::fromMjdTdb(mjdTdb, time::Scale::TT), mjdUtc, observer);
	}

	namespace {
		ObserverState groundObserverStateAtTimes(double mjdTt, double mjdUtc, const ObserverData& observer)
		{
			double mjdUt1 = mjdUtc;
			(void)orbital_core::earthOrientation().ttToUt1Mjd(mjdTt, mjdUt1);
			const double jdTt = time::jdFromMjd(mjdTt);
			const double gast = EarthRotation::Iau1976_1980::gmstRadiansMjd(mjdUt1) +
			                    EarthRotation::Iau1976_1980::equationOfEquinoxes(jdTt);
			const double lon = Convertors::degreesToRadians(observer.longitude);

			Vec3D obsEq(Constants::EARTH_RADIUS_IN_AU * observer.cosPhi * std::cos(lon),
			            Constants::EARTH_RADIUS_IN_AU * observer.cosPhi * std::sin(lon),
			            Constants::EARTH_RADIUS_IN_AU * observer.sinPhi);
			Vec3D velEq(0.0, 0.0, 0.0);

			{
				const double c = std::cos(gast), sg = std::sin(gast);
				const double x = obsEq.x, y = obsEq.y;
				obsEq.x = c * x - sg * y;
				obsEq.y = sg * x + c * y;
				velEq.x = -kEarthRotationRadPerDay * obsEq.y;
				velEq.y = kEarthRotationRadPerDay * obsEq.x;
				velEq.z = 0.0;
			}

			EarthRotation::Iau1976_1980::applyPrecessionNutation(jdTt, obsEq, -1);
			EarthRotation::Iau1976_1980::applyPrecessionNutation(jdTt, velEq, -1);

			return ObserverState{eqToEcl(obsEq), eqToEcl(velEq)};
		}
	}  // namespace

	Vec3D groundObserverOffsetEclAu(double jdTdb, double mjdUtc, const ObserverData& observer)
	{
		return groundObserverStateEclAu(jdTdb, mjdUtc, observer).offsetEclAu;
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
