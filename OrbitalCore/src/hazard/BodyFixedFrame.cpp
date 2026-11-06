#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include "earth_orientation/EarthRotation.h"
#include "earth_orientation/EarthOrientation.h"
#include "IauRotationData.h"
#include <Eigen/Core>
#include <cmath>
#include <limits>

namespace orbital_core::hazard {

	std::string earthBodyFixedFrameModel()
	{
		const char* profile = "unknown";
		switch (orbital_core::earthOrientation().profile()) {
			case EarthOrientationProfile::IAU1980:
				profile = "IAU1980";
				break;
			case EarthOrientationProfile::IAU2000:
				profile = "IAU2000";
				break;
			case EarthOrientationProfile::FullC04:
				profile = "full-C04";
				break;
		}
		return std::string(
		           "IAU1976/1980; UT1/GAST (IAU1982, finals2000A with UTC fallback); "
		           "polar motion: not applied; EOP profile: ") +
		       profile;
	}

	namespace {
		Vec3D icrfToEcefEarth(const Vec3D& rGeoEclKm, double mjdTdb) noexcept
		{
			const double mjdTt = time::fromMjdTdb(mjdTdb, time::Scale::TT);
			const double jdTt = time::jdFromMjd(mjdTt);
			Vec3D rEq = Helpers::eclToEcv(rGeoEclKm);
			EarthRotation::Iau1976_1980::applyPrecessionNutation(jdTt, rEq, +1);
			double mjdUt1 = time::fromMjdTdb(mjdTdb, time::Scale::UTC);
			(void)orbital_core::earthOrientation().ttToUt1Mjd(mjdTt, mjdUt1);
			const double GAST_rad = EarthRotation::Iau1976_1980::gmstRadiansMjd(mjdUt1) +
			                        EarthRotation::Iau1976_1980::equationOfEquinoxes(jdTt);
			const double c = std::cos(-GAST_rad);
			const double s = std::sin(-GAST_rad);
			return Vec3D(c * rEq.x - s * rEq.y, s * rEq.x + c * rEq.y, rEq.z);
		}

		Eigen::Matrix3d rotationZ(double angle)
		{
			const double c = std::cos(angle), s = std::sin(angle);
			Eigen::Matrix3d matrix;
			matrix << c, s, 0, -s, c, 0, 0, 0, 1;
			return matrix;
		}

		Eigen::Matrix3d rotationX(double angle)
		{
			const double c = std::cos(angle), s = std::sin(angle);
			Eigen::Matrix3d matrix;
			matrix << 1, 0, 0, 0, c, s, 0, -s, c;
			return matrix;
		}

		BodyFixedTransform transform(const Eigen::Matrix3d& matrix, std::string_view frame, std::string_view model)
		{
			BodyFixedTransform result;
			for (int i = 0; i < 3; ++i)
				result.equatorialRows[i] = {matrix(i, 0), matrix(i, 1), matrix(i, 2)};
			result.frame = frame;
			result.model = model;
			result.valid = matrix.allFinite();
			return result;
		}

		double polynomial(const orientation::Polynomial& coefficients, double t)
		{
			return (coefficients[2] * t + coefficients[1]) * t + coefficients[0];
		}

		double radians(double degrees) { return std::remainder(degrees, 360.0) * Constants::DEG_TO_RAD; }
	}  // namespace

	BodyFixedTransform bodyFixedTransformMjd(Planet planet, double mjdTdb, IEphemerisSource* source) noexcept
	{
		if (!std::isfinite(mjdTdb))
			return {};
		try {
			if (planet == Planet::Earth) {
				Eigen::Matrix3d matrix;
				for (int i = 0; i < 3; ++i) {
					const Vec3D basis{i == 0 ? 1.0 : 0.0, i == 1 ? 1.0 : 0.0, i == 2 ? 1.0 : 0.0};
					const auto column = icrfToEcefEarth(Helpers::ecvToEcl(basis), mjdTdb);
					matrix.col(i) << column.x, column.y, column.z;
				}
				return transform(matrix, "EARTH_OPTICAL_UT1", "IAU1976/1980; UT1/GAST; no polar motion");
			}
			if (planet == Planet::Moon) {
				if (!source)
					source = getDefaultEphemerisSource();
				const auto angles = source ? source->lunarLibrationMjdTdb(mjdTdb) : std::nullopt;
				if (!angles || angles->ephemerisVersion != 440)
					return {};
				const auto& e = angles->angles;
				if (!std::isfinite(e.x) || !std::isfinite(e.y) || !std::isfinite(e.z))
					return {};
				Eigen::Matrix3d paToMe;
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						paToMe(i, j) = orientation::kMoonPaToMe[i][j];
				const Eigen::Matrix3d matrix = paToMe * rotationZ(e.z) * rotationX(e.y) * rotationZ(e.x);
				return transform(matrix, "MOON_ME_DE440_ME421", "DE440 librations; NAIF moon_de440_250416");
			}
			const double d = mjdTdb - 51544.5;
			const double t = d / Constants::JULIAN_CENTURY_DAYS;
			for (const auto& model : orientation::kModels) {
				if (model.planet != planet)
					continue;
				double ra = polynomial(model.ra, t), dec = polynomial(model.dec, t);
				double pm = polynomial(model.pm, d);
				for (const auto& term : model.periodic) {
					const double phase = radians(polynomial(term.phase, t));
					ra += term.ra * std::sin(phase);
					dec += term.dec * std::cos(phase);
					pm += term.pm * std::sin(phase);
				}
				const Eigen::Matrix3d matrix =
				    rotationZ(radians(pm)) * rotationX(radians(90.0 - dec)) * rotationZ(radians(ra + 90.0));
				return transform(matrix, model.frame, "IAU WGCCRE; NAIF pck00011");
			}
		} catch (...) {
			return {};
		}
		return {};
	}

	Vec3D icrfToBodyFixedMjd(Planet planet, const Vec3D& r, double mjdTdb, IEphemerisSource* source) noexcept
	{
		if (planet == Planet::Earth && std::isfinite(mjdTdb))
			return icrfToEcefEarth(r, mjdTdb);
		const auto rotation = bodyFixedTransformMjd(planet, mjdTdb, source);
		if (rotation.valid)
			return rotation.applyEquatorial(Helpers::eclToEcv(r));
		const double nan = std::numeric_limits<double>::quiet_NaN();
		return {nan, nan, nan};
	}

	Vec3D icrfToBodyFixed(Planet planet, const Vec3D& rPlanetCentricEclKm, double jdTdb) noexcept
	{
		return icrfToBodyFixedMjd(planet, rPlanetCentricEclKm, time::mjdFromJd(jdTdb));
	}

	LatLon latLonFromBodyFixed(const Vec3D& rBodyFixedKm) noexcept
	{
		LatLon ll;
		const double rMag = rBodyFixedKm.magnitude();
		if (!std::isfinite(rMag) || rMag < 1.0e-6)
			return ll;
		ll.latitudeDeg = std::asin(rBodyFixedKm.z / rMag) * Constants::RAD_TO_DEG;
		const double lonRad = std::atan2(rBodyFixedKm.y, rBodyFixedKm.x);
		ll.longitudeDeg = Helpers::wrapPi(lonRad) * Constants::RAD_TO_DEG;
		ll.valid = true;
		return ll;
	}

	LatLon latLonFromBodyFixed(Planet planet, const Vec3D& rBodyFixedKm) noexcept
	{
		if (planet == Planet::Count_ || !std::isfinite(rBodyFixedKm.dot(rBodyFixedKm)))
			return {};
		const double flattening = planetFlattening(planet);
		if (flattening <= 0.0)
			return latLonFromBodyFixed(rBodyFixedKm);
		const SurfaceGeodetic g = geodeticOnReferenceEllipsoidKm(rBodyFixedKm, planetRadiusKm(planet), flattening);
		if (!g.valid)
			return {};
		return {g.latitudeRad * Constants::RAD_TO_DEG, Helpers::wrapPi(g.longitudeRad) * Constants::RAD_TO_DEG, true};
	}
}  // namespace orbital_core::hazard
