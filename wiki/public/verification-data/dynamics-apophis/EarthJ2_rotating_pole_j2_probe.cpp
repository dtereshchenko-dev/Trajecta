#include "orbital_core/forces/EarthJ2.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/Types.h"
#include "earth_orientation/EarthOrientation.h"
#include "earth_orientation/EarthRotation.h"
#include <Eigen/Geometry>

unsigned long long g_probeEarthJ2Evaluations = 0;

namespace orbital_core::forces::detail {
	struct EarthZonalParameters {
		double radiusAu;
		double j2;
		double j3;
		double j4;
	};

	inline constexpr EarthZonalParameters kDE440EarthZonals{6378.1366 / 149597870.7, 1.08262539e-3, 0.0, 0.0};

	inline constexpr double kDE440EarthGm = 8.9970113929473466e-10 * (81.300568221497215 / (1.0 + 81.300568221497215));

	struct EarthZonalResult {
		Eigen::Vector3d a = Eigen::Vector3d::Zero();
		Eigen::Matrix3d dadr = Eigen::Matrix3d::Zero();
	};

	[[nodiscard]] inline EarthZonalResult evaluateEarthZonals(const Eigen::Vector3d& rEq, double gm,
	                                                          const EarthZonalParameters& c,
	                                                          bool computePartials) noexcept
	{
		EarthZonalResult out;

		const double r2 = rEq.squaredNorm();
		if (r2 < 1e-24)
			return out;

		const double invR = 1.0 / std::sqrt(r2);
		const Eigen::Vector3d u = rEq * invR;

		const double s = u.z();
		const double s2 = s * s;
		const double s4 = s2 * s2;
		const double rho = c.radiusAu * invR;
		const double rho2 = rho * rho;
		const double g = gm / r2;

		const double q2 = g * c.j2 * rho2;
		const double q3 = g * c.j3 * rho2 * rho;
		const double q4 = g * c.j4 * rho2 * rho2;
		const double A2 = 1.5 * (5.0 * s2 - 1.0);
		const double B2 = 3.0 * s;
		const double A3 = 2.5 * s * (7.0 * s2 - 3.0);
		const double B3 = 1.5 * (5.0 * s2 - 1.0);
		const double A4 = 0.625 * (63.0 * s4 - 42.0 * s2 + 3.0);
		const double B4 = 2.5 * s * (7.0 * s2 - 3.0);

		const double radial = q2 * A2 + q3 * A3 + q4 * A4;
		const double axial = q2 * B2 + q3 * B3 + q4 * B4;

		out.a = radial * u;
		out.a.z() -= axial;

		if (!computePartials)
			return out;

		const double D2 = 15.0 * s;
		const double D3 = 7.5 * (7.0 * s2 - 1.0);
		const double D4 = 52.5 * s * (3.0 * s2 - 1.0);
		const double uu = -q2 * (5.0 * A2 + s * D2) - q3 * (6.0 * A3 + s * D3) - q4 * (7.0 * A4 + s * D4);
		const double uz = q2 * D2 + q3 * D3 + q4 * D4;
		const double zz = 3.0 * q2 + 15.0 * s * q3 + 7.5 * (7.0 * s2 - 1.0) * q4;

		out.dadr = uu * (u * u.transpose());
		out.dadr.diagonal().array() += radial;

		out.dadr.col(2) += uz * u;
		out.dadr.row(2) += uz * u.transpose();

		out.dadr(2, 2) -= zz;
		out.dadr *= invR;

		return out;
	}

}  // namespace orbital_core::forces::detail

namespace orbital_core::forces {
	void EarthJ2::evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d&) const noexcept
	{
		constexpr uint16_t kEarthBodyMask = static_cast<uint16_t>(1u << 2);
		if ((ctx.bodyMask & kEarthBodyMask) == 0 || !ctx.ephemerisSource || !ctx.ephemerisSource->isAvailable())
			return;

		const Vec3D rEarth = ctx.havePlanets
		                         ? ctx.planetPositions[2]
		                         : ctx.ephemerisSource->getPlanetPositionMjd(EphemerisPlanet::Earth, ctx.mjd);

		const Vec3D dGeoEcl = ctx.R - rEarth;
		const double rGeo = dGeoEcl.magnitude();
		if (rGeo >= Constants::EARTH_LOCAL_DYNAMICS_AU || rGeo < 1e-12)
			return;
		++g_probeEarthJ2Evaluations;

		const Eigen::Matrix3d rot = equatorialToEclipticAtMjd(ctx.mjd);
		const Eigen::Vector3d rEq = rot.transpose() * Eigen::Vector3d(dGeoEcl.x, dGeoEcl.y, dGeoEcl.z);

		const double gm = ctx.earthGm > 0.0 ? ctx.earthGm : detail::kDE440EarthGm;
		const auto zonals = detail::evaluateEarthZonals(rEq, gm, detail::kDE440EarthZonals, ctx.computePartials);

		const Eigen::Vector3d aEcl = rot * zonals.a;

		a.x += aEcl.x();
		a.y += aEcl.y();
		a.z += aEcl.z();

		if (ctx.computePartials)
			dadr += rot * zonals.dadr * rot.transpose();
	}

	Eigen::Matrix3d EarthJ2::equatorialToEclipticAtMjd(double mjdTdb) noexcept
	{
		thread_local bool cached = false;
		thread_local double cachedMjd = 0.0;
		thread_local Eigen::Matrix3d rotation;
		if (!cached || mjdTdb != cachedMjd) {
			const double mjdTt = time::fromMjdTdb(mjdTdb, time::Scale::TT);
			const auto eop = earthOrientation().lookupMjdTt(mjdTt);
			const Vec3D pole = Helpers::ecvToEcl(EarthRotation::WithEop::poleEquatorialMjd(mjdTt, eop));
			const Eigen::Vector3d z = Eigen::Vector3d(pole.x, pole.y, pole.z).normalized();
			const Eigen::Vector3d x = (Eigen::Vector3d::UnitX() - z.x() * z).normalized();
			rotation.col(0) = x;
			rotation.col(1) = z.cross(x);
			rotation.col(2) = z;
			cachedMjd = mjdTdb;
			cached = true;
		}
		return rotation;
	}
}  // namespace orbital_core::forces
