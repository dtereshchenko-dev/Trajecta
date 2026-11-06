#include "hazard/EncounterPlane.h"

#include "hazard/OpikTargetPlane.h"
#include "orbital_core/Constants.h"

#include <algorithm>
#include <cmath>

namespace orbital_core::hazard {
	namespace {
		constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
		constexpr double kAuPerDayToKmS = kAuKm / Constants::SECONDS_PER_DAY;

		[[nodiscard]] BPlaneBasis basisFromPlaneNormal(const Vec3D& normal) noexcept
		{
			BPlaneBasis basis;
			const double mag = normal.magnitude();
			if (!(mag > 1.0e-12) || !std::isfinite(mag))
				return basis;

			const Vec3D sHat = normal * (1.0 / mag);
			Vec3D reference{0.0, 0.0, 1.0};
			Vec3D xi = sHat.cross(reference);
			double xiMag = xi.magnitude();
			if (xiMag < 1.0e-12) {
				reference = Vec3D{1.0, 0.0, 0.0};
				xi = sHat.cross(reference);
				xiMag = xi.magnitude();
			}
			if (xiMag < 1.0e-12)
				return basis;

			basis.s = sHat;
			basis.xi = xi * (1.0 / xiMag);
			basis.zeta = sHat.cross(basis.xi);
			basis.valid = true;
			return basis;
		}

		[[nodiscard]] EncounterPlanePoint scaledBPlanePoint(const Vec3D& rGeoKm, const Vec3D& vGeoKmS,
		                                                    Planet planet) noexcept
		{
			EncounterPlanePoint out;
			out.kind = EncounterPlaneKind::ScaledBPlane;

			const ScaledBPlanePoint point = computeScaledBPlanePoint(rGeoKm, vGeoKmS, planet);
			if (!point.valid)
				return out;

			out.coordsKm = point.bScaledKm;
			out.basis = point.basis;
			out.gamma = point.gamma;
			out.vInfKmS = point.vInfKmS;
			out.valid = true;
			return out;
		}

		[[nodiscard]] EncounterPlanePoint opikTargetPlanePoint(const Vec3D& rGeoKm, const Vec3D& vGeoKmS,
		                                                       Planet planet) noexcept
		{
			EncounterPlanePoint out;
			out.kind = EncounterPlaneKind::OpikTargetPlane;

			const OpikTargetPlanePoint point = computeOpikTargetPlanePoint(rGeoKm, vGeoKmS, planet);
			if (!point.valid)
				return out;

			out.coordsKm = point.coordsKm;
			out.basis = point.basis;
			out.gamma = point.gamma;
			out.vInfKmS = point.vInfKmS;
			out.valid = true;
			return out;
		}

		[[nodiscard]] EncounterPlanePoint modifiedTargetPlanePoint(const Vec3D& rGeoKm, const Vec3D& vGeoKmS) noexcept
		{
			EncounterPlanePoint out;
			out.kind = EncounterPlaneKind::ModifiedTargetPlane;
			out.basis = basisFromPlaneNormal(vGeoKmS);
			if (!out.basis.valid)
				return out;

			out.coordsKm(0) = rGeoKm.dot(out.basis.xi);
			out.coordsKm(1) = rGeoKm.dot(out.basis.zeta);
			out.gamma = 1.0;
			out.vInfKmS = vGeoKmS.magnitude();
			out.valid = std::isfinite(out.coordsKm(0)) && std::isfinite(out.coordsKm(1)) && std::isfinite(out.vInfKmS);
			return out;
		}
	}  // namespace

	EncounterPlanePoint encounterPlaneAtHeliocentricState(const Vec3D& rTargetAu, const Vec3D& vTargetAuPerDay,
	                                                      const Vec3D& rPlanetAu, const Vec3D& vPlanetAuPerDay,
	                                                      Planet planet, EncounterPlaneKind kind) noexcept
	{
		const Vec3D rGeoKm = (rTargetAu - rPlanetAu) * kAuKm;
		const Vec3D vGeoKmS = (vTargetAuPerDay - vPlanetAuPerDay) * kAuPerDayToKmS;

		if (kind == EncounterPlaneKind::ModifiedTargetPlane)
			return modifiedTargetPlanePoint(rGeoKm, vGeoKmS);
		if (kind == EncounterPlaneKind::OpikTargetPlane)
			return opikTargetPlanePoint(rGeoKm, vGeoKmS, planet);

		return scaledBPlanePoint(rGeoKm, vGeoKmS, planet);
	}

	EncounterPlaneJacobian encounterPlaneJacobianAtHeliocentricState(const Vec3D& rTargetAu,
	                                                                 const Vec3D& vTargetAuPerDay,
	                                                                 const Vec3D& rPlanetAu,
	                                                                 const Vec3D& vPlanetAuPerDay, Planet planet,
	                                                                 EncounterPlaneKind kind) noexcept
	{
		EncounterPlaneJacobian out;
		const EncounterPlanePoint center =
		    encounterPlaneAtHeliocentricState(rTargetAu, vTargetAuPerDay, rPlanetAu, vPlanetAuPerDay, planet, kind);
		if (!center.valid)
			return out;

		out.coordsKm = center.coordsKm;
		out.basis = center.basis;
		out.gamma = center.gamma;
		out.kind = center.kind;

		if (kind == EncounterPlaneKind::OpikTargetPlane) {
			const Vec3D rGeoKm = (rTargetAu - rPlanetAu) * kAuKm;
			const Vec3D vGeoKmS = (vTargetAuPerDay - vPlanetAuPerDay) * kAuPerDayToKmS;
			const OpikTargetPlaneDifferential differential =
			    computeOpikTargetPlaneDifferential(rGeoKm, vGeoKmS, planet);
			if (!differential.valid)
				return out;

			out.dCoordsDTargetState.block<2, 3>(0, 0) = differential.dCoordsDGeoState.block<2, 3>(0, 0) * kAuKm;
			out.dCoordsDTargetState.block<2, 3>(0, 3) =
			    differential.dCoordsDGeoState.block<2, 3>(0, 3) * kAuPerDayToKmS;
			out.valid = true;
			return out;
		}

		Eigen::Matrix<double, 6, 1> state;
		state << rTargetAu.x, rTargetAu.y, rTargetAu.z, vTargetAuPerDay.x, vTargetAuPerDay.y, vTargetAuPerDay.z;

		for (int col = 0; col < 6; ++col) {
			const double base = state(col);
			const double step =
			    (col < 3) ? std::max(1.0e-10, std::abs(base) * 1.0e-8) : std::max(1.0e-12, std::abs(base) * 1.0e-8);

			auto plus = state;
			auto minus = state;
			plus(col) += step;
			minus(col) -= step;

			const EncounterPlanePoint p =
			    encounterPlaneAtHeliocentricState(Vec3D(plus(0), plus(1), plus(2)), Vec3D(plus(3), plus(4), plus(5)),
			                                      rPlanetAu, vPlanetAuPerDay, planet, kind);
			const EncounterPlanePoint m = encounterPlaneAtHeliocentricState(Vec3D(minus(0), minus(1), minus(2)),
			                                                                Vec3D(minus(3), minus(4), minus(5)),
			                                                                rPlanetAu, vPlanetAuPerDay, planet, kind);

			if (!p.valid || !m.valid)
				return out;

			out.dCoordsDTargetState.col(col) = (p.coordsKm - m.coordsKm) / (2.0 * step);
		}

		out.valid = true;
		return out;
	}

}  // namespace orbital_core::hazard
