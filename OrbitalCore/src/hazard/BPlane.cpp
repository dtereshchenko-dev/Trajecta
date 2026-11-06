#include "orbital_core/hazard/BPlane.h"
#include <span>

#include "orbital_core/Constants.h"
#include "orbital_core/hazard/PlanetData.h"

#include <algorithm>
#include <cmath>

namespace orbital_core::hazard {
	namespace {
		constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;
		constexpr double kAuPerDayToKmS = kAuKm / 86400.0;
		const Vec3D kEclipticNormal{0.0, 0.0, 1.0};

		BPlaneBasis computeBPlaneBasisFromAsymptote(const Vec3D& sHat) noexcept
		{
			BPlaneBasis basis;
			if (std::abs(sHat.magnitude() - 1.0) > 1.0e-6)
				return basis;

			Vec3D reference = kEclipticNormal;
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
	}  // namespace

	BPlaneBasis computeBPlaneBasisFromVrel(const Vec3D& vRelAuPerDay) noexcept
	{
		BPlaneBasis basis;
		const double vRelKmS = vRelAuPerDay.magnitude() * kAuPerDayToKmS;
		if (vRelKmS < 1.0e-9)
			return basis;

		const Vec3D sHat = vRelAuPerDay.normalized();
		const Vec3D sCrossZ = sHat.cross(kEclipticNormal);
		const double sCrossZmag = sCrossZ.magnitude();
		if (sCrossZmag < 1.0e-12)
			return basis;

		basis.s = sHat;
		basis.xi = sCrossZ * (1.0 / sCrossZmag);
		basis.zeta = sHat.cross(basis.xi);
		basis.valid = true;
		return basis;
	}

	BPlaneBasis computeBPlaneBasis(const CloseApproachEvent& event) noexcept
	{
		return computeBPlaneBasisFromVrel(event.vTarget - event.vPlanet);
	}

	BPlanePoint mapToBPlane(const Vec3D& missVecKm, const BPlaneBasis& basis) noexcept
	{
		BPlanePoint p;
		if (!basis.valid)
			return p;
		p.xi = missVecKm.dot(basis.xi);
		p.zeta = missVecKm.dot(basis.zeta);
		return p;
	}

	BPlaneCovariance fitBPlaneCovariance(std::span<const Vec3D> samplesKm, const BPlaneBasis& basis) noexcept
	{
		BPlaneCovariance cov;
		cov.numSamples = static_cast<int>(samplesKm.size());
		if (cov.numSamples < 2 || !basis.valid)
			return cov;

		double sumXi = 0.0, sumZeta = 0.0;
		for (const auto& m : samplesKm) {
			const BPlanePoint p = mapToBPlane(m, basis);
			sumXi += p.xi;
			sumZeta += p.zeta;
		}
		cov.meanXi = sumXi / cov.numSamples;
		cov.meanZeta = sumZeta / cov.numSamples;

		double varXi = 0.0, varZeta = 0.0, covXiZeta = 0.0;
		for (const auto& m : samplesKm) {
			const BPlanePoint p = mapToBPlane(m, basis);
			const double dXi = p.xi - cov.meanXi;
			const double dZeta = p.zeta - cov.meanZeta;
			varXi += dXi * dXi;
			varZeta += dZeta * dZeta;
			covXiZeta += dXi * dZeta;
		}
		const double denom = static_cast<double>(cov.numSamples - 1);
		varXi /= denom;
		varZeta /= denom;
		covXiZeta /= denom;

		cov.sigmaXi = std::sqrt(std::max(varXi, 0.0));
		cov.sigmaZeta = std::sqrt(std::max(varZeta, 0.0));

		const double sigmaProd = cov.sigmaXi * cov.sigmaZeta;
		cov.correlation = (sigmaProd > 0.0) ? std::clamp(covXiZeta / sigmaProd, -1.0, 1.0) : 0.0;
		return cov;
	}

	double gravitationalFocusing(double vInfKmS, Planet planet) noexcept
	{
		if (!(vInfKmS > 1.0e-6))
			return 1.0;
		const auto& pd = getPlanetData(planet);
		const double v2 = vInfKmS * vInfKmS;
		return std::sqrt(1.0 + 2.0 * pd.muKm3PerS2 / (pd.radiusEqKm * v2));
	}

	ScaledBPlanePoint computeScaledBPlanePoint(const Vec3D& rGeoKm, const Vec3D& vGeoKmS, Planet planet) noexcept
	{
		ScaledBPlanePoint out;
		const auto& pd = getPlanetData(planet);
		const double mu = pd.muKm3PerS2;
		if (!(mu > 0.0))
			return out;

		const double r = rGeoKm.magnitude();
		const double v2 = vGeoKmS.dot(vGeoKmS);
		if (!(r > 0.0) || !(v2 > 0.0))
			return out;

		const double energy = 0.5 * v2 - mu / r;
		if (!(energy > 0.0))
			return out;

		const double vInf = std::sqrt(2.0 * energy);
		const Vec3D hVec = rGeoKm.cross(vGeoKmS);
		const double h = hVec.magnitude();
		if (!(h > 0.0))
			return out;

		const Vec3D eVec = vGeoKmS.cross(hVec) * (1.0 / mu) - rGeoKm * (1.0 / r);
		const double e = eVec.magnitude();
		if (!(e > 1.0))
			return out;

		const Vec3D pHat = eVec * (1.0 / e);
		const Vec3D hHat = hVec * (1.0 / h);
		const Vec3D qHat = hHat.cross(pHat);
		const double asymptoteSlope = std::sqrt(std::max(e * e - 1.0, 0.0));
		const Vec3D incomingVInfHat = (pHat + qHat * asymptoteSlope).normalized();

		out.basis = computeBPlaneBasisFromAsymptote(incomingVInfHat);
		if (!out.basis.valid)
			return out;

		const Vec3D bVecKm = incomingVInfHat.cross(hVec) * (1.0 / vInf);
		out.gamma = gravitationalFocusing(vInf, planet);
		out.vInfKmS = vInf;
		out.bScaledKm(0) = bVecKm.dot(out.basis.xi) / out.gamma;
		out.bScaledKm(1) = bVecKm.dot(out.basis.zeta) / out.gamma;
		out.valid = std::isfinite(out.bScaledKm(0)) && std::isfinite(out.bScaledKm(1)) && std::isfinite(out.gamma) &&
		            out.gamma > 0.0;
		return out;
	}
}  // namespace orbital_core::hazard
