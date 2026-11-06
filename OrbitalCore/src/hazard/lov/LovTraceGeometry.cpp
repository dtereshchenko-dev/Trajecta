#include "hazard/lov/LovTraceDetail.h"
#include <span>
#include <vector>

#include "orbital_core/Constants.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>

namespace orbital_core::hazard::lov_detail {

	const LovCatalogPoint* catalogAnchorAtSigma(std::span<const LovCatalogPoint> catalog, double sigma) noexcept
	{
		if (!std::isfinite(sigma))
			return nullptr;

		const LovCatalogPoint* nearest = nullptr;
		double nearestDistance = std::numeric_limits<double>::infinity();
		for (const LovCatalogPoint& point : catalog) {
			if (!point.valid || !std::isfinite(point.sigma))
				continue;
			const double distance = std::abs(point.sigma - sigma);
			if (distance < nearestDistance ||
			    (nearest && distance == nearestDistance && point.sigma > nearest->sigma)) {
				nearest = &point;
				nearestDistance = distance;
			}
		}
		return nearest;
	}

	double angleDifferenceAbs(double a, double b) noexcept { return std::abs(std::remainder(a - b, Constants::TWOPI)); }

	static double safeAcos(double x) noexcept { return std::acos(std::clamp(x, -1.0, 1.0)); }

	static double velocityBodyRadiiPerDay(double vInfKmS, double bodyRadiusKm) noexcept
	{
		if (!(vInfKmS > 0.0) || !(bodyRadiusKm > 0.0))
			return 0.0;
		return vInfKmS * Constants::SECONDS_PER_DAY / bodyRadiusKm;
	}

	LovPairType classifyLovPair(const TraceRecord& left, const TraceRecord& right) noexcept
	{
		const TargetPlaneTracePoint& a = left.trace;
		const TargetPlaneTracePoint& b = right.trace;
		if (!a.valid || !b.valid || !(a.bodyRadiusKm > 0.0))
			return LovPairType::Unknown;

		const Eigen::Vector2d delta = (b.coordsKm - a.coordsKm) / a.bodyRadiusKm;
		const double chord = delta.norm();
		if (!(chord > 0.0))
			return LovPairType::NoConcern;

		const auto tangent = [](const TargetPlaneTracePoint& t) {
			return Eigen::Vector2d(-t.stretchLovBodyRadiiPerSigma * std::sin(t.alphaLovRad),
			                       t.stretchLovBodyRadiiPerSigma * std::cos(t.alphaLovRad));
		};

		const Eigen::Vector2d ta = tangent(a);
		const Eigen::Vector2d tb = tangent(b);
		const double taNorm = ta.norm();
		const double tbNorm = tb.norm();
		if (!(taNorm > 0.0) || !(tbNorm > 0.0))
			return LovPairType::Unknown;

		const double dalpha = angleDifferenceAbs(a.alphaLovRad, b.alphaLovRad);
		const double dbeta1 = safeAcos(delta.dot(ta) / (chord * taNorm));
		const double dbeta2 = safeAcos(delta.dot(tb) / (chord * tbNorm));

		constexpr double kBetaFactor = 20.0;
		constexpr double kSin20 = 0.34202014332566873304;
		const double ua = velocityBodyRadiiPerDay(left.event.vInfKmS, a.bodyRadiusKm);
		const double ub = velocityBodyRadiiPerDay(right.event.vInfKmS, b.bodyRadiusKm);
		const double betaLimit1 =
		    ua > 0.0 ? std::max(kBetaFactor * Constants::GAUSS_K * a.distanceBodyRadii / ua, kSin20) : kSin20;
		const double betaLimit2 =
		    ub > 0.0 ? std::max(kBetaFactor * Constants::GAUSS_K * b.distanceBodyRadii / ub, kSin20) : kSin20;

		const double cosDa = std::cos(dalpha);
		const bool derivativeMinimum = a.dDistance2DSigmaBodyRadii2 < 0.0 && b.dDistance2DSigmaBodyRadii2 > 0.0;
		const bool derivativeMaximum = a.dDistance2DSigmaBodyRadii2 > 0.0 && b.dDistance2DSigmaBodyRadii2 < 0.0;

		if (cosDa > 0.0) {
			if (std::abs(std::sin(dbeta1)) < betaLimit1 || std::abs(std::sin(dbeta2)) < betaLimit2) {
				if (derivativeMinimum)
					return LovPairType::SimpleMinimum;
				if (derivativeMaximum)
					return LovPairType::SimpleMaximum;
				return LovPairType::NoConcern;
			}

			return b.dDistance2DSigmaBodyRadii2 > 0.0 ? LovPairType::EntangledAscending
			                                          : LovPairType::EntangledDescending;
		}

		if (cosDa < 0.0) {
			if (derivativeMinimum)
				return LovPairType::InterruptedMinimum;
			if (derivativeMaximum)
				return LovPairType::InterruptedMaximum;
			return LovPairType::InterruptedUnknown;
		}

		return LovPairType::Unknown;
	}

	double closeApproachWindowDays(double encounterRadiusAu, double vInfKmS) noexcept
	{
		double windowDays = 365.25;
		const double vAuDay = vInfKmS * Constants::SECONDS_PER_DAY / kAuKm;
		if (vAuDay > 0.0 && std::isfinite(vAuDay))
			windowDays = 2.0 * encounterRadiusAu / vAuDay;

		return std::clamp(windowDays, 20.0, 365.25);
	}

	double endpointSigmaStepLimit(const std::vector<LovCatalogPoint>& catalog, int sampleIndex, bool increasingHead,
	                              double fallbackStep) noexcept
	{
		double deltaSigma = fallbackStep;
		const auto current = std::find_if(catalog.begin(), catalog.end(), [&](const LovCatalogPoint& point) {
			return point.sampleIndex == sampleIndex;
		});

		if (current != catalog.end()) {
			if (increasingHead && std::next(current) != catalog.end()) {
				deltaSigma = std::next(current)->sigma - current->sigma;
			} else if (!increasingHead && current != catalog.begin()) {
				deltaSigma = current->sigma - std::prev(current)->sigma;
			}
		}

		if (!(deltaSigma > 0.0) || !std::isfinite(deltaSigma))
			deltaSigma = fallbackStep;
		return 0.5 * deltaSigma;
	}

	double minimumPossibleBodyRadii(const TraceRecord& record, double sigmaMax) noexcept
	{
		const TargetPlaneTracePoint& trace = record.trace;
		double localMoid = trace.distanceBodyRadii;
		if (trace.bodyRadiusKm > 0.0) {
			if (std::isfinite(record.event.localMoidAu) && record.event.localMoidAu >= 0.0) {
				localMoid = record.event.localMoidAu * kAuKm / trace.bodyRadiusKm;
			} else if (std::isfinite(record.event.distAu) && record.event.distAu >= 0.0) {
				localMoid = record.event.distAu * kAuKm / trace.bodyRadiusKm;
			}
		}

		double dMin = localMoid - sigmaMax * trace.widthBodyRadii;
		if (dMin < 0.0)
			dMin = 0.0;
		dMin = std::min(dMin, trace.distanceBodyRadii);

		double along = trace.distanceBodyRadii;
		if (trace.dDistance2DSigmaBodyRadii2 > 0.0) {
			along -= (trace.sigma + sigmaMax) * trace.stretchBodyRadii;
		} else if (trace.dDistance2DSigmaBodyRadii2 < 0.0) {
			along -= (sigmaMax - trace.sigma) * trace.stretchBodyRadii;
		}

		return std::max(dMin, along);
	}

	double localMoidBodyRadii(const TraceRecord& record) noexcept
	{
		const TargetPlaneTracePoint& trace = record.trace;
		if (!(trace.bodyRadiusKm > 0.0))
			return std::numeric_limits<double>::quiet_NaN();

		if (std::isfinite(record.event.localMoidAu) && record.event.localMoidAu >= 0.0) {
			return record.event.localMoidAu * kAuKm / trace.bodyRadiusKm;
		}

		if (std::isfinite(record.event.distAu) && record.event.distAu >= 0.0) {
			return record.event.distAu * kAuKm / trace.bodyRadiusKm;
		}

		return std::numeric_limits<double>::quiet_NaN();
	}

}  // namespace orbital_core::hazard::lov_detail
