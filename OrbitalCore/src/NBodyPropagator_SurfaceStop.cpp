
#include "NBodyPropagator.h"
#include <span>
#include <vector>
#include "orbital_core/TimeScale.h"

#include "orbital_core/Constants.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include "integrators/IIntegrator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace {

	constexpr double kAuKm = Constants::AU_IN_METERS / 1000.0;

	[[nodiscard]]
	std::optional<orbital_core::hazard::Planet> hazardPlanetFromEphemeris(EphemerisPlanet planet) noexcept
	{
		using orbital_core::hazard::Planet;

		switch (planet) {
			case EphemerisPlanet::Mercury:
				return Planet::Mercury;
			case EphemerisPlanet::Venus:
				return Planet::Venus;
			case EphemerisPlanet::Earth:
				return Planet::Earth;
			case EphemerisPlanet::Mars:
				return Planet::Mars;
			case EphemerisPlanet::Jupiter:
				return Planet::Jupiter;
			case EphemerisPlanet::Saturn:
				return Planet::Saturn;
			case EphemerisPlanet::Uranus:
				return Planet::Uranus;
			case EphemerisPlanet::Neptune:
				return Planet::Neptune;
			case EphemerisPlanet::Moon:
				return Planet::Moon;
			case EphemerisPlanet::Sun:
			case EphemerisPlanet::Pluto:
				break;
		}
		return std::nullopt;
	}

	[[nodiscard]]
	double planetSurfaceOffsetKm(const Vec3D& rCloneAu, double mjd, EphemerisPlanet planet, double sphericalRadiusKm,
	                             const EphemerisRecord& rec, NBodyPropagator::SurfaceShape shape, double altitudeKm)
	{
		const Vec3D rPlanetCentricEclKm = (rCloneAu - rec.R) * kAuKm;

		if (shape == NBodyPropagator::SurfaceShape::ReferenceEllipsoid) {
			const auto hazardPlanet = hazardPlanetFromEphemeris(planet);
			if (!hazardPlanet)
				return rPlanetCentricEclKm.magnitude() - sphericalRadiusKm - altitudeKm;
			using namespace orbital_core::hazard;
			const Vec3D bodyFixedKm = icrfToBodyFixedMjd(*hazardPlanet, rPlanetCentricEclKm, mjd);
			if (!std::isfinite(bodyFixedKm.dot(bodyFixedKm)))
				throw std::runtime_error("Body orientation unavailable at surface-crossing epoch");
			return referenceSurfaceSignedAltitudeKm(bodyFixedKm, planetRadiusKm(*hazardPlanet),
			                                        planetFlattening(*hazardPlanet)) -
			       altitudeKm;
		}

		return rPlanetCentricEclKm.magnitude() - sphericalRadiusKm - altitudeKm;
	}

	struct CrossingState {
		Vec3D R;
		Vec3D V;
		double mjd;
	};

	constexpr double kCrossingToleranceDays = 1e-8;

	struct RadialSample {
		double derivative;
		double offsetKm;
		double speedKmPerDay;
	};

	template <class Sample>
	bool radialMinimumClearsSurface(double dt, double derivativeBefore, const Sample& sample)
	{
		const auto end = sample(1.0);
		if (!(derivativeBefore < 0 && end.derivative > 0 && end.offsetKm > 0))
			return false;

		double lo = 0, hi = 1;
		double flo = derivativeBefore, fhi = end.derivative;
		int previousSide = 0;
		for (int i = 0; i < 64; ++i) {
			double u = lo + (hi - lo) * (-flo / (fhi - flo));
			if (!(lo < u && u < hi))
				u = (lo + hi) * .5;
			const auto value = sample(u);
			if (!std::isfinite(value.derivative) || !(value.offsetKm > 0))
				return false;
			if (value.derivative == 0 || std::abs(dt) * (hi - lo) <= kCrossingToleranceDays) {
				return value.offsetKm > value.speedKmPerDay * kCrossingToleranceDays;
			}
			const int side = value.derivative < 0 ? -1 : 1;
			if (side < 0) {
				lo = u;
				flo = value.derivative;
				if (previousSide == side)
					fhi *= .5;
			} else {
				hi = u;
				fhi = value.derivative;
				if (previousSide == side)
					flo *= .5;
			}
			previousSide = side;
		}
		return false;
	}

	template <class Evaluate, class Offset>
	std::optional<CrossingState> findStepCrossing(const AcceptedStateStep& step, double offsetAfter,
	                                              const Evaluate& evaluate, const Offset& offset, double toleranceDays)
	{
		const double kTolMjd = toleranceDays;
		const double dt = step.endMjd - step.beginMjd;
		const auto at = [&](double u) { return evaluate(step.beginMjd + dt * u); };
		double lo = 0, hi = 1;
		if (offsetAfter > 0) {
			constexpr double ratio = 0.6180339887498948482;
			double a = 0, b = 1;
			double u = b - ratio * (b - a), v = a + ratio * (b - a);
			double fu = offset(at(u)), fv = offset(at(v));
			for (int i = 0; i < 64 && std::abs(dt) * (b - a) > kTolMjd; ++i) {
				if (fu < fv) {
					b = v;
					v = u;
					fv = fu;
					u = b - ratio * (b - a);
					fu = offset(at(u));
				} else {
					a = u;
					u = v;
					fu = fv;
					v = a + ratio * (b - a);
					fv = offset(at(v));
				}
			}
			hi = fu < fv ? u : v;
			if (std::min(fu, fv) > 0)
				return std::nullopt;
		}
		for (int i = 0; i < 64 && std::abs(dt) * (hi - lo) > kTolMjd; ++i) {
			const double mid = (lo + hi) * .5;
			if (offset(at(mid)) > 0)
				lo = mid;
			else
				hi = mid;
		}
		return at(hi);
	}

}  // namespace

NBodyPropagator::PropagationWithSurfaceStop NBodyPropagator::propagateUntilSurface(
    const Vec3D& R0, const Vec3D& V0, double epoch, double maxTargetJd, EphemerisPlanet planet, double rPlanetKm,
    IEphemerisSource& eph, const NonGravModel& ng, SurfaceShape shape, double altitudeKm) const
{
	const std::array<double, 1> targets{maxTargetJd};
	std::vector<StateVel> samples;
	return propagateUntilSurface(R0, V0, epoch, targets, samples, planet, rPlanetKm, eph, ng, shape, altitudeKm);
}

NBodyPropagator::PropagationWithSurfaceStop NBodyPropagator::propagateUntilSurface(
    const Vec3D& R0, const Vec3D& V0, double epoch, std::span<const double> targetJds, std::vector<StateVel>& samples,
    EphemerisPlanet planet, double rPlanetKm, IEphemerisSource& eph, const NonGravModel& ng, SurfaceShape shape,
    double altitudeKm) const
{
	std::vector<double> targets;
	targets.reserve(targetJds.size());
	for (const double jd : targetJds)
		targets.push_back(orbital_core::time::mjdFromJd(jd));
	return propagateUntilSurfaceMjd(R0, V0, orbital_core::time::mjdFromJd(epoch), targets, samples, planet, rPlanetKm,
	                                eph, ng, shape, altitudeKm);
}

NBodyPropagator::PropagationWithSurfaceStop NBodyPropagator::propagateUntilSurfaceMjd(
    const Vec3D& R0, const Vec3D& V0, double epochMjd, std::span<const double> targetMjds,
    std::vector<StateVel>& samples, EphemerisPlanet planet, double rPlanetKm, IEphemerisSource& eph,
    const NonGravModel& ng, SurfaceShape shape, double altitudeKm, std::stop_token cancel) const
{
	if (!std::isfinite(altitudeKm) || altitudeKm < 0.0)
		throw std::invalid_argument("Surface altitude must be finite and nonnegative");
	const double toleranceDays = altitudeKm > 0.0 ? 1e-11 : kCrossingToleranceDays;
	samples.clear();
	PropagationWithSurfaceStop out;
	out.R = R0;
	out.V = V0;
	out.jd = orbital_core::time::jdFromMjd(epochMjd);
	out.mjd = epochMjd;
	out.hitSurface = false;

	if (targetMjds.empty() || cancel.stop_requested())
		return out;

	double boundingRadiusKm = rPlanetKm;
	if (shape == SurfaceShape::ReferenceEllipsoid) {
		if (const auto body = hazardPlanetFromEphemeris(planet)) {
			const auto& data = orbital_core::hazard::getPlanetData(*body);
			boundingRadiusKm = std::max(data.radiusEqKm, data.radiusEqKm * (1.0 - data.flattening));
		}
	}

	EphemerisRecord cachedPlanet;
	double cachedMjd = epochMjd;
	cachedPlanet = eph.queryMjd(planet, epochMjd);
	const double offset0 = planetSurfaceOffsetKm(R0, epochMjd, planet, rPlanetKm, cachedPlanet, shape, altitudeKm);
	if (offset0 <= 0.0) {
		out.hitSurface = true;
		return out;
	}

	std::optional<CrossingState> crossing;
	const AcceptedStateStepObserver observer = [&](const AcceptedStateStep& step) {
		EphemerisRecord p0, p1;
		p0 = cachedMjd == step.beginMjd ? cachedPlanet : eph.queryMjd(planet, step.beginMjd);
		p1 = eph.queryMjd(planet, step.endMjd);
		cachedPlanet = p1;
		cachedMjd = step.endMjd;
		const Vec3D r0 = step.begin.R - p0.R, r1 = step.end.R - p1.R;
		const double radiusBeforeKm = r0.magnitude() * kAuKm;
		const double radiusAfterKm = r1.magnitude() * kAuKm;
		const double boundingOffsetBefore = radiusBeforeKm - boundingRadiusKm - altitudeKm;
		const double boundingOffsetAfter = radiusAfterKm - boundingRadiusKm - altitudeKm;
		const double dt = step.endMjd - step.beginMjd;
		const double derivativeBefore = r0.dot(step.begin.V - p0.V) * dt;
		const bool radialMinimum = derivativeBefore < 0 && r1.dot(step.end.V - p1.V) * dt > 0;
		if (!radialMinimum && std::min(radiusBeforeKm, radiusAfterKm) > 2 * (boundingRadiusKm + altitudeKm))
			return true;

		const double offsetBefore =
		    boundingOffsetBefore > 0
		        ? boundingOffsetBefore
		        : planetSurfaceOffsetKm(step.begin.R, step.beginMjd, planet, rPlanetKm, p0, shape, altitudeKm);
		const double offsetAfter = boundingOffsetAfter > 0 ? boundingOffsetAfter
		                                                   : planetSurfaceOffsetKm(step.end.R, step.endMjd, planet,
		                                                                           rPlanetKm, p1, shape, altitudeKm);
		if (offsetBefore <= 0)
			return true;
		const auto evaluate = [&](double mjd) {
			if (step.interpolateMjd) {
				const auto state = step.interpolateMjd(mjd);
				return CrossingState{state.R, state.V, mjd};
			}
			const auto state =
			    propagateAllWithVelocityMjd(step.begin.R, step.begin.V, step.beginMjd, {mjd}, ng, cancel).at(0);
			return CrossingState{state.R, state.V, mjd};
		};
		const auto offset = [&](const CrossingState& state) {
			return planetSurfaceOffsetKm(state.R, state.mjd, planet, rPlanetKm, eph.queryMjd(planet, state.mjd), shape,
			                             altitudeKm);
		};
		if (boundingOffsetBefore > 0 && boundingOffsetAfter > 0 && radialMinimum) {
			const auto sample = [&](double u) {
				const auto state = evaluate(step.beginMjd + dt * u);
				const auto p = eph.queryMjd(planet, state.mjd);
				const Vec3D r = state.R - p.R, v = state.V - p.V;
				return RadialSample{r.dot(v) * dt, r.magnitude() * kAuKm - boundingRadiusKm - altitudeKm,
				                    v.magnitude() * kAuKm};
			};
			if (radialMinimumClearsSurface(dt, derivativeBefore, sample))
				return true;
		}
		crossing = findStepCrossing(step, offsetAfter, evaluate, offset, toleranceDays);
		return !crossing.has_value();
	};

	auto options = integrationOptions();
	options.cancel = cancel;
	std::vector<IntegratedState> states(targetMjds.size());
	const IntegrationReport report =
	    m_integrator->integrateState(m_dynamics, ng, IntegratedState{R0, V0}, epochMjd, targetMjds,
	                                 std::span<IntegratedState>(states), options, observer);

	samples.reserve(report.completedOutputs);
	for (std::size_t i = 0; i < report.completedOutputs; ++i) {
		if (crossing && (targetMjds[i] - crossing->mjd) * (targetMjds.back() - epochMjd) > 0.0)
			break;
		samples.push_back({states[i].R, states[i].V});
	}

	if (crossing) {
		out.R = crossing->R;
		out.V = crossing->V;
		out.jd = orbital_core::time::jdFromMjd(crossing->mjd);
		out.mjd = crossing->mjd;
		out.hitSurface = true;
		return out;
	}

	out.R = report.finalState.R;
	out.V = report.finalState.V;
	out.jd = orbital_core::time::jdFromMjd(report.reachedMjd);
	out.mjd = report.reachedMjd;
	return out;
}

NBodyPropagator::PropagationWithSurfaceStop NBodyPropagator::propagateUntilSurface(
    const Vec3D& R0, const Vec3D& V0, double epoch, double maxTargetJd, EphemerisPlanet planet, double rPlanetKm,
    IEphemerisSource& eph, SurfaceShape shape, double altitudeKm) const
{
	return propagateUntilSurface(R0, V0, epoch, maxTargetJd, planet, rPlanetKm, eph, m_nonGrav, shape, altitudeKm);
}
