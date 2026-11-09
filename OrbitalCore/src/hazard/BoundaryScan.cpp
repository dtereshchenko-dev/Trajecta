#include "BoundaryScan.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/SpkReader.h"
#include <cmath>
#include <stdexcept>

namespace orbital_core::hazard::boundary {
	std::optional<NominalEarthEntry> nearbyNominalEarthEntry(const FitResult& fit, const NBodySettings& forces,
	                                                        std::stop_token cancel)
	{
		if (!fit.converged || !fit.orbit || !forces.enabled)
			return std::nullopt;
		auto* eph = getDefaultEphemerisSource();
		if (!eph || !eph->isAvailable())
			throw std::runtime_error("Planetary ephemerides are unavailable for nominal entry search");
		const auto& orbit = *fit.orbit;
		const double epoch = orbit.epochMjd();
		const auto position = orbit.getR(), velocity = orbit.getV() * Constants::GAUSS_K;
		const auto earth = eph->queryMjd(EphemerisPlanet::Earth, epoch);
		const auto r = position - earth.R, v = velocity - earth.V;
		const double distance = r.magnitude(), sunDistance = position.magnitude();
		if (!std::isfinite(distance) || distance <= 0 || !std::isfinite(sunDistance) || sunDistance <= 0)
			throw std::runtime_error("Invalid state for nominal entry search");
		if (r.dot(v) >= 0 || Constants::GM_EARTH_ALONE / std::pow(distance, 2.5) <=
		                         Constants::GM_SUN / std::pow(sunDistance, 2.5))
			return std::nullopt;
		constexpr double altitudeKm = 50.0;
		const auto h = r.cross(v);
		const auto eccentricity = v.cross(h) / Constants::GM_EARTH_ALONE - r / distance;
		const double periapsis = h.dot(h) / (Constants::GM_EARTH_ALONE * (1 + eccentricity.magnitude()));
		const double boundaryAu = (planetRadiusKm(Planet::Earth) + altitudeKm) / (Constants::AU_IN_METERS / 1000.0);
		if (!std::isfinite(periapsis) || periapsis > boundaryAu)
			return std::nullopt;
		if (!std::isfinite(forces.stepSize) || forces.stepSize <= 0)
			throw std::invalid_argument("Invalid accepted propagation step");
		if (forces.useSb441 && (!getDefaultAsteroidSource() || !getDefaultAsteroidSource()->isAvailable()))
			throw std::runtime_error("SB441 is unavailable for nominal entry search");
		NBodyPropagator propagator(IntegratorKind::RA15);
		propagator.setStepSize(forces.stepSize);
		propagator.setPlanetMask(forces.planetMask);
		propagator.setUseSb441(forces.useSb441);
		propagator.setTargetNaifId(forces.targetNaifId);
		propagator.setInactiveInteriorMassPolicy(forces.inactiveInteriorMassPolicy);
		propagator = propagator.forEncounterMonitoring(*eph);
		const std::array<double, 6> state{position.x, position.y, position.z, velocity.x, velocity.y, velocity.z};
		const auto result = search(state, epoch, fit.finalNonGrav, Planet::Earth, epoch, epoch + 2.0, altitudeKm,
		                           propagator, *eph, cancel);
		if (!result.entry)
			return std::nullopt;
		return NominalEarthEntry{altitudeKm, *result.entry};
	}

	SearchResult search(const std::array<double, 6>& state, double epochMjdTdb, const NonGravModel& nonGrav,
	                    Planet planet, double startMjdTdb, double endMjdTdb, double altitudeKm,
	                    NBodyPropagator& propagator, IEphemerisSource& ephemeris, std::stop_token cancel)
	{
		constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
		const auto target = toEphemerisPlanet(planet);
		std::vector<double> targets;
		if (startMjdTdb > epochMjdTdb)
			targets.push_back(startMjdTdb);
		targets.push_back(endMjdTdb);
		std::vector<NBodyPropagator::StateVel> outputs;
		const auto& x = state;
		const auto terminal = propagator.propagateUntilSurfaceMjd(
		    {x[0], x[1], x[2]}, {x[3], x[4], x[5]}, epochMjdTdb, targets, outputs, target,
		    planetRadiusKm(planet), ephemeris, nonGrav,
		    planet == Planet::Earth ? NBodyPropagator::SurfaceShape::ReferenceEllipsoid
		                                   : NBodyPropagator::SurfaceShape::Sphere,
		    altitudeKm, cancel);
		if (cancel.stop_requested())
			throw std::runtime_error("MC cancelled");
		if (terminal.hitSurface && terminal.mjd == epochMjdTdb)
			throw std::runtime_error("MC clone starts inside the requested boundary");
		if (!terminal.hitSurface && terminal.mjd != endMjdTdb)
			throw std::runtime_error("Incomplete MC interval propagation");
		if (!terminal.hitSurface)
			return {};
		if (terminal.mjd < startMjdTdb)
			return {.beforeInterval = true};
		const auto body = ephemeris.queryMjd(target, terminal.mjd);
		const auto r = (terminal.R - body.R) * auKm;
		const auto v = (terminal.V - body.V) * (auKm / 86400.0);
		if (!std::isfinite(r.dot(r)) || !std::isfinite(v.dot(v)))
			throw std::runtime_error("Non-finite MC boundary state");
		const auto fixed = icrfToBodyFixedMjd(planet, r, terminal.mjd, &ephemeris);
		const auto ll = latLonFromBodyFixed(planet, fixed);
		if (!ll.valid)
			throw std::runtime_error("MC body-fixed conversion failed");
		return {.entry = BoundaryCrossing{terminal.mjd, r, v, ll.latitudeDeg, ll.longitudeDeg}};
	}
}  // namespace orbital_core::hazard::boundary
