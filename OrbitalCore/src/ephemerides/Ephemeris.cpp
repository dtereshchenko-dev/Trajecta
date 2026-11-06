#include "orbital_core/Ephemeris.h"

#include "EphemeridesContainer.h"
#include "TwoBodyPropagation.h"
#include "NBodyPropagator.h"
#include "obs/LightTime.h"
#include "obs/ObservatoryCatalog.h"
#include "obs/ObserverGeometry.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core {
	namespace {

		struct ObserverSite {
			bool ground{false};
			obs::ObserverData data;
		};

		bool finite(const Vec3D& v) noexcept { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

		ObserverSite resolveObserver(const EphemerisRequest& request)
		{
			ObserverSite site;
			if (request.center != EphemerisCenter::Station)
				return site;
			const auto& catalog = obs::observatoryCatalog();
			const auto station = catalog.find(request.station);
			if (station == catalog.end())
				throw std::invalid_argument("Unknown observatory code: " + request.station);
			const auto& record = station->second;
			if (record.code != "500" && record.cosPhi == 0.0 && record.sinPhi == 0.0)
				throw std::invalid_argument("Observatory " + request.station +
				                            " has no fixed geodetic position (spacecraft)");
			site.ground = true;
			site.data = obs::ObserverData{record.longitude, record.cosPhi, record.sinPhi, record.code};
			return site;
		}

		NBodyPropagator ephemerisPropagator(const EphemerisSolution& solution)
		{
			NBodyPropagator propagator;
			propagator.setStepSize(solution.forces.stepSize);
			propagator.setPlanetMask(solution.forces.enabled ? solution.forces.planetMask : 0);
			propagator.setUseSb441(solution.forces.useSb441);
			propagator.setTargetNaifId(solution.forces.targetNaifId);
			propagator.setInactiveInteriorMassPolicy(solution.forces.inactiveInteriorMassPolicy);
			if (auto* ephemeris = getDefaultEphemerisSource())
				propagator = propagator.forEncounterMonitoring(*ephemeris);
			propagator.nonGrav() = solution.nonGrav;
			return propagator;
		}

		std::vector<NBodyPropagator::StateVel> receptionStates(const EphemerisSolution& solution,
		                                                       const std::vector<double>& mjdTdb, bool nBody,
		                                                       const NBodyPropagator& propagator,
		                                                       std::stop_token cancel)
		{
			if (nBody) {
				return propagator.propagateAllWithVelocityMjd(solution.position, solution.velocity,
				                                              solution.epochMjdTdb, mjdTdb, solution.nonGrav, cancel);
			}

			std::vector<NBodyPropagator::StateVel> states(mjdTdb.size());
			const StateVectors initial(solution.position, solution.velocity / Constants::GAUSS_K);
			for (std::size_t i = 0; i < mjdTdb.size(); ++i) {
				if (cancel.stop_requested())
					return {};
				const StateVectors state =
				    Calculations::propagateKeplerState(initial, mjdTdb[i] - solution.epochMjdTdb);
				if (!state.isValid())
					throw std::runtime_error("Two-body ephemeris propagation failed");
				states[i] = {state.first, state.second * Constants::GAUSS_K};
			}
			return states;
		}

		Vec3D toPlane(const Vec3D& ecliptic, EphemerisPlane plane, double cosObliquity, double sinObliquity) noexcept
		{
			if (plane == EphemerisPlane::Ecliptic)
				return ecliptic;
			return Vec3D(ecliptic.x, ecliptic.y * cosObliquity - ecliptic.z * sinObliquity,
			             ecliptic.y * sinObliquity + ecliptic.z * cosObliquity);
		}

		double angleBetween(const Vec3D& a, const Vec3D& b) noexcept
		{
			const double norms = a.magnitude() * b.magnitude();
			if (!(norms > 0.0))
				return std::numeric_limits<double>::quiet_NaN();
			return std::acos(std::clamp(a.dot(b) / norms, -1.0, 1.0));
		}

	}  // namespace

	std::size_t ephemerisEpochCount(double startMjd, double endMjd, double stepDays)
	{
		if (!std::isfinite(startMjd) || !std::isfinite(endMjd) || !std::isfinite(stepDays))
			throw std::invalid_argument("Ephemeris epochs must be finite");
		if (!(stepDays > 0.0))
			throw std::invalid_argument("Ephemeris step must be positive");
		if (endMjd < startMjd)
			throw std::invalid_argument("Ephemeris end precedes its start");
		if (startMjd == endMjd)
			return 1;
		const double largestEpoch = std::max(std::abs(startMjd), std::abs(endMjd));
		const double resolution = std::nextafter(largestEpoch, std::numeric_limits<double>::infinity()) - largestEpoch;
		if (stepDays < resolution)
			throw std::invalid_argument("Ephemeris step is below epoch resolution");
		std::size_t first = 1, end = kMaxEphemerisEpochs + 1;
		while (first < end) {
			const std::size_t middle = first + (end - first) / 2;
			if (startMjd + static_cast<double>(middle) * stepDays <= endMjd)
				first = middle + 1;
			else
				end = middle;
		}
		return first;
	}

	EphemerisResult computeEphemeris(const EphemerisRequest& request, const FitRunControl& control)
	{
		const auto& solution = request.solution;
		if (!finite(solution.position) || !finite(solution.velocity) || !std::isfinite(solution.epochMjdTdb))
			throw std::invalid_argument("Ephemeris requires a finite solution state and epoch");
		const std::size_t count = ephemerisEpochCount(request.startMjd, request.endMjd, request.stepDays);
		if (count > kMaxEphemerisEpochs)
			throw std::invalid_argument("Too many ephemeris epochs");

		const ObserverSite site = resolveObserver(request);
		const bool nBody = solution.forces.enabled || solution.nonGrav.enabled();
		if (nBody) {
			auto* ephemeris = getDefaultEphemerisSource();
			if (!ephemeris || !ephemeris->isAvailable())
				throw std::runtime_error("Planetary ephemerides are unavailable");
			if (solution.forces.useSb441 && (!getDefaultAsteroidSource() || !getDefaultAsteroidSource()->isAvailable()))
				throw std::runtime_error("SB441 was requested but is unavailable");
		}

		EphemerisResult result;
		result.nBody = nBody;
		result.nonGravitational = solution.nonGrav.enabled();
		result.lightTimeApplied = request.lightTime && request.center != EphemerisCenter::Heliocentric;
		const auto cancelled = [&] {
			result.cancelled = true;
			return result;
		};
		if (control.cancel.stop_requested())
			return cancelled();

		std::vector<EphemerisRow> rows(count);
		std::vector<double> mjdTdb(count);
		for (std::size_t i = 0; i < count; ++i) {
			if (control.cancel.stop_requested())
				return cancelled();
			auto& row = rows[i];
			row.mjd = request.startMjd + static_cast<double>(i) * request.stepDays;
			row.mjdTdb = time::toMjdTdb(row.mjd, request.scale);
			row.mjdTt = request.scale == time::Scale::TT ? row.mjd : time::fromMjdTdb(row.mjdTdb, time::Scale::TT);
			row.mjdUtc = request.scale == time::Scale::UTC ? row.mjd : time::fromMjdTdb(row.mjdTdb, time::Scale::UTC);
			mjdTdb[i] = row.mjdTdb;
		}

		if (control.cancel.stop_requested())
			return cancelled();
		if (control.progress)
			control.progress("Propagating ephemeris");
		const auto propagator = ephemerisPropagator(solution);
		const auto states = receptionStates(solution, mjdTdb, nBody, propagator, control.cancel);
		if (control.cancel.stop_requested())
			return cancelled();
		if (states.size() != count)
			throw std::runtime_error("Incomplete ephemeris propagation");

		EphemeridesContainer earth(EphemeridesContainer::Planet::Earth);
		const double cosObliquity = std::cos(Constants::EARTH_INCL);
		const double sinObliquity = std::sin(Constants::EARTH_INCL);
		for (std::size_t i = 0; i < count; ++i) {
			if (control.cancel.stop_requested())
				return cancelled();
			auto& row = rows[i];
			Vec3D target = states[i].R;
			Vec3D targetVelocity = states[i].V;
			Vec3D relative = target;
			Vec3D relativeVelocity = targetVelocity;

			if (request.center != EphemerisCenter::Heliocentric) {
				const auto earthState = earth.getMjd(row.mjdTdb);
				if (earth.error() != EphemeridesContainer::Error::NoError)
					throw std::runtime_error(earth.errorString());
				obs::ObserverState observer;
				if (site.ground)
					observer = obs::groundObserverStateEclAuMjd(row.mjdTdb, row.mjdUtc, site.data);

				if (result.lightTimeApplied) {
					const Vec3D observerPosition = earthState.R + observer.offsetEclAu;
					const Vec3D observerVelocity = earthState.V + observer.velocityEclAuPerDay;
					const Vec3D acceleration =
					    nBody ? propagator.accelerationAtMjd(target, targetVelocity, row.mjdTdb, solution.nonGrav)
					          : target * (-Constants::GM_SUN / std::pow(target.magnitude(), 3));
					const StateVectors initial(target, targetVelocity / Constants::GAUSS_K);
					const StateVectors initialRate(targetVelocity, acceleration / Constants::GAUSS_K);
					const Vec3D ray = target - observerPosition;
					const double firstDelay = ray.magnitude() / Constants::AU_PER_DAY;
					const double firstRate =
					    ray.dot(targetVelocity - observerVelocity) / (ray.magnitude() * Constants::AU_PER_DAY);
					const auto first = Calculations::propagateKeplerState(initial, -firstDelay);
					if (!first.isValid())
						throw std::runtime_error("Light-time first iteration failed");
					const Vec3D firstPositionRate =
					    Calculations::propagateKeplerPositionRate(initial, -firstDelay, initialRate, -firstRate);
					const Vec3D secondRay = first.first - observerPosition;
					const double secondDelay = secondRay.magnitude() / Constants::AU_PER_DAY;
					const double secondRate = secondRay.dot(firstPositionRate - observerVelocity) /
					                          (secondRay.magnitude() * Constants::AU_PER_DAY);
					const auto apparent = obs::solveLightTimeTwoPass(
					    time::jdFromMjd(row.mjdTdb), target, targetVelocity, earthState.R + observer.offsetEclAu);
					if (!apparent.valid)
						throw std::runtime_error("Light-time solution failed");
					target = apparent.position;
					targetVelocity =
					    Calculations::propagateKeplerPositionRate(initial, -secondDelay, initialRate, -secondRate);
				}

				relative = target - earthState.R - observer.offsetEclAu;
				relativeVelocity = targetVelocity - earthState.V - observer.velocityEclAuPerDay;
				const Vec3D observerHelio = earthState.R + observer.offsetEclAu;
				row.phase = angleBetween(target, relative);
				row.elongation = angleBetween(observerHelio * -1.0, relative);
			}

			row.position = toPlane(relative, request.plane, cosObliquity, sinObliquity);
			row.velocity = toPlane(relativeVelocity, request.plane, cosObliquity, sinObliquity);
			row.distanceAu = row.position.magnitude();
			row.heliocentricAu = target.magnitude();

			const Vec3D& p = row.position;
			const Vec3D& v = row.velocity;
			double longitude = std::atan2(p.y, p.x);
			if (longitude < 0.0)
				longitude += 2.0 * Constants::PI;
			row.longitude = longitude;
			row.latitude = std::asin(std::clamp(p.z / row.distanceAu, -1.0, 1.0));
			const double rho2 = p.x * p.x + p.y * p.y;
			if (rho2 > 0.0) {
				const double rho = std::sqrt(rho2);
				const double r2 = row.distanceAu * row.distanceAu;
				row.longitudeRate = (p.x * v.y - p.y * v.x) / rho2 * std::cos(row.latitude);
				row.latitudeRate = (v.z * r2 - p.z * p.dot(v)) / (r2 * rho);
			}

			if (solution.absoluteMagnitude && request.center != EphemerisCenter::Heliocentric &&
			    solution.magnitudeType == EphemerisMagnitudeType::Asteroid)
				row.visualMagnitude = photometry::visualMagnitude(*solution.absoluteMagnitude, solution.slopeG,
				                                                  row.heliocentricAu, row.distanceAu, row.phase);
			else if (solution.absoluteMagnitude && request.center != EphemerisCenter::Heliocentric)
				row.visualMagnitude =
				    photometry::cometMagnitude(*solution.absoluteMagnitude, solution.cometSlope, row.heliocentricAu,
				                               row.distanceAu, row.phase, solution.cometPhaseCoefficient);
		}

		if (control.cancel.stop_requested())
			return cancelled();
		result.rows = std::move(rows);
		return result;
	}

}  // namespace orbital_core
