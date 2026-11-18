#include "orbital_core/RelativeMotion.h"

#include "NBodyPropagator.h"
#include "RelativeMotionInternal.h"
#include "TwoBodyPropagation.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "orbital_core/perturbations/PlanetaryPerturbations.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core {
	NBodyPropagator detail::relativeMotionPropagator(const NBodySettings& forces, IEphemerisSource& source)
	{
		NBodyPropagator propagator;
		propagator.setStepSize(forces.stepSize);
		propagator.setPlanetMask(forces.planetMask);
		propagator.setUseSb441(forces.useSb441);
		propagator.setForceModel(forces.forceModel);
		propagator.setTargetNaifId(forces.targetNaifId);
		propagator.setInactiveInteriorMassPolicy(forces.inactiveInteriorMassPolicy);
		return propagator.forEncounterMonitoring(source);
	}

	namespace {
		double planetMu(EphemerisPlanet planet)
		{
			for (const auto& body : forces::detail::kPerturbers) {
				if (body.planet == planet && planet != EphemerisPlanet::Moon && planet != EphemerisPlanet::Pluto)
					return body.gm;
			}
			throw std::invalid_argument("Relative motion requires a planet from Mercury through Neptune");
		}

		bool finite(const Vec3D& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
	}  // namespace

	std::vector<RelativeMotionSample> computeRelativeMotion(const RelativeMotionRequest& request,
	                                                        const FitRunControl& control)
	{
		const auto& solution = request.solution;
		const double mu = planetMu(request.planet);
		if (!finite(solution.position) || !finite(solution.velocity) || !std::isfinite(solution.epochMjdTdb) ||
		    request.epochsMjdTdb.size() < 2)
			throw std::invalid_argument("Invalid relative motion state or number of epochs");
		for (std::size_t i = 0; i < request.epochsMjdTdb.size(); ++i) {
			if (!std::isfinite(request.epochsMjdTdb[i]) ||
			    (i && request.epochsMjdTdb[i] <= request.epochsMjdTdb[i - 1]))
				throw std::invalid_argument("Relative motion epochs must be finite and strictly increasing");
		}
		auto* source = getDefaultEphemerisSource();
		if (!source || !source->isAvailable())
			throw std::runtime_error("Planetary ephemeris is unavailable");
		if (control.cancel.stop_requested())
			return {};

		std::vector<NBodyPropagator::StateVel> states;
		if (solution.forces.enabled) {
			auto propagator = detail::relativeMotionPropagator(solution.forces, *source);
			states = propagator.propagateAllWithVelocityMjd(solution.position, solution.velocity, solution.epochMjdTdb,
			                                                request.epochsMjdTdb, solution.nonGrav, control.cancel);
		} else {
			const StateVectors initial(solution.position, solution.velocity / Constants::GAUSS_K);
			for (double epoch : request.epochsMjdTdb) {
				if (control.cancel.stop_requested())
					return {};
				const auto state = Calculations::propagateKeplerState(initial, epoch - solution.epochMjdTdb);
				if (!state.isValid())
					throw std::runtime_error("Two-body propagation failed");
				states.push_back({state.first, state.second * Constants::GAUSS_K});
			}
		}
		if (control.cancel.stop_requested())
			return {};
		if (states.size() != request.epochsMjdTdb.size())
			throw std::runtime_error("Relative motion propagation returned incomplete states");

		std::vector<RelativeMotionSample> result;
		result.reserve(states.size());
		for (std::size_t i = 0; i < states.size(); ++i) {
			if (control.cancel.stop_requested())
				return {};
			const double epoch = request.epochsMjdTdb[i];
			const auto planet = source->queryMjd(request.planet, epoch);
			const double radius = planet.R.magnitude();
			const Vec3D x = planet.R.normalized();
			const Vec3D z = planet.R.cross(planet.V).normalized();
			const Vec3D y = z.cross(x).normalized();
			if (!finite(planet.R) || !finite(planet.V) || radius <= 0 || z.isNull() || y.isNull() ||
			    !finite(states[i].R) || !finite(states[i].V))
				throw std::runtime_error("Invalid planet or object state in relative motion");
			const auto rotate = [&](const Vec3D& v) { return Vec3D(v.dot(x), v.dot(y), v.dot(z)); };
			RelativeMotionSample sample;
			sample.mjdTdb = epoch;
			sample.positionAu = states[i].R - planet.R;
			sample.velocityAuPerDay = states[i].V - planet.V;
			sample.rotatingPositionAu = rotate(sample.positionAu);
			sample.sunDirection = x * -1;
			sample.distanceAu = sample.positionAu.magnitude();
			if (!(sample.distanceAu > 0))
				throw std::runtime_error("Relative position is zero");
			sample.hillRadiusAu = radius * std::cbrt(mu / (3.0 * Constants::GM_SUN));
			sample.specificEnergyAu2PerDay2 =
			    sample.velocityAuPerDay.dot(sample.velocityAuPerDay) / 2 - mu / sample.distanceAu;
			if (request.planet == EphemerisPlanet::Earth) {
				const auto moon = source->queryMjd(EphemerisPlanet::Moon, epoch);
				if (!finite(moon.R) || moon.R.isNull())
					throw std::runtime_error("Lunar ephemeris unavailable");
				sample.moonPositionAu = moon.R - planet.R;
				sample.rotatingMoonPositionAu = rotate(sample.moonPositionAu);
				sample.hasMoon = true;
			}
			result.push_back(sample);
		}
		return result;
	}

	std::vector<RelativeMotionSample> computeRefinedRelativeMotion(const RelativeMotionRequest& request,
	                                                               const FitRunControl& control)
	{
		constexpr int maxRounds = 16;
		constexpr double maxTravelOverDistance = 0.08;
		const auto limit = request.maxRefinedSamples.value_or(request.epochsMjdTdb.max_size());
		if (request.epochsMjdTdb.size() > limit)
			throw std::length_error("Trajectory exceeds the point limit; shorten the interval or increase the step");
		auto samples = computeRelativeMotion(request, control);
		for (int round = 0; round < maxRounds && !samples.empty(); ++round) {
			if (control.cancel.stop_requested())
				return {};
			RelativeMotionRequest refined;
			refined.solution = request.solution;
			refined.planet = request.planet;
			if (samples.size() > refined.epochsMjdTdb.max_size() / 2)
				throw std::length_error("Too many trajectory points");
			refined.epochsMjdTdb.reserve(std::min(samples.size() * 2, limit));
			refined.epochsMjdTdb.push_back(samples.front().mjdTdb);
			const auto appendEpoch = [&](double epoch) {
				if (refined.epochsMjdTdb.size() >= limit)
					throw std::length_error("Trajectory refinement exceeds the point limit; shorten the interval");
				refined.epochsMjdTdb.push_back(epoch);
			};
			bool needsRefinement = false;
			for (std::size_t i = 1; i < samples.size(); ++i) {
				const auto& a = samples[i - 1];
				const auto& b = samples[i];
				const Vec3D chord = b.positionAu - a.positionAu;
				const double chord2 = chord.dot(chord);
				const double f = chord2 > 0 ? std::clamp(-a.positionAu.dot(chord) / chord2, 0.0, 1.0) : 0.0;
				const double distance = std::max((a.positionAu + chord * f).magnitude(), 1e-12);
				const double travel =
				    (b.mjdTdb - a.mjdTdb) * std::max(a.velocityAuPerDay.magnitude(), b.velocityAuPerDay.magnitude());
				if (travel > maxTravelOverDistance * distance) {
					const double midpoint = a.mjdTdb + (b.mjdTdb - a.mjdTdb) / 2;
					if (midpoint > a.mjdTdb && midpoint < b.mjdTdb) {
						appendEpoch(midpoint);
						needsRefinement = true;
					}
				}
				appendEpoch(b.mjdTdb);
			}
			if (!needsRefinement)
				return samples;
			samples = computeRelativeMotion(refined, control);
		}
		if (control.cancel.stop_requested())
			return {};
		throw std::runtime_error("Trajectory refinement did not finish; shorten the interval");
	}

	std::vector<HillSpherePassage> findHillSpherePassages(const std::vector<RelativeMotionSample>& samples)
	{
		std::vector<HillSpherePassage> result;
		if (samples.size() < 2)
			return result;
		for (std::size_t i = 0; i < samples.size(); ++i) {
			const auto& s = samples[i];
			if (!std::isfinite(s.mjdTdb) || !std::isfinite(s.distanceAu) || s.distanceAu <= 0 ||
			    !std::isfinite(s.hillRadiusAu) || s.hillRadiusAu <= 0 || (i && s.mjdTdb <= samples[i - 1].mjdTdb))
				throw std::invalid_argument("Invalid Hill-sphere samples");
		}
		HillSpherePassage current;
		bool active = samples.front().distanceAu < samples.front().hillRadiusAu;
		if (active) {
			current.entryMjdTdb = samples.front().mjdTdb;
			current.entryClipped = true;
		}
		for (std::size_t i = 1; i < samples.size(); ++i) {
			const auto& a = samples[i - 1];
			const auto& b = samples[i];
			const double x = a.distanceAu - a.hillRadiusAu;
			const double y = b.distanceAu - b.hillRadiusAu;
			if (!active && x >= 0 && y < 0) {
				current = {};
				current.entryMjdTdb = a.mjdTdb + (b.mjdTdb - a.mjdTdb) * x / (x - y);
				current.entryClipped = current.entryMjdTdb == samples.front().mjdTdb;
				active = true;
			} else if (active && x <= 0 && y > 0) {
				current.exitMjdTdb = a.mjdTdb + (b.mjdTdb - a.mjdTdb) * x / (x - y);
				result.push_back(current);
				active = false;
			}
		}
		if (active) {
			current.exitMjdTdb = samples.back().mjdTdb;
			current.exitClipped = true;
			result.push_back(current);
		}
		return result;
	}

	std::vector<TemporaryCaptureInterval> findTemporaryCaptures(const std::vector<RelativeMotionSample>& samples)
	{
		std::vector<TemporaryCaptureInterval> result;
		if (samples.size() < 2)
			return result;
		for (std::size_t i = 0; i < samples.size(); ++i) {
			const auto& s = samples[i];
			if (!std::isfinite(s.mjdTdb) || !std::isfinite(s.specificEnergyAu2PerDay2) ||
			    !std::isfinite(s.distanceAu) || s.distanceAu <= 0 || !std::isfinite(s.hillRadiusAu) ||
			    s.hillRadiusAu <= 0 || !finite(s.rotatingPositionAu) || (i && s.mjdTdb <= samples[i - 1].mjdTdb))
				throw std::invalid_argument("Invalid temporary-capture samples");
		}
		TemporaryCaptureInterval current;
		bool active = false;
		const auto finish = [&] {
			if (active) {
				current.endClipped = current.endMjdTdb == samples.back().mjdTdb;
				result.push_back(current);
				active = false;
			}
		};
		for (std::size_t i = 1; i < samples.size(); ++i) {
			const auto& a = samples[i - 1];
			const auto& b = samples[i];
			double lo = 0, hi = 1;
			const auto intersectNegative = [&](double x, double y) {
				if (x >= 0 && y >= 0)
					return false;
				if (x >= 0)
					lo = std::max(lo, x / (x - y));
				if (y >= 0)
					hi = std::min(hi, x / (x - y));
				return hi > lo;
			};
			if (!intersectNegative(a.specificEnergyAu2PerDay2, b.specificEnergyAu2PerDay2) ||
			    !intersectNegative(a.distanceAu - 3 * a.hillRadiusAu, b.distanceAu - 3 * b.hillRadiusAu)) {
				finish();
				continue;
			}
			if (active && lo > 0)
				finish();
			const auto epoch = [&](double f) { return a.mjdTdb + f * (b.mjdTdb - a.mjdTdb); };
			if (!active) {
				current = {};
				current.revolutionsValid = true;
				current.startMjdTdb = epoch(lo);
				current.startClipped = i == 1 && lo == 0;
				active = true;
			}
			current.endMjdTdb = epoch(hi);
			const auto addTurns = [&](const Vec3D& from, const Vec3D& to, double& turns, double& travel, bool& valid) {
				const Vec3D p = from * (1 - lo) + to * lo;
				const Vec3D q = from * (1 - hi) + to * hi;
				if (!finite(p) || !finite(q) || std::hypot(p.x, p.y) <= p.magnitude() * 1e-10 ||
				    std::hypot(q.x, q.y) <= q.magnitude() * 1e-10) {
					valid = false;
					return;
				}
				const double angle = std::atan2(p.x * q.y - p.y * q.x, p.x * q.x + p.y * q.y);
				if (std::abs(angle) >= Constants::PI / 2)
					valid = false;
				turns += angle / (2 * Constants::PI);
				travel += std::abs(angle) / (2 * Constants::PI);
			};
			addTurns(a.rotatingPositionAu, b.rotatingPositionAu, current.signedRevolutions,
			         current.angularTravelRevolutions, current.revolutionsValid);
			addTurns(a.positionAu, b.positionAu, current.inertialSignedRevolutions,
			         current.inertialAngularTravelRevolutions, current.inertialRevolutionsValid);
			if (hi < 1)
				finish();
		}
		finish();
		return result;
	}

}  // namespace orbital_core
