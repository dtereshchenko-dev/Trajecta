#include "NBodyPropagator.h"
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <cmath>
#include <cstdint>
#include "integrators/IAS15Integrator.h"
#include "integrators/RA15Integrator.h"
#include "integrators/BarycentricIntegrator.h"
#include "propagator/EncounterStepControl.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/TimeScale.h"

#include <algorithm>
#include <array>
#include <limits>
#include <vector>
#include <stdexcept>
#include <type_traits>

namespace {
	std::vector<double> mjdsFromJds(const std::vector<double>& jds)
	{
		std::vector<double> mjds;
		mjds.reserve(jds.size());
		for (double jd : jds)
			mjds.push_back(orbital_core::time::mjdFromJd(jd));
		return mjds;
	}

	void accumulateIntegrationReport(IntegrationReport* aggregate, const IntegrationReport& segment)
	{
		if (!aggregate)
			return;

		aggregate->finalState = segment.finalState;
		aggregate->reachedMjd = segment.reachedMjd;
		aggregate->completedOutputs += segment.completedOutputs;
		aggregate->acceptedSteps += segment.acceptedSteps;
		aggregate->rejectedSteps += segment.rejectedSteps;
		aggregate->rhsEvaluations += segment.rhsEvaluations;
		aggregate->correctorIterations += segment.correctorIterations;
		aggregate->encounterSteps += segment.encounterSteps;
		aggregate->encounterRestarts += segment.encounterRestarts;
		aggregate->unconvergedEncounterSteps += segment.unconvergedEncounterSteps;
		aggregate->stoppedByObserver = aggregate->stoppedByObserver || segment.stoppedByObserver;
		aggregate->cancelled = aggregate->cancelled || segment.cancelled;
	}

	template <typename InternalState, typename PublicState, typename Integrate, typename Convert>
	std::vector<PublicState> integrateTargetBatches(const InternalState& initial, double epoch,
	                                                const std::vector<double>& targetJds,
	                                                IntegrationReport* aggregateReport, Integrate&& integrate,
	                                                Convert&& convert)
	{
		if (aggregateReport)
			*aggregateReport = {};

		std::vector<PublicState> results(targetJds.size());
		if (targetJds.empty())
			return results;

		std::vector<std::size_t> forward;
		std::vector<std::size_t> backward;
		forward.reserve(targetJds.size());
		backward.reserve(targetJds.size());
		for (std::size_t index = 0; index < targetJds.size(); ++index) {
			(targetJds[index] >= epoch ? forward : backward).push_back(index);
		}

		std::stable_sort(forward.begin(), forward.end(),
		                 [&targetJds](std::size_t lhs, std::size_t rhs) { return targetJds[lhs] < targetJds[rhs]; });
		std::stable_sort(backward.begin(), backward.end(),
		                 [&targetJds](std::size_t lhs, std::size_t rhs) { return targetJds[lhs] > targetJds[rhs]; });

		auto runDirection = [&](const std::vector<std::size_t>& indices) {
			if (indices.empty())
				return true;

			std::vector<double> orderedTargets;
			orderedTargets.reserve(indices.size());
			for (std::size_t index : indices)
				orderedTargets.push_back(targetJds[index]);

			std::vector<InternalState> integrated(indices.size());
			const IntegrationReport directionReport =
			    integrate(std::span<const double>(orderedTargets), std::span<InternalState>(integrated));
			accumulateIntegrationReport(aggregateReport, directionReport);
			if (directionReport.cancelled)
				return false;
			for (std::size_t offset = 0; offset < indices.size(); ++offset)
				results[indices[offset]] = convert(integrated[offset]);
			return true;
		};

		if (!runDirection(forward) || !runDirection(backward))
			return {};
		return results;
	}

	template <typename InternalState, typename PublicState, typename IntegrateSegment, typename Compose,
	          typename Convert>
	std::vector<PublicState> integrateChainedTargetBatches(const InternalState& initial, double epoch,
	                                                       const std::vector<double>& targetJds,
	                                                       IntegrationReport* aggregateReport,
	                                                       IntegrateSegment&& integrateSegment, Compose&& compose,
	                                                       Convert&& convert)
	{
		if (aggregateReport)
			*aggregateReport = {};

		std::vector<PublicState> results(targetJds.size());
		if (targetJds.empty())
			return results;

		std::vector<std::size_t> forward;
		std::vector<std::size_t> backward;
		forward.reserve(targetJds.size());
		backward.reserve(targetJds.size());
		for (std::size_t index = 0; index < targetJds.size(); ++index) {
			(targetJds[index] >= epoch ? forward : backward).push_back(index);
		}
		std::stable_sort(forward.begin(), forward.end(),
		                 [&targetJds](std::size_t lhs, std::size_t rhs) { return targetJds[lhs] < targetJds[rhs]; });
		std::stable_sort(backward.begin(), backward.end(),
		                 [&targetJds](std::size_t lhs, std::size_t rhs) { return targetJds[lhs] > targetJds[rhs]; });

		auto runDirection = [&](const std::vector<std::size_t>& indices) {
			InternalState accumulated = initial;
			double segmentEpoch = epoch;

			for (const std::size_t index : indices) {
				InternalState localInitial;
				localInitial.R = accumulated.R;
				localInitial.V = accumulated.V;

				const std::array<double, 1> target{targetJds[index]};
				std::array<InternalState, 1> localOutput;
				const IntegrationReport segmentReport = integrateSegment(
				    localInitial, segmentEpoch, std::span<const double>(target), std::span<InternalState>(localOutput));
				accumulateIntegrationReport(aggregateReport, segmentReport);

				accumulated = compose(localOutput.front(), accumulated);
				results[index] = convert(accumulated);
				segmentEpoch = target.front();
			}
		};

		runDirection(forward);
		runDirection(backward);
		return results;
	}
}  // namespace

NBodyPropagator::NBodyPropagator(IntegratorKind integrator) { setIntegratorKind(integrator); }

NBodyPropagator::NBodyPropagator(const NBodyPropagator& other) : NBodyPropagator(other.m_integratorKind)
{
	m_dynamics = other.m_dynamics;
	m_nonGrav = other.m_nonGrav;
	m_stepSize = other.m_stepSize;
	m_encounterEphemeris = other.m_encounterEphemeris;
	m_encounterPoints = other.m_encounterPoints;
	m_earthEncounterRadiusAu = other.m_earthEncounterRadiusAu;
}

NBodyPropagator& NBodyPropagator::operator=(const NBodyPropagator& other)
{
	if (this != &other)
		*this = NBodyPropagator(other);

	return *this;
}

NBodyPropagator NBodyPropagator::forEncounterMonitoring(IEphemerisSource& eph, int npoint,
                                                        double earthEncounterRadiusAu) const
{
	if (npoint < 0 || !(earthEncounterRadiusAu > 0.0))
		throw std::invalid_argument("Invalid encounter step configuration");
	NBodyPropagator out = *this;
	out.m_encounterEphemeris = &eph;
	out.m_encounterPoints = npoint;
	out.m_earthEncounterRadiusAu = earthEncounterRadiusAu;
	return out;
}

IntegrationOptions NBodyPropagator::integrationOptions() const
{
	IntegrationOptions options;
	options.initialStepDays = m_stepSize;
	if (!m_encounterEphemeris || m_encounterPoints <= 0)
		return options;
	options.restartOnEncounterExit = true;

	struct Body {
		EphemerisPlanet planet;
		double radius;
	};

	std::vector<Body> bodies;
	const std::array candidates{Body{EphemerisPlanet::Mercury, 0.05},
	                            Body{EphemerisPlanet::Venus, 0.05},
	                            Body{EphemerisPlanet::Earth, m_earthEncounterRadiusAu},
	                            Body{EphemerisPlanet::Mars, 0.05},
	                            Body{EphemerisPlanet::Jupiter, 0.7},
	                            Body{EphemerisPlanet::Saturn, 0.7},
	                            Body{EphemerisPlanet::Uranus, 0.7},
	                            Body{EphemerisPlanet::Neptune, 0.7},
	                            Body{EphemerisPlanet::Moon, 0.01},
	                            Body{EphemerisPlanet::Pluto, 0.05}};
	for (const auto& body : candidates) {
		if (planetEnabled(body.planet))
			bodies.push_back(body);
	}
	SpkReader* asteroids = useSb441() ? getDefaultAsteroidSource() : nullptr;
	auto asteroidIds = asteroids && asteroids->isAvailable() ? asteroids->targetIds() : std::vector<int>{};
	std::erase(asteroidIds, targetNaifId());
	auto controller = [eph = m_encounterEphemeris, points = m_encounterPoints, bodies = std::move(bodies), asteroids,
	                   asteroidIds, active = false, cap = 0.0](
	                      double step, double mjd, const IntegratedState& state) mutable -> std::optional<double> {
		bool found = false;
		const double direction = std::copysign(1.0, step);
		double next = std::numeric_limits<double>::infinity();
		const auto constrain = [&](const Vec3D& r, const Vec3D& v, double radius) {
			if (r.dot(r) >= radius * radius)
				return;
			found = true;
			next = std::min(next, std::abs(orbital_core::encounterStepDays(r, v, points, direction)));
		};
		for (const auto& body : bodies) {
			const auto p = eph->queryMjd(body.planet, mjd);
			constrain(state.R - p.R, state.V - p.V, body.radius);
		}
		for (const int id : asteroidIds) {
			Vec3D r, v;
			if (!asteroids->queryStateMjd(id, mjd, r, v))
				continue;
			constrain(state.R - r, state.V - v, 0.02);
		}
		if (!found) {
			active = false;
			return std::nullopt;
		}
		if (!active) {
			active = true;
			cap = std::abs(step);
		}
		return direction * std::min(cap, next);
	};
	options.nextEncounterStepOffset = [controller](double step, double epochMjd, double offset,
	                                               const IntegratedState& state) mutable {
		return controller(step, epochMjd + offset, state);
	};
	options.nextEncounterStep = [controller](double step, double epochMjd, const IntegratedState& state) mutable {
		return controller(step, epochMjd, state);
	};
	return options;
}

IntegratorKind NBodyPropagator::integratorKind() const noexcept { return m_integratorKind; }

void NBodyPropagator::setIntegratorKind(IntegratorKind kind)
{
	switch (kind) {
		case IntegratorKind::RA15:
			m_integrator = std::make_unique<RA15Integrator>();
			break;
		case IntegratorKind::IAS15:
			m_integrator = std::make_unique<IAS15Integrator>();
			break;
		default:
			throw std::invalid_argument("Unsupported integrator kind");
	}
	m_integrator = std::make_unique<BarycentricIntegrator>(std::move(m_integrator));
	m_integratorKind = kind;
}

double NBodyPropagator::fixedAccuracyParameter() const noexcept
{
	return m_integratorKind == IntegratorKind::IAS15 ? IAS15Integrator::kEpsilon : RA15Integrator::kAccuracy;
}

NBodyPropagator::PropagationWithStmSurfaceStop NBodyPropagator::propagateWithStmToSurface(
    const Vec3D& r0, const Vec3D& v0, double epoch, std::span<const double> checkpointJds, EphemerisPlanet planet,
    double radiusKm, IEphemerisSource& eph, const NonGravModel& ng, IntegrationReport* report) const
{
	if (!(radiusKm > 0.0) || !std::isfinite(radiusKm))
		throw std::invalid_argument("Invalid surface radius");
	if (report)
		*report = {};
	const double radiusAu = radiusKm * 1000.0 / Constants::AU_IN_METERS;
	auto run = [&]<class State>() {
		State accumulated;
		accumulated.R = r0;
		accumulated.V = v0;
		PropagationWithStmSurfaceStop result;
		result.mjd = orbital_core::time::mjdFromJd(epoch);
		result.jd = orbital_core::time::jdFromMjd(result.mjd);
		std::size_t completedCheckpoints = 0;
		const auto copyState = [](const IntegratedStateFull& full) {
			State state;
			state.R = full.R;
			state.V = full.V;
			state.Phi = full.Phi;
			if constexpr (std::is_same_v<State, IntegratedStateFull>)
				state.Sens = full.Sens;
			return state;
		};
		const auto compose = [](State local, const State& previous) {
			if constexpr (std::is_same_v<State, IntegratedStateFull>)
				local.Sens = (local.Phi * previous.Sens + local.Sens).eval();
			local.Phi = (local.Phi * previous.Phi).eval();
			return local;
		};
		const auto offset = [&](const State& state, double mjd) {
			return (state.R - eph.queryMjd(planet, mjd).R).magnitude() - radiusAu;
		};
		if (offset(accumulated, result.mjd) <= 0.0)
			throw std::invalid_argument("Surface entry requires an exterior initial state");
		for (double targetJd : checkpointJds) {
			const double target = orbital_core::time::mjdFromJd(targetJd);
			if (!std::isfinite(target) || target < result.mjd)
				throw std::invalid_argument("Surface checkpoints must be forward and monotonic");
			State initial;
			initial.R = accumulated.R;
			initial.V = accumulated.V;
			State outside = initial;
			double outsideMjd = result.mjd;
			const auto advance = [&](const State& start, double from, double to, const IntegrationOptions& options,
			                         State& end) {
				const std::array<double, 1> times{to};
				std::array<State, 1> states;
				IntegrationReport segment;
				if constexpr (std::is_same_v<State, IntegratedStateFull>)
					segment =
					    m_integrator->integrateStateWithStmAndSens(m_dynamics, ng, start, from, times, states, options);
				else
					segment = m_integrator->integrateStateWithStm(m_dynamics, ng, start, from, times, states, options);
				accumulateIntegrationReport(report, segment);
				if (!segment.stoppedByObserver && segment.completedOutputs != 1)
					throw std::runtime_error("Surface propagation did not reach checkpoint");
				if (!segment.stoppedByObserver)
					end = states.front();
				return segment;
			};
			auto options = integrationOptions();
			options.variationalObserverOffset = [&](double segmentEpochMjd, double offsetDays,
			                                        const IntegratedStateFull& state) {
				const double mjd = segmentEpochMjd + offsetDays;
				const auto p = eph.queryMjd(planet, mjd);
				const auto r = state.R - p.R;
				if (r.dot(r) >= radiusAu * radiusAu) {
					outside = copyState(state);
					outsideMjd = mjd;
					return true;
				}
				result.mjd = mjd;
				result.hitSurface = true;
				return false;
			};
			State local;
			const auto segment = advance(initial, result.mjd, target, options, local);
			if (result.hitSurface) {
				State start;
				start.R = outside.R;
				start.V = outside.V;
				const auto refineOptions = integrationOptions();
				const auto atOffset = [&](double dt) {
					State refined;
					advance(start, outsideMjd, outsideMjd + dt, refineOptions, refined);
					return compose(refined, outside);
				};
				double lo = 0.0, hi = result.mjd - outsideMjd;
				double flo = offset(outside, outsideMjd);
				local = atOffset(hi);
				double fhi = offset(local, outsideMjd + hi);
				for (int attempt = 0; fhi > 0.0 && attempt < 8; ++attempt) {
					result.mjd = std::nextafter(result.mjd, std::numeric_limits<double>::infinity());
					hi = result.mjd - outsideMjd;
					local = atOffset(hi);
					fhi = offset(local, result.mjd);
				}
				if (!(flo >= 0.0 && fhi <= 0.0))
					throw std::runtime_error("Surface refinement lost the incoming bracket");
				double bestOffset = fhi;
				if (std::abs(flo) < std::abs(fhi)) {
					local = outside;
					result.mjd = outsideMjd;
					bestOffset = flo;
				}
				int previousSide = 0;
				for (int iteration = 0; iteration < 64; ++iteration) {
					if (std::abs(bestOffset) * Constants::AU_IN_METERS / 1000.0 <= 1.0e-7)
						break;
					double trial = lo + (hi - lo) * flo / (flo - fhi);
					if (!(trial > lo && trial < hi))
						trial = lo + (hi - lo) * 0.5;
					double mjd = outsideMjd + trial;
					if (!(mjd > outsideMjd + lo && mjd < outsideMjd + hi)) {
						mjd = outsideMjd + lo + (hi - lo) * 0.5;
						if (!(mjd > outsideMjd + lo && mjd < outsideMjd + hi))
							break;
					}
					trial = mjd - outsideMjd;
					const State refined = atOffset(trial);
					const double value = offset(refined, mjd);
					if (std::abs(value) < std::abs(bestOffset)) {
						local = refined;
						result.mjd = mjd;
						bestOffset = value;
					}
					const int side = value > 0.0 ? 1 : -1;
					if (side > 0) {
						lo = trial;
						flo = value;
						if (previousSide == side)
							fhi *= 0.5;
					} else {
						hi = trial;
						fhi = value;
						if (previousSide == side)
							flo *= 0.5;
					}
					previousSide = side;
				}
				const double timeQuantum =
				    std::nextafter(result.mjd, std::numeric_limits<double>::infinity()) - result.mjd;
				const double speed = (local.V - eph.queryMjd(planet, result.mjd).V).magnitude();
				const double spatialFloor =
				    4.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, local.R.magnitude());
				if (std::abs(bestOffset) > speed * timeQuantum + spatialFloor)
					throw std::runtime_error("Surface crossing did not converge to MJD resolution");
			} else {
				++completedCheckpoints;
			}
			if constexpr (std::is_same_v<State, IntegratedStateFull>)
				accumulated.Sens = (local.Phi * accumulated.Sens + local.Sens).eval();
			accumulated.Phi = (local.Phi * accumulated.Phi).eval();
			accumulated.R = local.R;
			accumulated.V = local.V;
			if (!result.hitSurface)
				result.mjd = segment.reachedMjd;
			result.jd = orbital_core::time::jdFromMjd(result.mjd);
			if (result.hitSurface)
				break;
		}
		result.state.R = accumulated.R;
		result.state.V = accumulated.V;
		result.state.Phi = accumulated.Phi;
		result.state.SensA.setZero();
		if constexpr (std::is_same_v<State, IntegratedStateFull>)
			result.state.SensA = accumulated.Sens;
		if (report) {
			report->reachedMjd = result.mjd;
			report->finalState = {accumulated.R, accumulated.V};
			report->completedOutputs = completedCheckpoints;
		}
		return result;
	};
	if (ng.numFitParams() > 0)
		return run.template operator()<IntegratedStateFull>();
	return run.template operator()<IntegratedStateWithStm>();
}

Vec3D NBodyPropagator::accelerationAtMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng) const
{
	return m_dynamics.accelerationAtMjd(r, v, mjd, ng);
}

Vec3D NBodyPropagator::accelerationAt(const Vec3D& r, const Vec3D& v, double jd) const
{
	return m_dynamics.accelerationAtJd(r, v, jd, m_nonGrav);
}

Vec3D NBodyPropagator::accelerationAt(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng) const
{
	return m_dynamics.accelerationAtJd(r, v, jd, ng);
}

void NBodyPropagator::setPlanetEnabled(EphemerisPlanet planet, bool enabled)
{
	m_dynamics.setPlanetEnabled(planet, enabled);
}

bool NBodyPropagator::planetEnabled(EphemerisPlanet planet) const { return m_dynamics.planetEnabled(planet); }

void NBodyPropagator::setAllPlanetsEnabled(bool enabled) { m_dynamics.setAllPlanetsEnabled(enabled); }

uint16_t NBodyPropagator::planetMask() const noexcept { return m_dynamics.bodyMask(); }

void NBodyPropagator::setPlanetMask(uint16_t mask) { m_dynamics.setBodyMask(mask); }

void NBodyPropagator::setUseSb441(bool enabled) { m_dynamics.setUseSb441(enabled); }

bool NBodyPropagator::useSb441() const noexcept { return m_dynamics.useSb441(); }

void NBodyPropagator::setTargetNaifId(int id) noexcept { m_dynamics.setTargetNaifId(id); }

int NBodyPropagator::targetNaifId() const noexcept { return m_dynamics.targetNaifId(); }

void NBodyPropagator::setInactiveInteriorMassPolicy(orbital_core::InactiveInteriorMassPolicy policy) noexcept
{
	m_dynamics.setInactiveInteriorMassPolicy(policy);
}

orbital_core::InactiveInteriorMassPolicy NBodyPropagator::inactiveInteriorMassPolicy() const noexcept
{
	return m_dynamics.inactiveInteriorMassPolicy();
}

NonGravModel& NBodyPropagator::nonGrav() { return m_nonGrav; }

const NonGravModel& NBodyPropagator::nonGrav() const { return m_nonGrav; }

void NBodyPropagator::setNonGravEnabled(bool enabled) { m_nonGrav.setEnabled(enabled); }

bool NBodyPropagator::nonGravEnabled() const { return m_nonGrav.enabled(); }

double NBodyPropagator::A1() const noexcept { return m_nonGrav.A1(); }

double NBodyPropagator::A2() const noexcept { return m_nonGrav.A2(); }

double NBodyPropagator::A3() const noexcept { return m_nonGrav.A3(); }

Eigen::Matrix3d NBodyPropagator::accelerationPartialR(const Vec3D& r, const Vec3D& v, double jd) const
{
	return m_dynamics.accelerationPartialR(r, v, jd, m_nonGrav);
}

Eigen::Matrix3d NBodyPropagator::accelerationPartialV(const Vec3D& r, const Vec3D& v, double jd) const
{
	return m_dynamics.accelerationPartialV(r, v, jd, m_nonGrav);
}

void NBodyPropagator::accelerationAndPartials(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng,
                                              Vec3D& a_out, Eigen::Matrix3d& dAdR_out, Eigen::Matrix3d& dAdV_out) const
{
	m_dynamics.accelerationAndPartialsJd(r, v, jd, ng, a_out, dAdR_out, dAdV_out);
}

const char* NBodyPropagator::integratorName() const noexcept { return m_integrator ? m_integrator->name() : ""; }

std::vector<NBodyPropagator::PropagationResultFull> NBodyPropagator::propagateAllWithStmAndSens(
    const Vec3D& R0, const Vec3D& V0, double epoch, const std::vector<double>& targetJDs) const
{
	return propagateAllWithStmAndSens(R0, V0, epoch, targetJDs, m_nonGrav);
}

std::vector<NBodyPropagator::PropagationResultFull> NBodyPropagator::propagateAllWithStmAndSens(
    const Vec3D& R0, const Vec3D& V0, double epoch, const std::vector<double>& targetJDs, const NonGravModel& ng,
    IntegrationReport* report) const
{
	IntegratedStateFull initial;
	initial.R = R0;
	initial.V = V0;
	const double epochMjd = orbital_core::time::mjdFromJd(epoch);
	const std::vector<double> targetMjds = mjdsFromJds(targetJDs);

	return integrateTargetBatches<IntegratedStateFull, PropagationResultFull>(
	    initial, epochMjd, targetMjds, report,
	    [this, &ng, &initial, epochMjd](std::span<const double> orderedTargets,
	                                    std::span<IntegratedStateFull> integrated) {
		    return m_integrator->integrateStateWithStmAndSens(m_dynamics, ng, initial, epochMjd, orderedTargets,
		                                                      integrated, integrationOptions());
	    },
	    [](const IntegratedStateFull& state) {
		    return PropagationResultFull{state.R, state.V, state.Phi, state.Sens};
	    });
}

std::vector<NBodyPropagator::PropagationResult> NBodyPropagator::propagateAllWithStm(
    const Vec3D& R0, const Vec3D& V0, double epoch, const std::vector<double>& targetJDs) const
{
	return propagateAllWithStm(R0, V0, epoch, targetJDs, m_nonGrav);
}

std::vector<NBodyPropagator::PropagationResult> NBodyPropagator::propagateAllWithStm(
    const Vec3D& R0, const Vec3D& V0, double epoch, const std::vector<double>& targetJDs, const NonGravModel& ng,
    IntegrationReport* report) const
{
	IntegratedStateWithStm initial;
	initial.R = R0;
	initial.V = V0;
	const double epochMjd = orbital_core::time::mjdFromJd(epoch);
	const std::vector<double> targetMjds = mjdsFromJds(targetJDs);

	return integrateTargetBatches<IntegratedStateWithStm, PropagationResult>(
	    initial, epochMjd, targetMjds, report,
	    [this, &ng, &initial, epochMjd](std::span<const double> orderedTargets,
	                                    std::span<IntegratedStateWithStm> integrated) {
		    return m_integrator->integrateStateWithStm(m_dynamics, ng, initial, epochMjd, orderedTargets, integrated,
		                                               integrationOptions());
	    },
	    [](const IntegratedStateWithStm& state) { return PropagationResult{state.R, state.V, state.Phi}; });
}

std::vector<NBodyPropagator::PropagationResult> NBodyPropagator::propagateAllWithStmChained(
    const Vec3D& R0, const Vec3D& V0, double epoch, const std::vector<double>& targetJDs, const NonGravModel& ng,
    IntegrationReport* report) const
{
	IntegratedStateWithStm initial;
	initial.R = R0;
	initial.V = V0;
	const double epochMjd = orbital_core::time::mjdFromJd(epoch);
	const std::vector<double> targetMjds = mjdsFromJds(targetJDs);

	return integrateChainedTargetBatches<IntegratedStateWithStm, PropagationResult>(
	    initial, epochMjd, targetMjds, report,
	    [this, &ng](const IntegratedStateWithStm& localInitial, double segmentEpoch, std::span<const double> target,
	                std::span<IntegratedStateWithStm> output) {
		    return m_integrator->integrateStateWithStm(m_dynamics, ng, localInitial, segmentEpoch, target, output,
		                                               integrationOptions());
	    },
	    [](const IntegratedStateWithStm& local, const IntegratedStateWithStm& accumulated) {
		    IntegratedStateWithStm composed = local;
		    composed.Phi = local.Phi * accumulated.Phi;
		    return composed;
	    },
	    [](const IntegratedStateWithStm& state) { return PropagationResult{state.R, state.V, state.Phi}; });
}

std::vector<NBodyPropagator::PropagationResultFull> NBodyPropagator::propagateAllWithStmAndSensChained(
    const Vec3D& R0, const Vec3D& V0, double epoch, const std::vector<double>& targetJDs, const NonGravModel& ng,
    IntegrationReport* report) const
{
	IntegratedStateFull initial;
	initial.R = R0;
	initial.V = V0;
	const double epochMjd = orbital_core::time::mjdFromJd(epoch);
	const std::vector<double> targetMjds = mjdsFromJds(targetJDs);

	return integrateChainedTargetBatches<IntegratedStateFull, PropagationResultFull>(
	    initial, epochMjd, targetMjds, report,
	    [this, &ng](const IntegratedStateFull& localInitial, double segmentEpoch, std::span<const double> target,
	                std::span<IntegratedStateFull> output) {
		    return m_integrator->integrateStateWithStmAndSens(m_dynamics, ng, localInitial, segmentEpoch, target,
		                                                      output, integrationOptions());
	    },
	    [](const IntegratedStateFull& local, const IntegratedStateFull& accumulated) {
		    IntegratedStateFull composed = local;
		    composed.Sens = local.Phi * accumulated.Sens + local.Sens;
		    composed.Phi = local.Phi * accumulated.Phi;
		    return composed;
	    },
	    [](const IntegratedStateFull& state) {
		    return PropagationResultFull{state.R, state.V, state.Phi, state.Sens};
	    });
}

Vec3D NBodyPropagator::propagatePosition(const Vec3D& R0, const Vec3D& V0, double epoch, double targetJD) const
{
	return propagatePosition(R0, V0, epoch, targetJD, m_nonGrav);
}

Vec3D NBodyPropagator::propagatePosition(const Vec3D& R0, const Vec3D& V0, double epoch, double targetJD,
                                         const NonGravModel& ng) const
{
	const auto states = propagateAllWithVelocity(R0, V0, epoch, std::vector<double>{targetJD}, ng);
	return states.front().R;
}

std::vector<Vec3D> NBodyPropagator::propagateAll(const Vec3D& R0, const Vec3D& V0, double epoch,
                                                 const std::vector<double>& targetJDs) const
{
	return propagateAll(R0, V0, epoch, targetJDs, m_nonGrav);
}

std::vector<Vec3D> NBodyPropagator::propagateAll(const Vec3D& R0, const Vec3D& V0, double epoch,
                                                 const std::vector<double>& targetJDs, const NonGravModel& ng) const
{
	const auto states = propagateAllWithVelocity(R0, V0, epoch, targetJDs, ng);
	std::vector<Vec3D> results;
	results.reserve(states.size());
	for (const StateVel& state : states)
		results.push_back(state.R);
	return results;
}

std::vector<NBodyPropagator::StateVel> NBodyPropagator::propagateAllWithVelocity(
    const Vec3D& R0, const Vec3D& V0, double epoch, const std::vector<double>& targetJDs) const
{
	return propagateAllWithVelocity(R0, V0, epoch, targetJDs, m_nonGrav);
}

std::vector<NBodyPropagator::StateVel> NBodyPropagator::propagateAllWithVelocityMjd(
    const Vec3D& R0, const Vec3D& V0, double epochMjd, const std::vector<double>& targetMjds, const NonGravModel& ng,
    std::stop_token cancel) const
{
	if (cancel.stop_requested())
		return {};
	const IntegratedState initial{R0, V0};
	auto options = integrationOptions();
	options.cancel = cancel;
	return integrateTargetBatches<IntegratedState, StateVel>(
	    initial, epochMjd, targetMjds, nullptr,
	    [&](std::span<const double> targets, std::span<IntegratedState> output) {
		    return m_integrator->integrateState(m_dynamics, ng, initial, epochMjd, targets, output, options);
	    },
	    [](const IntegratedState& state) { return StateVel{state.R, state.V}; });
}

std::optional<IntegratedStateFull> NBodyPropagator::DenseTrajectory::at(double mjd) const
{
	if (mjd == epochMjd)
		return initial;
	const double tolerance = 2.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(mjd));
	const auto it = std::lower_bound(steps.begin(), steps.end(), mjd, [tolerance](const DenseIntegratedStep& step, double t) {
		return std::max(step.beginMjd, step.endMjd) + tolerance < t;
	});
	if (it == steps.end() || mjd < std::min(it->beginMjd, it->endMjd) - tolerance || !it->evaluate)
		return std::nullopt;
	return it->evaluate(mjd);
}

NBodyPropagator::DenseTrajectory NBodyPropagator::propagateDenseMjd(
    const Vec3D& R0, const Vec3D& V0, double epochMjd, double firstMjd, double lastMjd, const NonGravModel& ng,
    bool variational, bool sensitivities) const
{
	if (!std::isfinite(epochMjd) || !std::isfinite(firstMjd) || !std::isfinite(lastMjd) || firstMjd > lastMjd ||
	    (sensitivities && !variational))
		throw std::invalid_argument("Invalid dense trajectory interval");
	DenseTrajectory result;
	result.epochMjd = epochMjd;
	result.initial.R = R0;
	result.initial.V = V0;
	if (firstMjd < epochMjd)
		integrateDenseSegment(result, result.initial, epochMjd, firstMjd, ng, variational, sensitivities);
	if (lastMjd > epochMjd)
		integrateDenseSegment(result, result.initial, epochMjd, lastMjd, ng, variational, sensitivities);
	std::sort(result.steps.begin(), result.steps.end(), [](const DenseIntegratedStep& a, const DenseIntegratedStep& b) {
		return std::min(a.beginMjd, a.endMjd) < std::min(b.beginMjd, b.endMjd);
	});
	return result;
}

void NBodyPropagator::integrateDenseSegment(DenseTrajectory& trajectory, const IntegratedStateFull& initial,
                                            double epochMjd, double targetMjd, const NonGravModel& ng,
                                            bool variational, bool sensitivities) const
{
	auto options = integrationOptions();
	options.denseStepObserver = [&trajectory](DenseIntegratedStep step) { trajectory.steps.push_back(std::move(step)); };
	const std::array targets{targetMjd};
	if (sensitivities) {
		std::array<IntegratedStateFull, 1> output;
		m_integrator->integrateStateWithStmAndSens(m_dynamics, ng, initial, epochMjd, targets, output, options);
	} else if (variational) {
		IntegratedStateWithStm state{initial.R, initial.V, initial.Phi};
		std::array<IntegratedStateWithStm, 1> output;
		m_integrator->integrateStateWithStm(m_dynamics, ng, state, epochMjd, targets, output, options);
	} else {
		std::array<IntegratedState, 1> output;
		m_integrator->integrateState(m_dynamics, ng, {initial.R, initial.V}, epochMjd, targets, output, options);
	}
}

void NBodyPropagator::extendDenseBackwardMjd(DenseTrajectory& trajectory, double firstMjd, const NonGravModel& ng,
                                             bool variational, bool sensitivities) const
{
	const double current = trajectory.steps.empty() ? trajectory.epochMjd :
	    std::min(trajectory.steps.front().beginMjd, trajectory.steps.front().endMjd);
	if (!std::isfinite(firstMjd) || firstMjd >= current || (sensitivities && !variational))
		throw std::invalid_argument("Invalid dense trajectory extension");
	const auto state = trajectory.at(current);
	if (!state)
		throw std::runtime_error("Dense trajectory has no initial state for extension");
	integrateDenseSegment(trajectory, *state, current, firstMjd, ng, variational, sensitivities);
	std::sort(trajectory.steps.begin(), trajectory.steps.end(), [](const DenseIntegratedStep& a, const DenseIntegratedStep& b) {
		return std::min(a.beginMjd, a.endMjd) < std::min(b.beginMjd, b.endMjd);
	});
}

std::vector<NBodyPropagator::PropagationResult> NBodyPropagator::propagateAllWithStmMjd(
    const Vec3D& R0, const Vec3D& V0, double epochMjd, const std::vector<double>& targetMjds, const NonGravModel& ng,
    IntegrationReport* report) const
{
	IntegratedStateWithStm initial;
	initial.R = R0;
	initial.V = V0;
	return integrateTargetBatches<IntegratedStateWithStm, PropagationResult>(
	    initial, epochMjd, targetMjds, report,
	    [&](std::span<const double> targets, std::span<IntegratedStateWithStm> output) {
		    auto segment = m_integrator->integrateStateWithStm(m_dynamics, ng, initial, epochMjd, targets, output,
		                                                       integrationOptions());
		    return segment;
	    },
	    [](const IntegratedStateWithStm& state) { return PropagationResult{state.R, state.V, state.Phi}; });
}

std::vector<NBodyPropagator::PropagationResultFull> NBodyPropagator::propagateAllWithStmAndSensMjd(
    const Vec3D& R0, const Vec3D& V0, double epochMjd, const std::vector<double>& targetMjds, const NonGravModel& ng,
    IntegrationReport* report, std::stop_token cancel) const
{
	if (cancel.stop_requested())
		return {};
	IntegratedStateFull initial;
	initial.R = R0;
	initial.V = V0;
	auto options = integrationOptions();
	options.cancel = cancel;
	return integrateTargetBatches<IntegratedStateFull, PropagationResultFull>(
	    initial, epochMjd, targetMjds, report,
	    [&](std::span<const double> targets, std::span<IntegratedStateFull> output) {
		    auto segment =
		        m_integrator->integrateStateWithStmAndSens(m_dynamics, ng, initial, epochMjd, targets, output, options);
		    return segment;
	    },
	    [](const IntegratedStateFull& state) {
		    return PropagationResultFull{state.R, state.V, state.Phi, state.Sens};
	    });
}

std::vector<NBodyPropagator::StateVel> NBodyPropagator::propagateAllWithVelocity(const Vec3D& R0, const Vec3D& V0,
                                                                                 double epoch,
                                                                                 const std::vector<double>& targetJDs,
                                                                                 const NonGravModel& ng,
                                                                                 IntegrationReport* report) const
{
	IntegratedState initial{R0, V0};
	const double epochMjd = orbital_core::time::mjdFromJd(epoch);
	const std::vector<double> targetMjds = mjdsFromJds(targetJDs);
	return integrateTargetBatches<IntegratedState, StateVel>(
	    initial, epochMjd, targetMjds, report,
	    [this, &ng, &initial, epochMjd](std::span<const double> orderedTargets, std::span<IntegratedState> integrated) {
		    return m_integrator->integrateState(m_dynamics, ng, initial, epochMjd, orderedTargets, integrated,
		                                        integrationOptions());
	    },
	    [](const IntegratedState& state) { return StateVel{state.R, state.V}; });
}
