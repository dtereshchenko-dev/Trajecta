#include "strewn_core/TrajectoryIntegrator.h"
#include "VectorMath.h"
#include "RunDeadline.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace strewn {
	namespace {

		using State = std::array<double, 7>;
		constexpr double a[7][7] = {{},
		                            {0.2},
		                            {3.0 / 40, 9.0 / 40},
		                            {44.0 / 45, -56.0 / 15, 32.0 / 9},
		                            {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729},
		                            {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656},
		                            {35.0 / 384, 0.0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84}};
		constexpr double e[7] = {71.0 / 57600,      0.0,        -71.0 / 16695, 71.0 / 1920,
		                         -17253.0 / 339200, 22.0 / 525, -1.0 / 40};

		struct Evaluation {
			State derivative;
			GeodeticPosition position;
			double relativeSpeed;
			double density;
		};

		void positive(double value, const char* name)
		{
			if (!std::isfinite(value) || value <= 0.0)
				throw std::invalid_argument(name);
		}

		void nonnegative(double value, const char* name)
		{
			if (!std::isfinite(value) || value < 0.0)
				throw std::invalid_argument(name);
		}

		State pack(const FlightState& state)
		{
			return {state.positionEcefM[0],
			        state.positionEcefM[1],
			        state.positionEcefM[2],
			        state.velocityEcefMps[0],
			        state.velocityEcefMps[1],
			        state.velocityEcefMps[2],
			        state.massKg};
		}

		FlightState unpack(const State& y) { return {{y[0], y[1], y[2]}, {y[3], y[4], y[5]}, y[6]}; }

		Evaluation evaluate(const State& y, const RunParameters& parameters, Shape shape,
		                    const Environment& environment, const PhysicsConfig& physics)
		{
			Evaluation result;
			result.position = geodetic({y[0], y[1], y[2]});
			const auto air = environment.air(result.position);
			auto wind = enuToEcef(air.windEnuMps, result.position.latitudeRad, result.position.longitudeRad);
			Vector3 relative;
			for (std::size_t i = 0; i < 3; ++i)
				relative[i] = y[i + 3] - parameters.windScale * wind[i];
			const double speed = detail::norm(relative);
			const double area = projectedArea(y[6], parameters.densityKgM3, shape);
			const double cd = parameters.dragScale * dragCoefficient(speed / std::sqrt(401.874 * air.temperatureK));
			const double drag = -0.5 * air.densityKgM3 * cd * area * speed / y[6];
			const double r2 = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
			const double r = std::sqrt(r2), mu = detail::earthMu / (r2 * r);
			const double zr = y[2] * y[2] / r2;
			const double j2 =
			    1.5 * detail::earthJ2 * detail::earthMu * detail::gravityRadius * detail::gravityRadius / (r2 * r2 * r);
			const double w = detail::omega;
			result.derivative = {
			    y[3],
			    y[4],
			    y[5],
			    -mu * y[0] - j2 * (1.0 - 5.0 * zr) * y[0] + 2.0 * w * y[4] + w * w * y[0] + drag * relative[0],
			    -mu * y[1] - j2 * (1.0 - 5.0 * zr) * y[1] - 2.0 * w * y[3] + w * w * y[1] + drag * relative[1],
			    -mu * y[2] - j2 * (3.0 - 5.0 * zr) * y[2] + drag * relative[2],
			    physics.phase == FlightPhase::Entry && speed >= physics.luminousSpeedMps
			        ? -0.5 * parameters.ablationKgJ * air.densityKgM3 * area * speed * speed * speed
			        : 0.0};
			result.relativeSpeed = speed;
			result.density = air.densityKgM3;
			return result;
		}

		std::array<double, 3> eventValues(const State& y, const Evaluation& evaluation, double strength,
		                                  const Environment& environment, const PhysicsConfig& physics)
		{
			const double pressure = evaluation.density * evaluation.relativeSpeed * evaluation.relativeSpeed *
			                        (physics.ramPressure == RamPressure::RhoV2 ? 1.0 : 0.5);
			const bool canBreak =
			    physics.phase == FlightPhase::Entry && evaluation.relativeSpeed >= physics.luminousSpeedMps;
			return {environment.clearance(evaluation.position),
			        canBreak ? strength - pressure : std::numeric_limits<double>::max(), y[6] - physics.minMassKg};
		}

		State interpolate(double fraction, double step, const State& y0, const State& f0, const State& y1,
		                  const State& f1)
		{
			const double s2 = fraction * fraction, s3 = s2 * fraction;
			State result;
			for (std::size_t i = 0; i < 7; ++i)
				result[i] = (2 * s3 - 3 * s2 + 1) * y0[i] + (s3 - 2 * s2 + fraction) * step * f0[i] +
				            (-2 * s3 + 3 * s2) * y1[i] + (s3 - s2) * step * f1[i];
			return result;
		}

	}  // namespace

	void validate(const PhysicsConfig& config, const IntegrationOptions& options)
	{
		if ((config.phase != FlightPhase::Entry && config.phase != FlightPhase::DarkFlight) ||
		    (config.ramPressure != RamPressure::RhoV2 && config.ramPressure != RamPressure::HalfRhoV2) ||
		    (config.fragmentShape != Shape::Sphere && config.fragmentShape != Shape::BrickAreaOnly &&
		     config.fragmentShape != Shape::ConeAreaOnly))
			throw std::invalid_argument("Unknown physical model");
		positive(config.minMassKg, "minMassKg must be positive");
		positive(config.luminousSpeedMps, "luminousSpeedMps must be positive");
		nonnegative(config.weibullAlpha, "weibullAlpha must be nonnegative");
		positive(config.strengthCapPa, "strengthCapPa must be positive");
		nonnegative(config.lateralCoefficient, "lateralCoefficient must be nonnegative");
		if (config.binsPerDecade[0] == 0 || config.binsPerDecade[1] == 0 || config.binsPerDecade[0] > 1000 ||
		    config.binsPerDecade[1] > 1000 || config.representativesPerBin == 0 || config.representativesPerBin > 1024)
			throw std::invalid_argument("Invalid fragmentation limits");
		positive(options.relativeTolerance, "relativeTolerance must be positive");
		positive(options.positionToleranceM, "positionToleranceM must be positive");
		positive(options.velocityToleranceMps, "velocityToleranceMps must be positive");
		positive(options.massToleranceKg, "massToleranceKg must be positive");
		positive(options.minStepS, "minStepS must be positive");
		positive(options.initialStepS, "initialStepS must be positive");
		positive(options.maxStepS, "maxStepS must be positive");
		positive(options.maxFlightTimeS, "maxFlightTimeS must be positive");
		if (options.minStepS > options.initialStepS || options.initialStepS > options.maxStepS ||
		    options.maxStepAttempts == 0)
			throw std::invalid_argument("Invalid integration step limits");
	}

	void validate(const RunParameters& p)
	{
		positive(p.massKg, "massKg must be positive");
		positive(p.densityKgM3, "densityKgM3 must be positive");
		positive(p.strengthPa, "strengthPa must be positive");
		nonnegative(p.ablationKgJ, "ablationKgJ must be nonnegative");
		nonnegative(p.windScale, "windScale must be nonnegative");
		nonnegative(p.dragScale, "dragScale must be nonnegative");
		if (!std::isfinite(p.dustFraction) || p.dustFraction < 0.0 || p.dustFraction > 1.0 ||
		    !std::isfinite(p.largestFraction) || p.largestFraction <= 0.0 || p.largestFraction >= 1.0 ||
		    !std::isfinite(p.beta) || p.beta <= 0.0 || p.beta >= 1.0)
			throw std::invalid_argument("Invalid fragment mass distribution");
	}

	double projectedArea(double mass, double density, Shape shape)
	{
		const double radius = std::pow(3.0 * mass / (4.0 * detail::pi * density), 1.0 / 3.0);
		if (shape == Shape::Sphere)
			return detail::pi * radius * radius;
		const double factor = shape == Shape::BrickAreaOnly ? 15.0 / std::pow(30.0, 2.0 / 3.0)
		                                                    : detail::pi * std::pow(3.0 / detail::pi, 2.0 / 3.0);
		return factor * std::pow(mass / density, 2.0 / 3.0);
	}

	double dragCoefficient(double mach)
	{
		constexpr double speeds[] = {0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4};
		constexpr double coefficients[] = {0.66, 0.66, 0.70, 0.78, 0.88, 1, 1.26, 1.24, 1.16};
		for (std::size_t i = 1; i < std::size(speeds); ++i)
			if (mach <= speeds[i])
				return coefficients[i - 1] +
				       (mach - speeds[i - 1]) / (speeds[i] - speeds[i - 1]) * (coefficients[i] - coefficients[i - 1]);
		return coefficients[8];
	}

	TrajectoryIntegrator::TrajectoryIntegrator(const Environment& environment, const PhysicsConfig& physics,
	                                           const IntegrationOptions& options)
	    : m_environment(environment), m_physics(physics), m_options(options)
	{
		validate(physics, options);
	}

	FlightResult TrajectoryIntegrator::advance(const Fragment& fragment, const RunParameters& parameters,
	                                           std::stop_token cancel) const
	{
		return advance(fragment, parameters, cancel, {});
	}

	FlightResult TrajectoryIntegrator::advance(const Fragment& fragment, const RunParameters& parameters,
	                                           std::stop_token cancel, const FlightObserver& observer,
	                                           std::chrono::steady_clock::time_point deadline) const
	{
		validate(parameters);
		if (!detail::finite(fragment.state.positionEcefM) || !detail::finite(fragment.state.velocityEcefMps) ||
		    !std::isfinite(fragment.timeS) || fragment.timeS < 0.0 || !std::isfinite(fragment.multiplicity) ||
		    fragment.multiplicity <= 0.0 || !std::isfinite(fragment.state.massKg) || fragment.state.massKg <= 0.0 ||
		    !std::isfinite(fragment.strengthPa) || fragment.strengthPa <= 0.0)
			throw std::invalid_argument("Invalid initial fragment");
		IntegrationStats stats;
		State y = pack(fragment.state);
		double time = fragment.timeS, step = m_options.initialStepS;
		const auto finish = [&](FlightEvent event) {
			if (observer)
				observer(time, unpack(y));
			if (const auto interrupted = detail::interruption(cancel, deadline))
				event = *interrupted;
			return FlightResult{unpack(y), time, event, stats};
		};
		if (observer)
			observer(time, fragment.state);
		if (const auto interrupted = detail::interruption(cancel, deadline))
			return finish(*interrupted);
		auto current = evaluate(y, parameters, fragment.shape, m_environment, m_physics);
		if (!std::all_of(current.derivative.begin(), current.derivative.end(),
		                 [](double v) { return std::isfinite(v); }))
			throw std::runtime_error("Non-finite initial derivative");
		++stats.evaluations;
		auto previousEvents = eventValues(y, current, fragment.strengthPa, m_environment, m_physics);
		if (previousEvents[2] <= 0.0)
			return finish(FlightEvent::MassCutoff);
		if (previousEvents[0] <= 0.0)
			return finish(FlightEvent::Ground);
		if (previousEvents[1] <= 0.0)
			return finish(FlightEvent::Breakup);
		const State atol{m_options.positionToleranceM,   m_options.positionToleranceM,   m_options.positionToleranceM,
		                 m_options.velocityToleranceMps, m_options.velocityToleranceMps, m_options.velocityToleranceMps,
		                 m_options.massToleranceKg};
		std::array<State, 7> stages;
		for (std::size_t attempt = 0; attempt < m_options.maxStepAttempts; ++attempt) {
			if (cancel.stop_requested())
				return finish(FlightEvent::Cancelled);
			if ((attempt & 63) == 0)
				if (const auto interrupted = detail::interruption(cancel, deadline))
					return finish(*interrupted);
			const double remaining = m_options.maxFlightTimeS - (time - fragment.timeS);
			if (remaining <= m_options.minStepS)
				return finish(FlightEvent::Timeout);
			step = std::min(step, remaining);
			if (step < m_options.minStepS || time + step == time)
				throw std::runtime_error("Integrator step underflow");
			stages[0] = current.derivative;
			State next;
			Evaluation candidate;
			bool valid = true;
			for (std::size_t i = 1; i < 7; ++i) {
				for (std::size_t c = 0; c < 7; ++c) {
					double sum = 0.0;
					for (std::size_t j = 0; j < i; ++j)
						sum += a[i][j] * stages[j][c];
					next[c] = y[c] + step * sum;
					valid = valid && std::isfinite(next[c]);
				}
				if (!valid || next[6] <= 0.0) {
					valid = false;
					break;
				}
				candidate = evaluate(next, parameters, fragment.shape, m_environment, m_physics);
				++stats.evaluations;
				stages[i] = candidate.derivative;
				if (!std::all_of(stages[i].begin(), stages[i].end(), [](double v) { return std::isfinite(v); })) {
					valid = false;
					break;
				}
			}
			if (!valid) {
				step *= 0.25;
				++stats.rejectedSteps;
				continue;
			}
			double error = 0.0;
			for (std::size_t c = 0; c < 7; ++c) {
				double sum = 0.0;
				for (std::size_t i = 0; i < 7; ++i)
					sum += e[i] * stages[i][c];
				const double component =
				    std::abs(step * sum) / (atol[c] + m_options.relativeTolerance * std::abs(y[c]));
				if (!std::isfinite(component))
					throw std::runtime_error("Non-finite integration error");
				error = std::max(error, component);
			}
			if (error > 1.0) {
				step *= std::max(0.2, 0.9 * std::pow(error, -0.2));
				++stats.rejectedSteps;
				continue;
			}
			++stats.acceptedSteps;
			const auto newEvents = eventValues(next, candidate, fragment.strengthPa, m_environment, m_physics);
			double earliest = 2.0;
			std::size_t eventIndex = 0;
			for (std::size_t event = 0; event < newEvents.size(); ++event) {
				if (!(previousEvents[event] > 0.0 && newEvents[event] <= 0.0))
					continue;
				double low = 0.0, high = 1.0;
				for (unsigned iteration = 0; iteration < 40; ++iteration) {
					const double mid = 0.5 * (low + high);
					const auto interpolated = interpolate(mid, step, y, stages[0], next, stages[6]);
					double value;
					if (event == 2) {
						value = interpolated[6] - m_physics.minMassKg;
					} else if (event == 0) {
						value = m_environment.clearance(geodetic({interpolated[0], interpolated[1], interpolated[2]}));
					} else {
						const auto eval = evaluate(interpolated, parameters, fragment.shape, m_environment, m_physics);
						++stats.evaluations;
						value = eventValues(interpolated, eval, fragment.strengthPa, m_environment, m_physics)[event];
					}
					if (value <= 0.0)
						high = mid;
					else
						low = mid;
				}
				if (high < earliest) {
					earliest = high;
					eventIndex = event;
				}
			}
			if (earliest <= 1.0) {
				y = interpolate(earliest, step, y, stages[0], next, stages[6]);
				time += earliest * step;
				constexpr FlightEvent events[]{FlightEvent::Ground, FlightEvent::Breakup, FlightEvent::MassCutoff};
				return finish(events[eventIndex]);
			}
			time += step;
			y = next;
			if (observer)
				observer(time, unpack(y));
			current = candidate;
			previousEvents = newEvents;
			step = std::min(step * std::min(5.0, 0.9 * std::pow(std::max(error, 1e-10), -0.2)), m_options.maxStepS);
		}
		throw std::runtime_error("Integration attempt limit exceeded");
	}

}  // namespace strewn
