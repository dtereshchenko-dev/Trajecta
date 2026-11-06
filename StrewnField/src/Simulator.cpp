#include "strewn_core/Simulator.h"
#include "strewn_core/AtmosphericInput.h"
#include "strewn_core/FragmentationModel.h"
#include "VectorMath.h"
#include "Random.h"
#include "Convergence.h"
#include "RunDeadline.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <optional>
#include <queue>
#include <stdexcept>
#include <thread>
#include <type_traits>

namespace strewn {
	namespace {

		void validateEntry(const EntryState& entry)
		{
			if (!detail::finite(entry.positionEcefM) || !detail::finite(entry.velocityEcefMps) ||
			    !std::isfinite(entry.epochMjdUtc) || !std::isfinite(entry.weight) || entry.weight < 0.0)
				throw std::invalid_argument("Invalid entry state, epoch or weight");
			geodetic(entry.positionEcefM);
		}

		void validateRanges(const ParameterRanges& ranges)
		{
			const Range all[] = {ranges.densityKgM3,  ranges.strongStrengthPa, ranges.weakStrengthPa,
			                     ranges.dustFraction, ranges.largestFraction,  ranges.beta,
			                     ranges.windScale,    ranges.dragScale,        ranges.ablationKgJ};
			for (auto range : all)
				if (!std::isfinite(range.min) || !std::isfinite(range.max) || range.min > range.max)
					throw std::invalid_argument("Invalid physical parameter range");
			for (bool upper : {false, true}) {
				const auto value = [upper](Range r) { return upper ? r.max : r.min; };
				RunParameters p;
				p.densityKgM3 = value(ranges.densityKgM3);
				p.strengthPa = value(ranges.strongStrengthPa);
				p.dustFraction = value(ranges.dustFraction);
				p.largestFraction = value(ranges.largestFraction);
				p.beta = value(ranges.beta);
				p.windScale = value(ranges.windScale);
				p.dragScale = value(ranges.dragScale);
				p.ablationKgJ = value(ranges.ablationKgJ);
				validate(p);
				p.strengthPa = value(ranges.weakStrengthPa);
				validate(p);
			}
			if (!std::isfinite(ranges.weakProbability) || ranges.weakProbability < 0 || ranges.weakProbability > 1)
				throw std::invalid_argument("weakProbability must be in [0, 1]");
			std::visit(
			    [](const auto& size) {
				    using T = std::decay_t<decltype(size)>;
				    if constexpr (std::is_same_v<T, PhotometricSize>) {
					    if (!std::isfinite(size.absoluteMagnitude) || !std::isfinite(size.albedoMin) ||
					        !std::isfinite(size.albedoMode) || !std::isfinite(size.albedoMax) || size.albedoMin <= 0 ||
					        size.albedoMax > 1 || size.albedoMin > size.albedoMode || size.albedoMode > size.albedoMax)
						    throw std::invalid_argument("Invalid photometric size");
				    } else {
					    Range range;
					    if constexpr (std::is_same_v<T, MassSize>)
						    range = size.kg;
					    else
						    range = size.m;
					    if (!std::isfinite(range.min) || !std::isfinite(range.max) || range.min <= 0 ||
					        range.max < range.min)
						    throw std::invalid_argument("Invalid mass or diameter range");
				    }
			    },
			    ranges.size);
		}

		RunParameters sample(const ParameterRanges& ranges, detail::Random& random)
		{
			RunParameters result;
			result.densityKgM3 = random.uniform(ranges.densityKgM3);
			result.massKg = std::visit(
			    [&](const auto& size) {
				    using T = std::decay_t<decltype(size)>;
				    if constexpr (std::is_same_v<T, MassSize>) {
					    return random.uniform(size.kg);
				    } else {
					    double diameter;
					    if constexpr (std::is_same_v<T, DiameterSize>) {
						    diameter = random.uniform(size.m);
					    } else {
						    const double u = random.uniform();
						    double albedo = size.albedoMin;
						    if (size.albedoMax > size.albedoMin) {
							    const double width = size.albedoMax - size.albedoMin;
							    albedo =
							        u < (size.albedoMode - size.albedoMin) / width
							            ? size.albedoMin + std::sqrt(u * width * (size.albedoMode - size.albedoMin))
							            : size.albedoMax -
							                  std::sqrt((1.0 - u) * width * (size.albedoMax - size.albedoMode));
						    }
						    diameter = 1329e3 * std::pow(10.0, -size.absoluteMagnitude / 5.0) / std::sqrt(albedo);
					    }
					    return detail::pi * result.densityKgM3 * diameter * diameter * diameter / 6.0;
				    }
			    },
			    ranges.size);
			const bool weak = random.uniform() < ranges.weakProbability;
			result.strengthPa = random.uniform(weak ? ranges.weakStrengthPa : ranges.strongStrengthPa);
			result.dustFraction = random.uniform(ranges.dustFraction);
			result.largestFraction = random.uniform(ranges.largestFraction);
			result.beta = random.uniform(ranges.beta);
			result.windScale = random.uniform(ranges.windScale);
			result.dragScale = random.uniform(ranges.dragScale);
			result.ablationKgJ = random.uniform(ranges.ablationKgJ);
			validate(result);
			return result;
		}

		struct LargerFirst {
			bool operator()(const Fragment& a, const Fragment& b) const
			{
				return a.state.massKg != b.state.massKg ? a.state.massKg < b.state.massKg : a.id > b.id;
			}
		};

	}  // namespace

	void validate(const MonteCarloConfig& options)
	{
		detail::validateRunTimeout(options.runTimeoutS);
		if (options.runs == 0 || options.workers == 0 || options.workers > 1024)
			throw std::invalid_argument("Runs and workers must be nonzero; workers <= 1024");
		validateRanges(options.parameters);
		if (options.convergence)
			detail::validateConvergence(*options.convergence);
	}

	Simulator::Simulator(const Environment& environment, PhysicsConfig physics, IntegrationOptions integration)
	    : m_environment(environment), m_physics(std::move(physics)), m_integration(std::move(integration))
	{
		validate(m_physics, m_integration);
	}

	RunResult Simulator::runOne(const EntryState& entry, const RunParameters& parameters, std::uint64_t seed,
	                            std::uint64_t runId, std::stop_token cancel, double runTimeoutS) const
	{
		return runOne(entry, parameters, seed, runId, cancel, FragmentObserver{}, runTimeoutS);
	}

	RunResult Simulator::runOne(const EntryState& entry, const RunParameters& parameters, std::uint64_t seed,
	                            std::uint64_t runId, std::stop_token cancel, const FragmentObserver& observer,
	                            double runTimeoutS) const
	{
		const auto deadline = detail::runDeadline(runTimeoutS);
		validateEntry(entry);
		validate(parameters);
		RunResult result;
		result.runId = runId;
		result.entry = entry;
		result.parameters = parameters;
		result.mass.initialKg = parameters.massKg;
		detail::Random random(seed, runId, 2);
		TrajectoryIntegrator integrator(m_environment, m_physics, m_integration);
		FragmentationModel fragmentation(m_physics);
		std::priority_queue<Fragment, std::vector<Fragment>, LargerFirst> queue;
		Fragment root;
		root.state = {entry.positionEcefM, entry.velocityEcefMps, parameters.massKg};
		root.strengthPa = parameters.strengthPa;
		root.shape = m_physics.phase == FlightPhase::Entry ? Shape::Sphere : m_physics.fragmentShape;
		queue.push(root);
		std::uint64_t nextId = 1;
		while (!queue.empty()) {
			if (const auto interrupted = detail::interruption(cancel, deadline)) {
				const detail::RunInterrupted reason(*interrupted);
				result.status = reason.status();
				result.error = reason.what();
				break;
			}
			if (m_physics.maxParticles && result.particles >= m_physics.maxParticles) {
				result.status = RunStatus::ParticleLimit;
				break;
			}
			auto fragment = queue.top();
			queue.pop();
			++result.particles;
			try {
				FlightObserver trace;
				if (observer)
					trace = [&](double time, const FlightState& state) { observer(fragment, time, state); };
				const auto flight = integrator.advance(fragment, parameters, cancel, trace, deadline);
				result.integration.evaluations += flight.stats.evaluations;
				result.integration.acceptedSteps += flight.stats.acceptedSteps;
				result.integration.rejectedSteps += flight.stats.rejectedSteps;
				result.mass.ablatedKg +=
				    std::max(0.0, fragment.state.massKg - flight.state.massKg) * fragment.multiplicity;
				fragment.state = flight.state;
				fragment.timeS = flight.timeS;
				const auto position = geodetic(fragment.state.positionEcefM);
				if (flight.event == FlightEvent::Ground) {
					result.impacts.push_back({fragment, position.latitudeRad / detail::deg,
					                          position.longitudeRad / detail::deg, position.heightWgs84M});
					result.mass.landedKg += fragment.state.massKg * fragment.multiplicity;
				} else if (flight.event == FlightEvent::Breakup) {
					const auto groups = fragmentation.groups(fragment, parameters, cancel, deadline);
					std::vector<double> draws(groups.size());
					for (std::size_t i = 0; i < draws.size(); ++i) {
						if ((i & 255) == 0)
							detail::checkInterruption(cancel, deadline);
						draws[i] = random.uniform();
					}
					auto batch = fragmentation.split(fragment, parameters, m_environment.air(position), groups, draws,
					                                 cancel, deadline);
					result.breakups.push_back({fragment, position.heightWgs84M, batch.children.size()});
					result.mass.dustKg += batch.dustKg;
					result.mass.residualKg += batch.residualKg;
					for (auto& child : batch.children) {
						child.id = nextId++;
						queue.push(std::move(child));
					}
				} else if (flight.event == FlightEvent::MassCutoff) {
					result.mass.belowCutoffKg += fragment.state.massKg * fragment.multiplicity;
				} else {
					result.mass.unresolvedKg += fragment.state.massKg * fragment.multiplicity;
					if (flight.event == FlightEvent::Timeout) {
						result.status = RunStatus::Timeout;
						result.error = "Simulated flight-time limit exceeded";
					} else {
						const detail::RunInterrupted reason(flight.event);
						result.status = reason.status();
						result.error = reason.what();
					}
					break;
				}
			} catch (const detail::RunInterrupted& error) {
				result.status = error.status();
				result.error = error.what();
				result.mass.unresolvedKg += fragment.state.massKg * fragment.multiplicity;
				break;
			} catch (const std::exception& error) {
				result.status = RunStatus::Failed;
				result.error = error.what();
				result.mass.unresolvedKg += fragment.state.massKg * fragment.multiplicity;
				break;
			}
		}
		while (!queue.empty()) {
			result.mass.unresolvedKg += queue.top().state.massKg * queue.top().multiplicity;
			queue.pop();
		}
		return result;
	}

	RunSummary Simulator::run(std::span<const EntryState> entries, const MonteCarloConfig& options,
	                          const RunControl& control, const RunConsumer& consume) const
	{
		detail::validateRunTimeout(options.runTimeoutS);
		if (entries.empty() || options.runs == 0 || options.workers == 0 || options.workers > 1024)
			throw std::invalid_argument("Entries, runs and workers must be nonzero; workers <= 1024");
		validateRanges(options.parameters);
		std::vector<double> cumulative;
		cumulative.reserve(entries.size());
		double totalWeight = 0.0;
		for (const auto& entry : entries) {
			validateEntry(entry);
			totalWeight += entry.weight;
			cumulative.push_back(totalWeight);
		}
		if (!std::isfinite(totalWeight) || totalWeight <= 0.0)
			throw std::invalid_argument("Entry weights must have a positive finite sum");
		RunSummary summary;
		summary.requested = options.runs;
		summary.workers = static_cast<unsigned>(std::min(options.runs, static_cast<std::size_t>(options.workers)));
		if (control.cancel.stop_requested()) {
			summary.cancelled = true;
			return summary;
		}
		const auto start = std::chrono::steady_clock::now();
		std::optional<detail::ConvergenceMonitor> convergence;
		if (options.convergence)
			convergence.emplace(*options.convergence);
		const std::size_t window = 2 * summary.workers;
		std::vector<std::optional<RunResult>> slots(window);
		std::mutex mutex;
		std::condition_variable changed;
		std::size_t next = 0, committed = 0;
		std::stop_source stop;
		std::exception_ptr workerException;
		std::stop_callback propagateStop(control.cancel, [&] {
			stop.request_stop();
			changed.notify_all();
		});
		std::vector<std::jthread> workers;
		workers.reserve(summary.workers);
		const auto worker = [&] {
			try {
				for (;;) {
					std::size_t id;
					{
						std::unique_lock lock(mutex);
						changed.wait(lock, [&] {
							return stop.stop_requested() || next >= options.runs || next < committed + window;
						});
						if (stop.stop_requested() || next >= options.runs)
							return;
						id = next++;
					}
					detail::Random entryRandom(options.seed, id, 0), parameterRandom(options.seed, id, 1);
					const double choice =
					    std::min(entryRandom.uniform() * totalWeight, std::nextafter(totalWeight, 0.0));
					const auto index = static_cast<std::size_t>(
					    std::upper_bound(cumulative.begin(), cumulative.end(), choice) - cumulative.begin());
					auto result = runOne(entries[index], sample(options.parameters, parameterRandom), options.seed, id,
					                     stop.get_token(), options.runTimeoutS);
					{
						std::lock_guard lock(mutex);
						slots[id % window] = std::move(result);
					}
					changed.notify_all();
				}
			} catch (...) {
				{
					std::lock_guard lock(mutex);
					if (!workerException)
						workerException = std::current_exception();
				}
				stop.request_stop();
				changed.notify_all();
			}
		};
		const auto join = [&] {
			stop.request_stop();
			changed.notify_all();
			workers.clear();
		};
		try {
			for (unsigned i = 0; i < summary.workers; ++i)
				workers.emplace_back(worker);
			for (std::size_t id = 0; id < options.runs; ++id) {
				RunResult result;
				{
					std::unique_lock lock(mutex);
					changed.wait(lock, [&] { return slots[id % window].has_value() || stop.stop_requested(); });
					if (workerException)
						std::rethrow_exception(workerException);
					if (control.cancel.stop_requested() || !slots[id % window]) {
						summary.cancelled = true;
						break;
					}
					result = std::move(*slots[id % window]);
					slots[id % window].reset();
					committed = id + 1;
				}
				changed.notify_all();
				if (result.status == RunStatus::Cancelled) {
					summary.cancelled = true;
					break;
				}
				if (consume)
					consume(result);
				++summary.completed;
				summary.truncated += result.status == RunStatus::ParticleLimit;
				summary.failed += result.status == RunStatus::Failed || result.status == RunStatus::Timeout ||
				                  result.status == RunStatus::WallTimeout;
				summary.wallTimeouts += result.status == RunStatus::WallTimeout;
				summary.runsWithImpacts += !result.impacts.empty();
				summary.trackedImpacts += result.impacts.size();
				if (convergence)
					convergence->update(result, summary);
				summary.elapsedS = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
				if (control.progress)
					control.progress(summary);
				if (control.cancel.stop_requested()) {
					summary.cancelled = true;
					summary.converged = false;
					break;
				}
				if (summary.converged)
					break;
			}
		} catch (...) {
			join();
			throw;
		}
		join();
		summary.elapsedS = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
		return summary;
	}

}  // namespace strewn
