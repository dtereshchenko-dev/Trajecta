#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include "orbital_core/SpkReader.h"
#include "McBoundaryScan.h"
#include "McObservationSampler.h"
#include "GaussianSampler.h"
#include "orbital_core/hazard/PlanetData.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <mutex>
#include <stdexcept>

namespace orbital_core::hazard {
	namespace {
		void initialize(McClone& clone, const Orbit& orbit, const NonGravModel& ng)
		{
			clone.epochMjdTdb = orbit.epochMjd();
			const auto r = orbit.getR(), v = orbit.getV() * Constants::GAUSS_K;
			clone.state = {r.x, r.y, r.z, v.x, v.y, v.z};
			clone.nonGrav = ng;
		}

		Eigen::MatrixXd covarianceFactor(const FitParameterCovariance& covariance, const NonGravModel& nonGrav)
		{
			const int dimension = 6 + nonGrav.fitParameterCount();
			if (covariance.rows() != dimension || covariance.cols() != dimension || !covariance.allFinite() ||
			    (covariance.diagonal().array() <= 0).any())
				throw std::invalid_argument("Full fit-epoch covariance is unavailable");
			const Eigen::VectorXd scale = covariance.diagonal().cwiseSqrt();
			const Eigen::MatrixXd correlation =
			    scale.cwiseInverse().asDiagonal() * covariance * scale.cwiseInverse().asDiagonal();
			if (!correlation.allFinite() || !correlation.isApprox(correlation.transpose(), 1e-10))
				throw std::invalid_argument("Fit-epoch covariance is not symmetric");
			const Eigen::LLT<Eigen::MatrixXd> decomposition(correlation);
			if (decomposition.info() != Eigen::Success)
				throw std::invalid_argument("Fit-epoch covariance is not positive definite");
			return scale.asDiagonal() * Eigen::MatrixXd(decomposition.matrixL());
		}
	}  // namespace

	static IntervalMcResult runIntervalMonteCarloImpl(const Orbit& orbit, const NonGravModel& nonGrav,
	                                                  const FitParameterCovariance& covariance,
	                                                  const NBodySettings& forces, const FitPipelineResult* source,
	                                                  const IntervalMcConfig& config, std::stop_token cancel,
	                                                  const std::function<void(int, int)>& progress)
	{
		if (!forces.enabled)
			throw std::invalid_argument("MC requires an N-body accepted solution");
		if (config.method != McMethod::Covariance && config.method != McMethod::Observations)
			throw std::invalid_argument("Unknown MC method");
		if (config.method == McMethod::Observations && !source)
			throw std::invalid_argument("Observation MC requires observations and a fitted solution");
		if (config.integrator != IntegratorKind::RA15 && config.integrator != IntegratorKind::IAS15)
			throw std::invalid_argument("Unsupported MC integrator");
		if (config.target != Planet::Earth && config.target != Planet::Moon && config.target != Planet::Mars)
			throw std::invalid_argument("Unsupported MC target");
		if (!std::isfinite(config.startMjdTdb) || !std::isfinite(config.endMjdTdb) ||
		    !std::isfinite(orbit.epochMjd()) || config.startMjdTdb < orbit.epochMjd() ||
		    config.endMjdTdb <= config.startMjdTdb || !std::isfinite(config.altitudeKm) || config.altitudeKm < 0 ||
		    config.samples < 1 || config.samples > 1000000 || !config.workers || config.workers > 1024)
			throw std::invalid_argument("Invalid MC interval or sampling settings");
		IntervalMcResult result;
		result.config = config;
		initialize(result.nominal, orbit, nonGrav);
		result.clones.resize(config.samples);
		for (int i = 0; i < config.samples; ++i)
			result.clones[i].index = i;
		if (cancel.stop_requested()) {
			result.cancelled = true;
			mc::summarize(result);
			return result;
		}
		Eigen::MatrixXd factor;
		if (config.method == McMethod::Covariance)
			factor = covarianceFactor(covariance, nonGrav);
		auto* eph = getDefaultEphemerisSource();
		if (!eph || !eph->isAvailable())
			throw std::runtime_error("Planetary ephemerides are unavailable");
		const auto rotation = bodyFixedTransformMjd(config.target, config.startMjdTdb, eph);
		if (!rotation.valid || !bodyFixedTransformMjd(config.target, config.endMjdTdb, eph).valid)
			throw std::runtime_error("Body-fixed orientation is unavailable for the requested interval");
		result.bodyFixedFrame = rotation.frame;
		result.bodyFixedModel =
		    config.target == Planet::Earth ? earthBodyFixedFrameModel() : std::string(rotation.model);
		if (!std::isfinite(forces.stepSize) || forces.stepSize <= 0)
			throw std::invalid_argument("Invalid accepted propagation step");
		if (forces.useSb441 && (!getDefaultAsteroidSource() || !getDefaultAsteroidSource()->isAvailable()))
			throw std::runtime_error("SB441 is unavailable");
		NBodyPropagator base(config.integrator);
		base.setStepSize(forces.stepSize);
		base.setPlanetMask(forces.planetMask);
		base.setUseSb441(forces.useSb441);
		base.setTargetNaifId(forces.targetNaifId);
		base.setInactiveInteriorMassPolicy(forces.inactiveInteriorMassPolicy);
		base = base.forEncounterMonitoring(*eph);
		const auto propagate = [&](McClone& clone, std::stop_token stop) {
			try {
				if (!std::all_of(clone.state.begin(), clone.state.end(), [](double x) { return std::isfinite(x); }) ||
				    !std::isfinite(clone.nonGrav.A1()) || !std::isfinite(clone.nonGrav.A2()) ||
				    !std::isfinite(clone.nonGrav.A3()) || !std::isfinite(clone.nonGrav.DT()))
					throw std::runtime_error("Non-finite MC initial state");
				if (clone.epochMjdTdb > config.startMjdTdb)
					throw std::runtime_error("Refit epoch is later than the scan start");
				auto propagator = base;
				propagator.nonGrav() = clone.nonGrav;
				mc::scan(clone, config, propagator, *eph, stop);
				clone.status = McCloneStatus::Complete;
			} catch (const std::exception& e) {
				clone.status = stop.stop_requested() ? McCloneStatus::Cancelled : McCloneStatus::PropagationFailed;
				clone.error = e.what();
				clone.entry.reset();
			}
		};
		std::unique_ptr<mc::ObservationSampler> sampler;
		if (!cancel.stop_requested() && config.method == McMethod::Observations)
			sampler = std::make_unique<mc::ObservationSampler>(*source, cancel);
		propagate(result.nominal, cancel);
		std::stop_source stop;
		std::stop_callback external(cancel, [&] { stop.request_stop(); });
		std::mutex progressMutex;
		int done = 0;
		std::exception_ptr callbackError;
#pragma omp parallel for schedule(dynamic, 1) num_threads(config.workers)
		for (int i = 0; i < config.samples; ++i) {
			if (stop.stop_requested())
				continue;
			auto& clone = result.clones[i];
			try {
				if (config.method == McMethod::Covariance) {
					initialize(clone, orbit, nonGrav);
					Eigen::VectorXd normal(factor.cols());
					sampling::fillStandardNormal(normal, config.seed, i);
					const Eigen::VectorXd offset = factor * normal;
					for (int j = 0; j < 6; ++j)
						clone.state[j] += offset(j);
					for (int j = 0; j < clone.nonGrav.fitParameterCount(); ++j)
						clone.nonGrav.setFitParameterValue(j, clone.nonGrav.fitParameterValue(j) + offset(6 + j));
				} else {
					const auto fit = sampler->refit(i, config.seed, stop.get_token());
					if (fit.cancelled || stop.stop_requested())
						throw std::runtime_error("MC cancelled");
					if (!fit.converged || !fit.fitEpochOrbit)
						throw std::runtime_error("Observation MC fit did not converge");
					initialize(clone, *fit.fitEpochOrbit, fit.finalNonGrav);
				}
				propagate(clone, stop.get_token());
			} catch (const std::exception& e) {
				clone.status = stop.stop_requested()                     ? McCloneStatus::Cancelled
				               : config.method == McMethod::Observations ? McCloneStatus::FitFailed
				                                                         : McCloneStatus::PropagationFailed;
				clone.error = e.what();
			}
			std::lock_guard lock(progressMutex);
			++done;
			if (progress && !callbackError) {
				try {
					progress(done, config.samples);
				} catch (...) {
					callbackError = std::current_exception();
					stop.request_stop();
				}
			}
		}
		if (callbackError)
			std::rethrow_exception(callbackError);
		result.cancelled = cancel.stop_requested();
		mc::summarize(result);
		return result;
	}

	IntervalMcResult runIntervalMonteCarlo(const FitPipelineResult& source, const IntervalMcConfig& config,
	                                       std::stop_token cancel, std::function<void(int, int)> progress)
	{
		if (!source.fit.fitEpochOrbit || source.fit.cancelled || !source.fit.converged)
			throw std::invalid_argument("MC requires a converged accepted solution");
		return runIntervalMonteCarloImpl(*source.fit.fitEpochOrbit, source.fit.finalNonGrav,
		                                 source.fit.uncertainty.fitParameterCovariance, source.effectiveConfig.nbody,
		                                 &source, config, cancel, progress);
	}

	IntervalMcResult runIntervalMonteCarlo(const OrbitSolution& solution, const NBodySettings& forces,
	                                       const IntervalMcConfig& config, std::stop_token cancel,
	                                       std::function<void(int, int)> progress)
	{
		return runIntervalMonteCarloImpl(solution.nominal, solution.nonGrav, solution.covariance, forces, nullptr,
		                                 config, cancel, progress);
	}
}  // namespace orbital_core::hazard
