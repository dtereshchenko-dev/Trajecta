#include "McObservationSampler.h"
#include "fit/FitObservations.h"
#include "fit/OrbitFitModel.h"
#include "hazard/GaussianSampler.h"
#include "orbital_core/Types.h"
#include "logger/Logger.h"
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace orbital_core::hazard::mc {

	ObservationSampler::ObservationSampler(const FitPipelineResult& source, std::stop_token cancel)
	    : source_(source),
	      included_(source.observations.all.size(), source.fit.keptObsIndices.empty() ? 1 : 0),
	      factors_(included_.size()),
	      cosDec_(included_.size())
	{
		if (included_.empty())
			throw std::invalid_argument("Observation MC requires observations");
		for (const auto index : source.fit.keptObsIndices) {
			if (index >= included_.size() || included_[index])
				throw std::invalid_argument("Invalid accepted observation mask");
			included_[index] = 1;
		}
		OrbitOptimizerConfig config;
		config.nbody = source.effectiveConfig.nbody;
		config.nonGrav = source.effectiveConfig.nonGrav;
		config.useWeights = source.effectiveConfig.useWeights;
		config.weightProfile = source.effectiveConfig.weightProfile;
		config.obsWeighting = source.effectiveConfig.obsWeighting;
		config.logDiagnostics = false;
		config.logCacheStats = false;
		auto prepared = fit::prepareFitObservations(source.observations.all, config);
		if (cancel.stop_requested())
			return;
		fit::OrbitFitModel model(source.fit.fitEpochOrbit, std::move(prepared), config);
		model.propagator().nonGrav() = source.fit.finalNonGrav;
		auto params = fit::OrbitFitModel::stateVectorFromOrbit(*source.fit.fitEpochOrbit);
		std::vector<double> steps(params.size(), 0.0);
		model.nonGravProfile().extendParams(params, steps);
		const auto evaluation = model.computeResidualsAndStateJac(params);
		if (cancel.stop_requested())
			return;
		if (evaluation.transforms.size() != included_.size())
			throw std::runtime_error("Observation MC measurement covariance is unavailable");
		for (std::size_t i = 0; i < included_.size(); ++i) {
			if (!included_[i])
				continue;
			const auto& observation = model.observations()[i];
			if (observation.radar) {
				const auto validate = [](const auto& channel) {
					if (channel && (!(channel->sigma > 0.0) || !std::isfinite(channel->sigma)))
						throw std::invalid_argument("Invalid radar measurement uncertainty");
				};
				validate(observation.radar->delaySeconds);
				validate(observation.radar->dopplerHz);
				continue;
			}
			const auto& t = evaluation.transforms[i];
			Eigen::Matrix2d whitening;
			whitening << t.m00, t.m01, t.m10, t.m11;
			const Eigen::FullPivLU<Eigen::Matrix2d> factor(whitening);
			if (!t.valid || !whitening.allFinite() || !factor.isInvertible())
				throw std::invalid_argument("Singular optical measurement covariance");
			factors_[i] = factor.inverse();
			cosDec_[i] = observation.raMetricCosDec;
			if (!std::isfinite(cosDec_[i]) || cosDec_[i] == 0)
				throw std::invalid_argument("Singular RA measurement metric");
		}
	}

	FitResult ObservationSampler::refit(int index, std::uint64_t seed, std::stop_token cancel) const
	{
		if (cancel.stop_requested()) {
			FitResult result;
			result.cancelled = true;
			return result;
		}
		auto observations = source_.observations.all;
		Eigen::VectorXd normal(static_cast<Eigen::Index>(observations.size() * 2));
		sampling::fillStandardNormal(normal, seed, index);
		for (std::size_t i = 0; i < observations.size(); ++i) {
			if (!included_[i])
				continue;
			auto& observation = observations[i];
			const double first = normal(2 * i), second = normal(2 * i + 1);
			if (observation.radar) {
				auto& radar = *observation.radar;
				if (radar.delaySeconds)
					radar.delaySeconds->value += first * radar.delaySeconds->sigma;
				if (radar.dopplerHz)
					radar.dopplerHz->value += second * radar.dopplerHz->sigma;
				continue;
			}
			const Eigen::Vector2d noise = factors_[i] * Eigen::Vector2d(first, second);
			observation.RA += noise(0) / cosDec_[i];
			observation.DEC += noise(1);
			const Vec3D equatorial(std::cos(observation.DEC) * std::cos(observation.RA),
			                       std::cos(observation.DEC) * std::sin(observation.RA), std::sin(observation.DEC));
			observation.RA = std::atan2(equatorial.y, equatorial.x);
			if (observation.RA < 0)
				observation.RA += 2 * Constants::PI;
			observation.DEC = std::asin(std::clamp(equatorial.z, -1.0, 1.0));
			observation.vector = {
			    equatorial.x,
			    equatorial.y * std::cos(Constants::EARTH_INCL) + equatorial.z * std::sin(Constants::EARTH_INCL),
			    -equatorial.y * std::sin(Constants::EARTH_INCL) + equatorial.z * std::cos(Constants::EARTH_INCL)};
		}
		ScopedLogSilence silence(LogLevel::Warning);
		return FitSession(source_.effectiveConfig)
		    .optimize({.initialOrbit = std::make_shared<Orbit>(*source_.fit.fitEpochOrbit),
		               .observations = observations,
		               .initialIncluded = included_,
		               .warmStartNonGrav = &source_.fit.finalNonGrav,
		               .options = {.selection = FitSelectionMode::Fixed, .uncertainty = FitUncertaintyMode::Skip}},
		              {.cancel = cancel});
	}

}  // namespace orbital_core::hazard::mc
