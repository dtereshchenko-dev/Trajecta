#include "fit/FitObservations.h"

#include "logger/Logger.h"
#include "obs/ObservationPipeline.h"
#include "obs/OverObservationReweight.h"

#include <iomanip>
#include <cmath>
#include <stdexcept>

namespace orbital_core::fit {

	std::vector<obs::PreparedObservation> prepareFitObservations(std::span<const ObsNode> observations,
	                                                             const OrbitOptimizerConfig& config)
	{
		const obs::ObservationPipeline pipeline(config.weightProfile, config.obsWeighting);
		auto prepared = pipeline.prepare(observations, true);
		if (config.gaiaPhotocenterRadiusKm &&
		    (!std::isfinite(*config.gaiaPhotocenterRadiusKm) || *config.gaiaPhotocenterRadiusKm <= 0.0))
			throw std::invalid_argument("Gaia photocenter radius must be positive and finite");
		bool gaiaWithoutRadius = false;
		for (auto& observation : prepared) {
			if (observation.astrometricConvention == obs::AstrometricConvention::GaiaFpr &&
			    config.gaiaPhotocenterRadiusKm)
				observation.photocenterRadiusAu = *config.gaiaPhotocenterRadiusKm / (Constants::AU_IN_METERS / 1000.0);
			else if (observation.astrometricConvention == obs::AstrometricConvention::GaiaFpr)
				gaiaWithoutRadius = true;
		}
		if (gaiaWithoutRadius && config.logDiagnostics)
			LOG_WARN << "[Gaia FPR] object radius unavailable; photocenter correction not applied";
		else if (config.gaiaPhotocenterRadiusKm && config.logDiagnostics)
			LOG_INFO << "[Gaia FPR] photocenter radius=" << *config.gaiaPhotocenterRadiusKm << " km";
		relaxSameStationBatches(prepared, config);
		return prepared;
	}

	void relaxSameStationBatches(std::vector<obs::PreparedObservation>& observations,
	                             const OrbitOptimizerConfig& config)
	{
		const auto factors = config.useWeights
		                         ? obs::computeSameStationOverObservationFactors(observations, config.weightProfile)
		                         : std::vector<obs::SameStationReweightFactor>{};
		if (factors.size() == observations.size()) {
			for (std::size_t i = 0; i < observations.size(); ++i) {
				const auto& f = factors[i];
				if (f.ra > 0.0)
					observations[i].sigmaRa /= f.ra;
				if (f.dec > 0.0)
					observations[i].sigmaDec /= f.dec;
			}
		}
		if (!config.logDiagnostics)
			return;
		if (factors.empty()) {
			LOG_INFO << "[Overobs] inactive"
			         << (config.useWeights ? ": no same-station nearby downweighting"
			                               : ": observation weighting disabled");
			return;
		}
		const auto stats = obs::summarizeSameStationOverObservationFactors(factors);
		LOG_INFO << "[Overobs] vfcc17 same-station batch gap=8h ceiling=4"
		         << " downweighted=" << stats.downweightedCount << "/" << stats.observationCount
		         << " maxNearby=" << stats.maxEffectiveNearby << " factor[min/mean/max]=" << std::fixed
		         << std::setprecision(6) << stats.minFactor << "/" << stats.meanFactor << "/" << stats.maxFactor;
	}

	std::vector<obs::PreparedObservation> selectObservations(std::span<const obs::PreparedObservation> observations,
	                                                         std::span<const std::uint8_t> included)
	{
		if (included.size() != observations.size())
			throw std::invalid_argument("selectObservations: mask size mismatch");
		std::vector<obs::PreparedObservation> selected;
		selected.reserve(observations.size());
		for (std::size_t i = 0; i < observations.size(); ++i) {
			if (included[i])
				selected.push_back(observations[i]);
		}
		return selected;
	}

}  // namespace orbital_core::fit
