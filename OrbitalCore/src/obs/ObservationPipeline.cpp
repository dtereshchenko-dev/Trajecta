#include "ObservationPipeline.h"
#include <span>
#include <string>
#include <vector>
#include "AstrometricSigma.h"
#include "ObserverGeometry.h"

#include "EphemeridesContainer.h"
#include "ObservatoryCatalog.h"
#include "logger/Logger.h"
#include "orbital_core/Constants.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/Types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace orbital_core::obs {
	namespace {
		constexpr double kDefaultSurfaceRadiusKm = 1.0;

		[[nodiscard]] const char* weightProfileName(AstrometricWeightProfile profile) noexcept
		{
			switch (profile) {
				case AstrometricWeightProfile::Vfcc17:
					return "vfcc17.rules";
				case AstrometricWeightProfile::Uniform:
					return "uniform 0.5 arcsec";
			}
			return "unknown";
		}

	}  // namespace

	ObservationPipeline::ObservationPipeline(AstrometricWeightProfile weightProfile,
	                                         ObservationWeightingOptions options)
	    : m_weightProfile(weightProfile), m_options(options)
	{
		LOG_INFO << "ObservationPipeline: ensure sigma table";
		ensureSigmaTable(m_weightProfile);
		LOG_INFO << "ObservationPipeline: ready";
	}

	std::vector<PreparedObservation> ObservationPipeline::prepare(std::span<const ObsNode> rawObs,
	                                                              bool assumeSortedInput) const
	{
		std::vector<PreparedObservation> out;
		out.reserve(rawObs.size());

		EphemeridesContainer eph(EphemeridesContainer::Planet::Earth);

		std::size_t satCount = 0;
		std::size_t groundCount = 0;
		std::size_t defaultSurfaceRadii = 0;
		std::size_t sigmaResolved = 0;
		std::size_t observerSigmaUsed = 0;
		std::size_t observerCorrUsed = 0;
		std::size_t timingFromRms = 0;
		std::size_t timingFromPrecision = 0;
		std::size_t timingFromDefault = 0;
		std::size_t timingUnstated = 0;
		std::size_t timingVideo = 0;
		std::size_t timingExcluded = 0;

		const ObservatoryCatalog& catalog = observatoryCatalog();
		for (const auto& obs : rawObs) {
			const auto station = catalog.find(obs.obsCode);
			if (!obs.radar && station == catalog.end())
				throw std::runtime_error("Observer code not found: " + obs.obsCode);
			const auto& receiver = obs.radar ? radarObservatory(obs.obsCode) : station->second;
			const ObserverData observer =
			    obs.rovingObserver ? rovingObserverData(*obs.rovingObserver, obs.obsCode)
			                       : ObserverData{receiver.longitude, receiver.cosPhi, receiver.sinPhi, receiver.code};

			double jd_utc = obs.JD;
			double mjd_utc = jd_utc - 2400000.5;
			const double mjd_tdb = obs.epochMjdTdb();
			if (obs.hasSplitUtc) {
				jd_utc = 2400000.5 + (static_cast<double>(obs.utcMjdDay) + obs.utcSecondsOfDay / 86400.0);
				mjd_utc = static_cast<double>(obs.utcMjdDay) + obs.utcSecondsOfDay / 86400.0;
			}

			const double jd_tdb = 2400000.5 + mjd_tdb;
			const auto& earthData = eph.getMjd(mjd_tdb);

			out.emplace_back();
			PreparedObservation& prep = out.back();
			prep.jdTdb = jd_tdb;
			prep.mjdTdb = mjd_tdb;
			prep.jdUtc = jd_utc;
			prep.mjdUtc = mjd_utc;
			prep.raObs = obs.RA;
			prep.decObs = obs.DEC;
			prep.raMetricCosDec = std::cos(obs.submittedDecRad());
			prep.earthR = earthData.R;
			prep.earthV = earthData.V;
			prep.observer = observer;
			if (obs.isSatellite) {
				prep.isSatellite = true;
				prep.observerOffsetEclAu = obs.satelliteOffsetEclAu;
				++satCount;
			}

			else {
				const ObserverState observerState = groundObserverStateEclAuMjd(mjd_tdb, mjd_utc, observer);
				prep.observerOffsetEclAu = observerState.offsetEclAu;
				prep.observerVelocityEclAuPerDay = observerState.velocityEclAuPerDay;
				++groundCount;
			}

			prep.note1 = obs.note1;
			prep.note2 = obs.note2;
			prep.catalogCode = obs.catalogCode;
			prep.astrometricConvention = obs.astrometricConvention;
			if (obs.radar) {
				if (!obs.radar->valid())
					throw std::invalid_argument("Invalid radar observation");
				const auto& tx = radarObservatory(obs.radar->transmitter);
				prep.radar = obs.radar;
				prep.radar->transmitter = tx.code;
				if (prep.radar->returnType == RadarReturn::Surface && !prep.radar->surfaceRadiusKm) {
					prep.radar->surfaceRadiusKm = kDefaultSurfaceRadiusKm;
					++defaultSurfaceRadii;
				}
				prep.transmitter = {tx.longitude, tx.cosPhi, tx.sinPhi, tx.code};
				if (std::hypot(observer.cosPhi, observer.sinPhi) == 0.0 ||
				    std::hypot(prep.transmitter.cosPhi, prep.transmitter.sinPhi) == 0.0)
					throw std::invalid_argument("Radar requires terrestrial station coordinates");
				continue;
			}
			const ObservationUncertainty uncertainty = lookupUncertainty(obs, jd_utc, m_weightProfile);
			prep.sigmaRa = uncertainty.sigmaRaRad;
			prep.sigmaDec = uncertainty.sigmaDecRad;
			prep.forceWeightRa = uncertainty.forceWeightRa;
			prep.forceWeightDec = uncertainty.forceWeightDec;
			prep.mag = obs.MAG;
			prep.magBand = obs.magBand;
			prep.nuclearMagnitude = obs.nuclearMagnitude;
			prep.magSigma = photometry::resolveSigma(obs.magString, obs.rmsMag);
			if (uncertainty.resolved) {
				++sigmaResolved;
			}

			if (usesObserverSigma(
			        obs, m_options.useObserverSigmas || obs.astrometricConvention == AstrometricConvention::GaiaFpr)) {
				prep.sigmaRa = obs.rmsRaArcsec * Constants::ARCSEC_RAD;
				prep.sigmaDec = obs.rmsDecArcsec * Constants::ARCSEC_RAD;
				prep.observerSigma = true;
				++observerSigmaUsed;
				if (std::isfinite(obs.rmsCorr) && std::abs(obs.rmsCorr) < 1.0) {
					prep.sigmaCorr = obs.rmsCorr;
					if (std::abs(obs.rmsCorr) > 0.05)
						++observerCorrUsed;
				}
			}

			if (m_options.timingEnabled()) {
				const bool excluded = obs.isSatellite || obs.note2 == 'E' || obs.obsCode == "258";
				if (excluded) {
					prep.timeSigmaDays = 0.0;
					++timingExcluded;
				} else {
					double sigmaDays = 0.0;
					if (obs.hasRmsTime && obs.rmsTimeSec > 0.0) {
						sigmaDays = obs.rmsTimeSec / 86400.0;
						++timingFromRms;
					} else if (obs.hasTimeAccuracy && obs.timeAccuracyDays > 0.0) {
						sigmaDays = obs.timeAccuracyDays;
						++timingFromPrecision;
					} else if (m_options.timingSigmaDefaultSeconds > 0.0) {
						sigmaDays = m_options.timingSigmaDefaultSeconds / 86400.0;
						++timingFromDefault;
					} else {
						++timingUnstated;
					}

					if (obs.note2 == 'n') {
						sigmaDays = 0.010 / 86400.0;
						++timingVideo;
					}

					prep.timeSigmaDays = sigmaDays;
				}
			}
		}

		if (sigmaResolved > 0) {
			LOG_INFO << "applyWeights: " << sigmaResolved << "/" << out.size() << " obs got sigma from "
			         << weightProfileName(m_weightProfile);
		}

		if (m_options.timingEnabled()) {
			LOG_INFO << "applyWeights: timing sigma of " << out.size() << " obs -- rmsTime " << timingFromRms
			         << ", stated precision " << timingFromPrecision << ", default "
			         << m_options.timingSigmaDefaultSeconds << "s " << timingFromDefault << ", unstated "
			         << timingUnstated << ", video 10ms " << timingVideo << ", excluded (spacecraft/occ/Gaia) "
			         << timingExcluded;
		}

		if (m_options.useObserverSigmas || observerSigmaUsed > 0) {
			LOG_INFO << "applyWeights: observer sigma used for " << observerSigmaUsed << "/" << out.size()
			         << " obs; the rest kept " << weightProfileName(m_weightProfile);
			if (observerCorrUsed > 0) {
				LOG_INFO << "applyWeights: " << observerCorrUsed
				         << " obs carry |rmsCorr| > 0.05, applied as a tilt of"
				            " the error ellipse";
			}
		}

		if (defaultSurfaceRadii > 0) {
			LOG_WARN << "prepareObservations: " << defaultSurfaceRadii
			         << " surface radar returns have no body radius; using " << kDefaultSurfaceRadiusKm << " km";
		}

		if (satCount > 0) {
			LOG_INFO << "prepareObservations: total=" << out.size() << " ground=" << groundCount
			         << " satellite=" << satCount;
		}

		if (!assumeSortedInput) {
			std::sort(out.begin(), out.end(),
			          [](const PreparedObservation& a, const PreparedObservation& b) { return a.mjdTdb < b.mjdTdb; });
		}

		return out;
	}
}  // namespace orbital_core::obs
