#include "orbital_core/ObservationSet.h"

#include "AstrometryBias.h"
#include "FileReader.h"
#include "logger/Logger.h"
#include "obs/IObservationReader.h"
#include "obs/ObservatoryCatalog.h"
#include "orbital_core/Types.h"
#include "orbital_core/perturbations/AsteroidPerturbations.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace orbital_core {
	namespace {

		int numberedBody(std::string_view text, bool packed)
		{
			const auto fold = [](char c) { return c >= 'A' && c <= 'Z' ? c + ('a' - 'A') : c; };
			const auto matches = [&](std::string_view name) {
				return text.size() == name.size() && std::equal(text.begin(), text.end(), name.begin(),
				                                                [&](char a, char b) { return fold(a) == fold(b); });
			};
			if (matches("Moon") || (packed && text == "E001S"))
				return 301;
			if (matches("Pluto"))
				return 999;
			for (const auto& body : forces::detail::kSb441Asteroids) {
				if (matches(body.name))
					return body.naifId;
			}
			int number = 0;
			if (packed && text.size() == 5 &&
			    ((text.front() >= 'A' && text.front() <= 'Z') || (text.front() >= 'a' && text.front() <= 'z'))) {
				const int prefix = text.front() <= 'Z' ? text.front() - 'A' + 10 : text.front() - 'a' + 36;
				const auto digits = text.substr(1);
				const auto [end, error] = std::from_chars(digits.data(), digits.data() + digits.size(), number);
				if (error != std::errc{} || end != digits.data() + digits.size())
					return 0;
				number += prefix * 10000;
			} else {
				const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), number);
				if (error != std::errc{} || end != text.data() + text.size())
					return 0;
			}
			if (number == 134340)
				return 999;
			return number > 0 && number <= 999999 ? 2000000 + number : 0;
		}

		std::string readFile(const std::filesystem::path& path)
		{
			FileReader reader;
			reader.read(path);
			if (reader.error() != FileReader::Error::NoError)
				throw std::runtime_error(reader.errorString() + ": " + path.string());
			return reader.data();
		}

		std::vector<ObsNode> makeObservations(std::span<const obs::ParsedObservation> parsed,
		                                      const obs::ObservatoryCatalog& catalog)
		{
			std::vector<ObsNode> result;
			result.reserve(parsed.size());
			for (const auto& observation : parsed) {
				const auto station = catalog.find(observation.obsCode);
				if (!observation.radar && station == catalog.end())
					throw std::runtime_error("Unknown observatory code: " + observation.obsCode);
				try {
					result.emplace_back(
					    observation, observation.radar ? obs::radarObservatory(observation.obsCode) : station->second);
				} catch (const std::string& message) {
					throw std::runtime_error(message);
				}
			}
			return result;
		}

		void selectIodPool(std::vector<obs::ParsedObservation>& parsed)
		{
			std::erase_if(parsed, [](const auto& observation) { return observation.radar.has_value(); });
			if (parsed.size() < 3)
				throw std::runtime_error("Insufficient data");
			const auto end =
			    std::unique(parsed.begin(), parsed.end(), [](const auto& a, const auto& b) { return a == b; });
			parsed.erase(end, parsed.end());
			const auto satellites = std::ranges::count(parsed, true, &obs::ParsedObservation::isSatellite);
			LOG_INFO << "findBestObsNodes: IOD pool ready — " << parsed.size() << " astrometric obs over "
			         << (parsed.back().jd - parsed.front().jd) << "d total span";
			if (satellites > 0) {
				LOG_INFO << "findBestObsNodes: retained " << satellites << " spacecraft obs in IOD pool";
			}
			if (parsed.size() < 3)
				throw std::runtime_error("Insufficient data");
		}

		void restoreSubmittedCoordinates(ObservationSet& observations)
		{
			for (auto* pool : {&observations.all, &observations.iod}) {
				for (auto& o : *pool) {
					if (o.radar || !o.submittedCoordinates)
						continue;
					o.RA = o.submittedRaRad();
					o.DEC = o.submittedDecRad();
					o.vector = o.submittedCoordinates->directionEcl;
					o.submittedCoordinates.reset();
					o.biasRaArcsec = 0.0;
					o.biasDecArcsec = 0.0;
				}
			}
			observations.preparation.reset();
		}

	}  // namespace

	ObservationSet loadObservations(const std::filesystem::path& path)
	{
		if (path.empty())
			throw std::invalid_argument("Observation input path is empty");
		const auto payload = readFile(path);
		const auto* reader = obs::selectObservationReader(payload);
		if (!reader)
			throw std::runtime_error("unrecognised observation format: " + path.string());

		std::string error;
		auto parsed = reader->read(payload, error);
		if (!error.empty())
			throw std::runtime_error(error);
		for (auto& observation : parsed) {
			observation.id = newObservationId();
			if (observation.radar) {
				observation.obsCode = obs::radarObservatory(observation.obsCode).code;
				observation.radar->transmitter = obs::radarObservatory(observation.radar->transmitter).code;
			}
		}
		const auto withSigma = std::ranges::count(parsed, true, &obs::ParsedObservation::hasObserverSigma);
		LOG_INFO << "Observation input: " << reader->formatName() << ", " << parsed.size() << " records, " << withSigma
		         << " carrying observer sigmas";

		ObservationSet result;
		result.source = path;
		result.format = reader->formatName();
		bool hasGaia = false;
		for (const auto& observation : parsed) {
			if (observation.astrometricConvention == obs::AstrometricConvention::GaiaFpr) {
				if (hasGaia && observation.gaiaPhotocenterRadiusKm != result.gaiaPhotocenterRadiusKm)
					throw std::runtime_error("Inconsistent Gaia photocenter radius in observations");
				hasGaia = true;
				result.gaiaPhotocenterRadiusKm = observation.gaiaPhotocenterRadiusKm;
			}
			const bool mpc = result.format == "MPC 80-column";
			const auto& designation = observation.permanentDesignation.empty() && mpc
			                              ? observation.objectName
			                              : observation.permanentDesignation;
			const int id = numberedBody(designation, mpc && observation.permanentDesignation.empty());
			if (id && result.targetNaifId && result.targetNaifId != id)
				throw std::runtime_error("Observation input contains conflicting numbered bodies");
			if (id)
				result.targetNaifId = id;
		}
		const auto named =
		    std::ranges::find_if(parsed, [](const auto& observation) { return !observation.objectName.empty(); });
		if (named != parsed.end())
			result.objectName = named->objectName;
		const auto designated = std::ranges::find_if(
		    parsed, [](const auto& observation) { return !observation.classificationDesignation.empty(); });
		if (designated != parsed.end())
			result.classificationDesignation = designated->classificationDesignation;
		std::sort(parsed.begin(), parsed.end(), [](const auto& a, const auto& b) { return a.jd < b.jd; });

		const auto& catalog = obs::observatoryCatalog();
		result.all = makeObservations(parsed, catalog);
		result.coincidentObservations = findCoincidentObservations(result.all);
		result.visibilityDiagnostics = diagnoseObservationVisibility(result.all);
		selectIodPool(parsed);
		result.iod = makeObservations(parsed, catalog);
		return result;
	}

	ObservationSet withGaiaFprObservations(const ObservationSet& source, std::span<const GaiaFprObservation> gaia,
	                                       std::optional<double> photocenterRadiusKm)
	{
		if (gaia.empty())
			throw std::invalid_argument("Gaia FPR returned no observations");
		if (photocenterRadiusKm && (!std::isfinite(*photocenterRadiusKm) || *photocenterRadiusKm <= 0.0))
			throw std::invalid_argument("Gaia photocenter radius must be positive and finite");
		std::vector<obs::ParsedObservation> parsed;
		parsed.reserve(gaia.size());
		for (const auto& row : gaia) {
			if (!std::isfinite(row.epochMjdUtc) || !std::isfinite(row.raDeg) || !std::isfinite(row.decDeg) ||
			    row.raDeg < 0.0 || row.raDeg >= 360.0 || std::abs(row.decDeg) > 90.0 ||
			    !std::isfinite(row.sigmaRaArcsec) || !std::isfinite(row.sigmaDecArcsec) || !(row.sigmaRaArcsec > 0.0) ||
			    !(row.sigmaDecArcsec > 0.0) || !std::isfinite(row.correlation) || std::abs(row.correlation) >= 1.0 ||
			    !std::isfinite(row.observerOffsetEquatorialAu.x) || !std::isfinite(row.observerOffsetEquatorialAu.y) ||
			    !std::isfinite(row.observerOffsetEquatorialAu.z))
				throw std::invalid_argument("Invalid Gaia FPR observation");
			obs::ParsedObservation o;
			o.id = newObservationId();
			o.objectName = source.objectName;
			o.obsCode = "258";
			o.jd = row.epochMjdUtc + 2400000.5;
			o.utcMjdDay = static_cast<int>(std::floor(row.epochMjdUtc));
			o.utcSecondsOfDay = (row.epochMjdUtc - o.utcMjdDay) * 86400.0;
			o.hasSplitUtc = true;
			o.raRad = row.raDeg * Constants::PI / 180.0;
			o.decRad = row.decDeg * Constants::PI / 180.0;
			o.technique = obs::Technique::Ccd;
			o.catalogCode = 'W';
			o.astCat = "Gaia3";
			o.astrometricConvention = obs::AstrometricConvention::GaiaFpr;
			o.hasObserverSigma = true;
			o.rmsRaArcsec = row.sigmaRaArcsec;
			o.rmsDecArcsec = row.sigmaDecArcsec;
			o.rmsCorr = row.correlation;
			o.isSatellite = true;
			o.offsetUnits = obs::OffsetUnits::Au;
			o.offsetX = row.observerOffsetEquatorialAu.x;
			o.offsetY = row.observerOffsetEquatorialAu.y;
			o.offsetZ = row.observerOffsetEquatorialAu.z;
			parsed.push_back(std::move(o));
		}
		const auto [first, last] = std::ranges::minmax_element(gaia, {}, &GaiaFprObservation::epochMjdUtc);
		for (const auto& existing : source.all) {
			if (existing.astrometricConvention == obs::AstrometricConvention::GaiaFpr)
				throw std::invalid_argument("Input already contains Gaia FPR observations");
			if (existing.obsCode != "258")
				continue;
			// A transit may span multiple CCD measurements and UTC dates. Without
			// transit identifiers in MPC/ADES, reject potentially overlapping data.
			if (existing.timeScale != time::Scale::UTC)
				throw std::invalid_argument("Cannot check existing Gaia observation time against FPR");
			const double epoch = existing.hasSplitUtc ? existing.utcMjdDay + existing.utcSecondsOfDay / 86400.0
			                                          : existing.JD - 2400000.5;
			const double margin = existing.hasTimeAccuracy ? std::max(1.0, existing.timeAccuracyDays) : 1.0;
			if (!std::isfinite(epoch) || !std::isfinite(margin) ||
			    (epoch >= first->epochMjdUtc - margin && epoch <= last->epochMjdUtc + margin))
				throw std::invalid_argument(
				    "Existing Gaia observations may overlap FPR; resolve duplicates before adding it");
		}
		auto result = source;
		restoreSubmittedCoordinates(result);
		auto added = makeObservations(parsed, obs::observatoryCatalog());
		result.all.insert(result.all.end(), added.begin(), added.end());
		result.iod.insert(result.iod.end(), added.begin(), added.end());
		const auto byTime = [](const ObsNode& a, const ObsNode& b) { return a.JD < b.JD; };
		std::stable_sort(result.all.begin(), result.all.end(), byTime);
		std::stable_sort(result.iod.begin(), result.iod.end(), byTime);
		result.gaiaPhotocenterRadiusKm = photocenterRadiusKm;
		result.format += " + Gaia FPR";
		LOG_INFO << "Gaia FPR: added " << gaia.size() << " CCD observations; transit systematics "
		         << "use a per-observation approximation";
		result.coincidentObservations = findCoincidentObservations(result.all);
		result.visibilityDiagnostics = diagnoseObservationVisibility(result.all);
		return result;
	}

	void prepareObservations(ObservationSet& observations, const ObservationPreparationOptions& options)
	{
		if (observations.preparation)
			throw std::logic_error("Observations are already prepared; reload before changing preparation");
		ObservationPreparation result;
		result.debiasRequested = options.debias;
		result.weightProfile = options.weightProfile;
		if (options.debias) {
			const auto& directory = options.runtimeDirectory.empty() ? getExecutableDir() : options.runtimeDirectory;
			result.biasTable = Calculations::astrometryBiasPathForProfile(directory, options.weightProfile);
			result.correctedAll = Calculations::debiasObservations(observations.all, result.biasTable);
			result.correctedIod = Calculations::debiasObservations(observations.iod, result.biasTable);
		}
		observations.preparation = std::move(result);
	}

	ObservationSet prepareObservationCopy(const ObservationSet& source, const ObservationPreparationOptions& options)
	{
		auto result = source;
		restoreSubmittedCoordinates(result);
		prepareObservations(result, options);
		return result;
	}

}  // namespace orbital_core
