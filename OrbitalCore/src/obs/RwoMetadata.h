#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/Constants.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>

namespace orbital_core::obs {

	inline time::Scale rwoTimeScale(int value)
	{
		if (value < static_cast<int>(time::Scale::UTC) || value > static_cast<int>(time::Scale::UT1))
			throw std::invalid_argument("Invalid RWO time scale metadata");
		return static_cast<time::Scale>(value);
	}

	// Plain RWO columns follow the legacy UT1/UTC convention. Metadata retains
	// the submitted epoch and scale when ADES or API input uses a different one.
	inline double rwoEpochJd(double jd, time::Scale sourceScale)
	{
		const auto target = time::legacyObservationScale(jd);
		if (sourceScale == target)
			return jd;
		const double tt = time::toMjdTt(time::mjdFromJd(jd), sourceScale);
		const double converted = time::jdFromMjd(time::fromMjdTt(tt, target));
		if (!std::isfinite(converted))
			throw std::invalid_argument("RWO epoch time conversion unavailable");
		return converted;
	}

	// Standard fixed columns remain authoritative for other readers. These versioned
	// comments retain input precision and ADES metadata that RWO cannot represent.
	inline nlohmann::json rwoInputMetadata(const ObsNode& obs, std::optional<double> gaiaPhotocenterRadiusKm = {})
	{
		nlohmann::json result{
		    {"v", 1},
		    {"station", obs.obsCode},
		    {"jd", obs.JD},
		    {"splitUtc", {obs.hasSplitUtc, obs.utcMjdDay, obs.utcSecondsOfDay}},
		    {"ra", obs.submittedRaRad()},
		    {"dec", obs.submittedDecRad()},
		    {"astCat", obs.astCat},
		    {"program", obs.adesProgramCode},
		    {"notes", obs.notes},
		    {"note2", std::string(1, obs.note2)},
		    {"magnitude", obs.magString},
		    {"rmsMag", obs.rmsMag},
		    {"precision", {obs.hasTimeAccuracy, obs.timeAccuracyDays, obs.raAccuracyArcsec, obs.decAccuracyArcsec}},
		    {"sigma", {obs.hasObserverSigma, obs.rmsRaArcsec, obs.rmsDecArcsec, obs.rmsCorr}},
		    {"rmsTime", {obs.hasRmsTime, obs.rmsTimeSec}},
		    {"uncTime", {obs.hasUncTime, obs.uncTimeSec}}};
		if (obs.timeScale != time::legacyObservationScale(obs.JD))
			result["timeScale"] = static_cast<int>(obs.timeScale);
		if (obs.astrometricConvention == AstrometricConvention::Absolute)
			result["astrometricConvention"] = "absolute";
		if (obs.astrometricConvention == AstrometricConvention::GaiaFpr) {
			result["astrometricConvention"] = "gaia_fpr";
			if (gaiaPhotocenterRadiusKm) {
				if (!std::isfinite(*gaiaPhotocenterRadiusKm) || *gaiaPhotocenterRadiusKm <= 0)
					throw std::invalid_argument("Invalid Gaia photocenter radius");
				result["gaiaPhotocenterRadiusKm"] = *gaiaPhotocenterRadiusKm;
			}
		}
		if (obs.nuclearMagnitude)
			result["nuclearMagnitude"] = *obs.nuclearMagnitude;
		if (obs.rovingObserver)
			result["roving"] = {obs.rovingObserver->longitudeDeg, obs.rovingObserver->latitudeDeg,
			                    obs.rovingObserver->heightMeters};
		return result;
	}

	inline double finiteRwoNumber(const nlohmann::json& value)
	{
		const double result = value.get<double>();
		if (!std::isfinite(result))
			throw std::invalid_argument("Non-finite RWO metadata");
		return result;
	}

	inline void restoreRwoInputMetadata(const nlohmann::json& data, ParsedObservation& obs)
	{
		if (data.at("v").get<int>() != 1 || data.at("station").get<std::string>() != obs.obsCode)
			throw std::invalid_argument("RWO input metadata version or station mismatch");
		const double jd = finiteRwoNumber(data.at("jd"));
		const auto scale = rwoTimeScale(data.value("timeScale", static_cast<int>(obs.timeScale)));
		const double ra = finiteRwoNumber(data.at("ra")), dec = finiteRwoNumber(data.at("dec"));
		if (std::abs(rwoEpochJd(jd, scale) - obs.jd) > 0.5 / 86400.0 || ra < 0 || ra >= 2 * Constants::PI ||
		    std::abs(dec) > Constants::PI / 2 || std::abs(std::remainder(ra - obs.raRad, 2 * Constants::PI)) > 4e-8 ||
		    std::abs(dec - obs.decRad) > 3e-8)
			throw std::invalid_argument("RWO input metadata disagrees with coordinates or epoch");
		obs.jd = jd;
		obs.timeScale = scale;
		obs.raRad = ra;
		obs.decRad = dec;
		const auto& utc = data.at("splitUtc");
		obs.hasSplitUtc = utc.at(0).get<bool>();
		obs.utcMjdDay = utc.at(1).get<int>();
		obs.utcSecondsOfDay = finiteRwoNumber(utc.at(2));
		if (obs.hasSplitUtc && (obs.utcSecondsOfDay < 0 || obs.utcSecondsOfDay >= 86400 ||
		                        std::abs((obs.utcMjdDay + 2400000.5 + obs.utcSecondsOfDay / 86400.0) - jd) > 1e-9))
			throw std::invalid_argument("Invalid split UTC in RWO metadata");
		obs.astCat = data.at("astCat").get<std::string>();
		if (data.contains("astrometricConvention")) {
			const auto convention = data.at("astrometricConvention").get<std::string>();
			if (convention == "absolute")
				obs.astrometricConvention = AstrometricConvention::Absolute;
			else if (convention == "gaia_fpr")
				obs.astrometricConvention = AstrometricConvention::GaiaFpr;
			else
				throw std::invalid_argument("Invalid RWO astrometric convention");
		}
		if (data.contains("gaiaPhotocenterRadiusKm")) {
			const double radius = finiteRwoNumber(data.at("gaiaPhotocenterRadiusKm"));
			if (obs.astrometricConvention != AstrometricConvention::GaiaFpr || radius <= 0)
				throw std::invalid_argument("Invalid RWO Gaia photocenter radius");
			obs.gaiaPhotocenterRadiusKm = radius;
		}
		obs.adesProgramCode = data.at("program").get<std::string>();
		obs.notes = data.at("notes").get<std::string>();
		if (data.contains("nuclearMagnitude")) {
			if (!data.at("nuclearMagnitude").is_boolean())
				throw std::invalid_argument("Invalid RWO comet magnitude type metadata");
			const bool nuclear = data.at("nuclearMagnitude").get<bool>();
			if (obs.nuclearMagnitude && *obs.nuclearMagnitude != nuclear)
				throw std::invalid_argument("RWO comet magnitude type disagrees with metadata");
			obs.nuclearMagnitude = nuclear;
		}
		if (data.contains("note2")) {
			const auto note = data.at("note2").get<std::string>();
			if (note.size() != 1)
				throw std::invalid_argument("Invalid RWO technique metadata");
			bool found = false;
			for (int i = 0; i <= static_cast<int>(Technique::ReplacedPhotographic); ++i) {
				const auto technique = static_cast<Technique>(i);
				if (mpcNote2(technique) == note.front()) {
					obs.technique = technique;
					found = true;
					break;
				}
			}
			if (!found &&
			    ((note == "S" && obs.isSatellite) || (note == "V" && (obs.rovingObserver || obs.isSatellite)))) {
				obs.technique = Technique::Ccd;
				found = true;
			}
			if (!found)
				throw std::invalid_argument("Unsupported RWO technique metadata");
		}
		if (data.contains("roving")) {
			const auto& position = data.at("roving");
			const RovingObserver rover{finiteRwoNumber(position.at(0)), finiteRwoNumber(position.at(1)),
			                           finiteRwoNumber(position.at(2))};
			if (!obs.rovingObserver || !rover.valid() ||
			    std::abs(std::remainder(rover.longitudeDeg - obs.rovingObserver->longitudeDeg, 360.0)) > 0.00000051 ||
			    std::abs(rover.latitudeDeg - obs.rovingObserver->latitudeDeg) > 0.00000051 ||
			    std::abs(rover.heightMeters - obs.rovingObserver->heightMeters) > 0.051)
				throw std::invalid_argument("Roving metadata disagrees with RWO position");
			obs.rovingObserver = rover;
		}
		obs.magnitudeText = data.at("magnitude").get<std::string>();
		obs.rmsMag = finiteRwoNumber(data.at("rmsMag"));
		const auto& precision = data.at("precision");
		obs.hasTimeAccuracy = precision.at(0).get<bool>();
		obs.timeAccuracyDays = finiteRwoNumber(precision.at(1));
		obs.raAccuracyArcsec = finiteRwoNumber(precision.at(2));
		obs.decAccuracyArcsec = finiteRwoNumber(precision.at(3));
		const auto& sigma = data.at("sigma");
		obs.hasObserverSigma = sigma.at(0).get<bool>();
		obs.rmsRaArcsec = finiteRwoNumber(sigma.at(1));
		obs.rmsDecArcsec = finiteRwoNumber(sigma.at(2));
		obs.rmsCorr = finiteRwoNumber(sigma.at(3));
		obs.hasRmsTime = data.at("rmsTime").at(0).get<bool>();
		obs.rmsTimeSec = finiteRwoNumber(data.at("rmsTime").at(1));
		obs.hasUncTime = data.at("uncTime").at(0).get<bool>();
		obs.uncTimeSec = finiteRwoNumber(data.at("uncTime").at(1));
		if ((obs.hasObserverSigma && (obs.rmsRaArcsec <= 0 || obs.rmsDecArcsec <= 0 || std::abs(obs.rmsCorr) >= 1)) ||
		    obs.rmsMag < 0 || obs.rmsTimeSec < 0 || obs.uncTimeSec < 0 ||
		    (obs.hasTimeAccuracy && obs.timeAccuracyDays <= 0))
			throw std::invalid_argument("Invalid RWO uncertainty metadata");
	}
}  // namespace orbital_core::obs
