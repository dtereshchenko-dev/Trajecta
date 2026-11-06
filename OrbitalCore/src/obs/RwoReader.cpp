#include "obs/RwoReader.h"
#include "obs/RwoMetadata.h"
#include "obs/ObservatoryCatalog.h"
#include "orbital_core/Types.h"
#include <chrono>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace orbital_core::obs {
	namespace {
		std::string trimmed(std::string value)
		{
			const auto first = value.find_first_not_of(" \t\r\n");
			return first == std::string::npos ? std::string{}
			                                  : value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
		}

		std::string field(const std::string& line, std::size_t offset, std::size_t width)
		{
			if (line.size() < offset + width)
				throw std::invalid_argument("Truncated RWO record");
			return trimmed(line.substr(offset, width));
		}

		double numeric(std::string token)
		{
			std::replace(token.begin(), token.end(), 'D', 'E');
			std::replace(token.begin(), token.end(), 'd', 'e');
			std::istringstream stream(token);
			stream.imbue(std::locale::classic());
			double result;
			if (!(stream >> result) || !stream.eof() || !std::isfinite(result))
				throw std::invalid_argument("Invalid RWO number: " + token);
			return result;
		}

		struct Record {
			const std::string& line;
			std::size_t shift;

			std::string field(std::size_t offset, std::size_t width) const
			{
				return obs::field(line, offset + shift, width);
			}

			double at(std::size_t offset, std::size_t width) const { return numeric(field(offset, width)); }

			char character(std::size_t offset) const { return line.at(offset + shift); }
		};

		bool selected(const std::string& value)
		{
			if (value != "0" && value != "1")
				throw std::invalid_argument("Invalid RWO selection flag");
			return value == "1";
		}

		void date(ParsedObservation& obs, const Record& record, bool radar)
		{
			using namespace std::chrono;
			const int y = static_cast<int>(record.at(17, 4)), m = static_cast<int>(record.at(22, 2));
			const double d = record.at(25, radar ? 2 : 13);
			const year_month_day calendar{year{y}, month{static_cast<unsigned>(m)}, day{static_cast<unsigned>(d)}};
			if (!calendar.ok())
				throw std::invalid_argument("Invalid RWO date");
			obs.utcMjdDay = static_cast<int>(sys_days(calendar).time_since_epoch().count()) + 40587;
			obs.utcSecondsOfDay = radar ? record.at(28, 2) * 3600 + record.at(31, 2) * 60 + record.at(34, 2)
			                            : (d - std::floor(d)) * 86400;
			if (!(obs.utcSecondsOfDay >= 0 && obs.utcSecondsOfDay < 86400))
				throw std::invalid_argument("Invalid RWO time of day");
			obs.hasSplitUtc = true;
			obs.jd = 2400000.5 + (obs.utcMjdDay + obs.utcSecondsOfDay / 86400.0);
			obs.timeScale = time::legacyObservationScale(obs.jd);
			obs.hasTimeAccuracy = true;
			obs.timeAccuracyDays = radar ? 1.0 / 86400.0 : record.at(39, 10);
		}

		Technique technique(char code)
		{
			for (int i = 0; i <= static_cast<int>(Technique::ReplacedPhotographic); ++i)
				if (mpcNote2(static_cast<Technique>(i)) == code)
					return static_cast<Technique>(i);
			if (code == 'S' || code == 'V')
				return Technique::Ccd;
			throw std::invalid_argument("Unsupported RWO technique");
		}

		struct RadarMetadata {
			double frequency, radius, observed, sigma;
			bool reference;
			std::optional<double> jd;
			std::optional<time::Scale> timeScale;
			bool splitUtc{};
			int utcMjdDay{};
			double utcSeconds{};
		};
	}  // namespace

	bool RwoReader::canRead(std::string_view text) const noexcept
	{
		return text.substr(0, 4096).find("version") != std::string_view::npos &&
		       text.substr(0, 4096).find("END_OF_HEADER") != std::string_view::npos;
	}

	std::vector<ParsedObservation> RwoReader::read(const std::string& text, std::string& error) const
	{
		error.clear();
		std::vector<ParsedObservation> observations;
		std::istringstream stream(text);
		std::size_t lineNumber = 0;
		try {
			bool header = false;
			int version = 0;
			std::optional<std::string> objectName;
			std::optional<RadarMetadata> metadata;
			std::optional<nlohmann::json> inputMetadata;
			for (std::string line; std::getline(stream, line);) {
				++lineNumber;
				if (!line.empty() && line.back() == '\r')
					line.pop_back();
				if (!header) {
					if (line.starts_with("! OC_OBJECT ")) {
						if (objectName)
							throw std::invalid_argument("Duplicate RWO object metadata");
						objectName = nlohmann::json::parse(line.substr(12)).get<std::string>();
					}
					if (line.starts_with("version")) {
						const auto eq = line.find('=');
						if (version || eq == std::string::npos)
							throw std::invalid_argument("Invalid RWO version header");
						const double value = numeric(trimmed(line.substr(eq + 1)));
						if (value != 2 && value != 3)
							throw std::invalid_argument("Only RWO versions 2 and 3 are supported");
						version = static_cast<int>(value);
					}
					if (trimmed(line) == "END_OF_HEADER")
						header = true;
					continue;
				}
				if (trimmed(line).empty())
					continue;
				if (line.starts_with("! OC_INPUT ")) {
					if (inputMetadata)
						throw std::invalid_argument("Duplicate RWO input metadata");
					inputMetadata = nlohmann::json::parse(line.substr(11));
					continue;
				}
				if (line.starts_with("! OC_RADAR ")) {
					std::istringstream info(line.substr(11));
					info.imbue(std::locale::classic());
					RadarMetadata data{};
					if (metadata ||
					    !(info >> data.frequency >> data.reference >> data.radius >> data.observed >> data.sigma))
						throw std::invalid_argument("Malformed radar metadata");
					info >> std::ws;
					if (!info.eof()) {
						double jd;
						if (!(info >> jd >> data.splitUtc >> data.utcMjdDay >> data.utcSeconds) || !std::isfinite(jd) ||
						    !std::isfinite(data.utcSeconds))
							throw std::invalid_argument("Malformed radar epoch metadata");
						data.jd = jd;
						info >> std::ws;
						if (!info.eof()) {
							int scale;
							if (!(info >> scale))
								throw std::invalid_argument("Malformed radar time scale");
							data.timeScale = rwoTimeScale(scale);
							info >> std::ws;
						}
						if (!info.eof())
							throw std::invalid_argument("Trailing radar metadata fields");
					}
					if (!std::isfinite(data.radius) || data.radius < 0 || !std::isfinite(data.frequency) ||
					    data.frequency <= 0)
						throw std::invalid_argument("Invalid radar frequency or radius");
					metadata = data;
					continue;
				}
				if (line.front() == '!')
					continue;
				if (!version)
					throw std::invalid_argument("Missing RWO version");
				const Record record{line, version == 3 ? 3u : 0u};
				ParsedObservation obs;
				obs.objectName = field(line, 1, 9 + record.shift);
				if (objectName)
					obs.objectName = *objectName;
				obs.sourceFormat = formatName();
				const char type = record.character(11), tech = record.character(13);
				obs.programCode = record.character(15);
				obs.importedFit.emplace();
				auto& archived = *obs.importedFit;
				if (type == 'R' || type == 'V') {
					if (inputMetadata)
						throw std::invalid_argument("Optical metadata before radar record");
					date(obs, record, true);
					RadarObservation radar;
					std::istringstream tail(line.substr(99 + record.shift));
					std::string transmitter, receiver;
					if (!(tail >> transmitter >> receiver))
						throw std::invalid_argument("Missing RWO radar stations");
					radar.transmitter = radarObservatory(transmitter).code;
					obs.obsCode = radarObservatory(receiver).code;
					std::vector<std::string> diagnostics;
					for (std::string token; tail >> token;)
						diagnostics.push_back(std::move(token));
					if (diagnostics.empty() || diagnostics.size() > 2)
						throw std::invalid_argument("Invalid RWO radar diagnostics");
					archived.selected = selected(diagnostics.back());
					if (diagnostics.size() == 2 && numeric(diagnostics.front()) < 0)
						throw std::invalid_argument("Invalid RWO radar chi");
					if (tech != 'c' && tech != 's')
						throw std::invalid_argument("Invalid radar return type");
					radar.returnType = tech == 's' ? RadarReturn::Surface : RadarReturn::CenterOfMass;
					radar.frequencyHz = metadata ? metadata->frequency : 1.0;
					radar.frequencyIsReference = metadata ? metadata->reference : true;
					if (metadata && metadata->radius > 0)
						radar.surfaceRadiusKm = metadata->radius;
					constexpr double cKmSec = Constants::AU_PER_DAY * Constants::AU_IN_METERS / 1000.0 / 86400.0;
					const double factor = type == 'R' ? cKmSec / 2.0 : -cKmSec * 86400 / (2 * radar.frequencyHz);
					RadarMeasurement channel{record.at(36, 16) / factor, record.at(63, 9) / std::abs(factor)};
					if (metadata) {
						if (std::abs((metadata->observed - channel.value) * factor) >
						        std::max(5.1e-6, std::abs(channel.value * factor) * 5e-6) ||
						    std::abs((metadata->sigma - channel.sigma) * factor) >
						        std::max(5.1e-6, std::abs(channel.sigma * factor) * 0.006))
							throw std::invalid_argument("Radar metadata disagrees with RWO measurement");
						channel = {metadata->observed, metadata->sigma};
						if (metadata->jd) {
							const auto scale = metadata->timeScale.value_or(obs.timeScale);
							if (std::abs(rwoEpochJd(*metadata->jd, scale) - obs.jd) > 0.51 / 86400.0 ||
							    (metadata->splitUtc &&
							     (metadata->utcSeconds < 0 || metadata->utcSeconds >= 86400 ||
							      std::abs(metadata->utcMjdDay + 2400000.5 + metadata->utcSeconds / 86400.0 -
							               *metadata->jd) > 1e-9)))
								throw std::invalid_argument("Radar metadata disagrees with RWO epoch");
							obs.jd = *metadata->jd;
							obs.timeScale = scale;
							obs.hasSplitUtc = metadata->splitUtc;
							obs.utcMjdDay = metadata->utcMjdDay;
							obs.utcSecondsOfDay = metadata->utcSeconds;
						}
					}
					if (type == 'R')
						radar.delaySeconds = channel;
					else
						radar.dopplerHz = channel;
					if (!radar.valid())
						throw std::invalid_argument("Invalid radar measurement");
					archived.forceWeight1 = record.character(73) == 'T';
					archived.bias1 = record.at(75, 11) / factor;
					if (!record.field(87, 11).empty())
						archived.residual1 = record.at(87, 11) / factor;
					obs.radar = radar;
					metadata.reset();
				} else if (type == 'O' || type == 'S') {
					if (metadata)
						throw std::invalid_argument("Radar metadata before optical record");
					date(obs, record, false);
					obs.technique = technique(tech);
					obs.raRad = Convertors::raToRadians(static_cast<int>(record.at(50, 2)),
					                                    static_cast<int>(record.at(53, 2)), record.at(56, 6));
					obs.decRad = Convertors::decToRadians(static_cast<int>(record.at(104, 2)),
					                                      static_cast<int>(record.at(107, 2)), record.at(110, 5));
					if (record.character(103) == '-')
						obs.decRad = -obs.decRad;
					obs.raAccuracyArcsec = record.at(63, 10);
					obs.decAccuracyArcsec = record.at(116, 10);
					obs.rmsRaArcsec = record.at(74, 8);
					obs.rmsDecArcsec = record.at(127, 8);
					obs.hasObserverSigma = obs.rmsRaArcsec > 0 && obs.rmsDecArcsec > 0;
					obs.magnitudeText = record.field(156, 5);
					obs.magBand = record.character(161);
					if (obs.magBand == 'N' || obs.magBand == 'T')
						obs.nuclearMagnitude = obs.magBand == 'N';
					if (!record.field(163, 5).empty())
						obs.rmsMag = record.at(163, 5);
					const auto catalog = record.field(177, 2);
					obs.catalogCode = catalog.empty() ? ' ' : catalog.front();
					obs.obsCode = record.field(180, 3);
					archived.selected = selected(record.field(194, 1));
					archived.bias1 = record.at(85, 8);
					archived.bias2 = record.at(138, 8);
					archived.forceWeight1 = record.character(83) == 'T';
					archived.forceWeight2 = record.character(136) == 'T';
					if (!record.field(93, 9).empty())
						archived.residual1 = record.at(93, 9);
					if (!record.field(146, 9).empty())
						archived.residual2 = record.at(146, 9);
					if (type == 'S') {
						std::string position;
						if (!std::getline(stream, position))
							throw std::invalid_argument("Missing satellite position");
						++lineNumber;
						const Record observer{position, record.shift};
						if (position.size() < 35 + record.shift || observer.character(13) != 's' ||
						    position.substr(0, 13 + record.shift) != line.substr(0, 13 + record.shift) ||
						    position.substr(14 + record.shift, 19) != line.substr(14 + record.shift, 19))
							throw std::invalid_argument("Satellite continuation does not match observation");
						obs.isSatellite = true;
						obs.offsetUnits = observer.character(34) == '1'   ? OffsetUnits::Km
						                  : observer.character(34) == '2' ? OffsetUnits::Au
						                                                  : OffsetUnits::None;
						if (obs.offsetUnits == OffsetUnits::None)
							throw std::invalid_argument("Invalid observer coordinate units");
						std::istringstream fields(position.substr(35 + record.shift));
						fields.imbue(std::locale::classic());
						std::string station;
						if (!(fields >> obs.offsetX >> obs.offsetY >> obs.offsetZ >> station) ||
						    station != obs.obsCode || !std::isfinite(obs.offsetX) || !std::isfinite(obs.offsetY) ||
						    !std::isfinite(obs.offsetZ))
							throw std::invalid_argument("Invalid satellite position");
					} else if (tech == 'V') {
						std::string position;
						if (!std::getline(stream, position))
							throw std::invalid_argument("Missing roving position");
						++lineNumber;
						const Record observer{position, record.shift};
						if (position.size() < 33 + record.shift || observer.character(13) != 'v' ||
						    position.substr(0, 13 + record.shift) != line.substr(0, 13 + record.shift) ||
						    position.substr(14 + record.shift, 19) != line.substr(14 + record.shift, 19))
							throw std::invalid_argument("Roving continuation does not match observation");
						RovingObserver rover;
						std::istringstream fields(position.substr(33 + record.shift));
						fields.imbue(std::locale::classic());
						std::string station, trailing;
						if (!(fields >> rover.longitudeDeg >> rover.latitudeDeg >> rover.heightMeters >> station) ||
						    station != obs.obsCode || (fields >> trailing) || !rover.valid())
							throw std::invalid_argument("Invalid roving position");
						obs.rovingObserver = rover;
					}
				} else {
					throw std::invalid_argument("Unsupported RWO record type");
				}
				if (inputMetadata) {
					restoreRwoInputMetadata(*inputMetadata, obs);
					inputMetadata.reset();
				}
				observations.push_back(std::move(obs));
			}
			if (!version || !header || metadata || inputMetadata)
				throw std::invalid_argument("Incomplete RWO file");
			if (observations.empty())
				throw std::invalid_argument("No observations in RWO file");
		} catch (const std::exception& e) {
			error = "RWO line " + std::to_string(lineNumber) + ": " + e.what();
			return {};
		}
		return observations;
	}
}  // namespace orbital_core::obs
