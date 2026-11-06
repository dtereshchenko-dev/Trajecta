#include "AdesParser.h"
#include <string>
#include <utility>
#include <vector>

#include "orbital_core/Constants.h"
#include "orbital_core/obs/ParsedObservation.h"
#include "orbital_core/obs/StarCatalog.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <string_view>
#include <unordered_map>

namespace orbital_core::ades {
	namespace {

		[[nodiscard]] obs::Technique techniqueFromMode(std::string_view mode) noexcept
		{
			if (mode == "CCD" || mode == "TDI")
				return obs::Technique::Ccd;
			if (mode == "CMO")
				return obs::Technique::Cmos;
			if (mode == "VID")
				return obs::Technique::Video;
			if (mode == "PHO")
				return obs::Technique::Photographic;
			if (mode == "ENC")
				return obs::Technique::Encoder;
			if (mode == "MIC")
				return obs::Technique::Micrometer;
			if (mode == "MER")
				return obs::Technique::Meridian;
			if (mode == "OCC")
				return obs::Technique::Occultation;
			if (mode == "PMT")
				return obs::Technique::Photometer;
			return obs::Technique::Unknown;
		}

		[[nodiscard]] char photometricBand(std::string_view band) noexcept
		{
			if (band.empty())
				return ' ';
			if (band.size() > 1 &&
			    (((band[0] == 'P' || band[0] == 'S') && std::string_view("grizwy").find(band[1]) != band.npos) ||
			     (band[0] == 'A' && (band[1] == 'o' || band[1] == 'c')) ||
			     (band[0] == 'G' && (band[1] == 'b' || band[1] == 'r'))))
				return band[1];
			return band.front();
		}

		[[nodiscard]] std::string trim(std::string_view text)
		{
			std::size_t begin = 0;
			std::size_t end = text.size();
			while (begin < end && std::isspace(static_cast<unsigned char>(text[begin])))
				++begin;
			while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1])))
				--end;
			return std::string(text.substr(begin, end - begin));
		}

		[[nodiscard]] std::string lower(std::string text)
		{
			std::transform(text.begin(), text.end(), text.begin(),
			               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
			return text;
		}

		[[nodiscard]] bool toDouble(const std::string& text, double& out) noexcept
		{
			if (text.empty())
				return false;
			char* end = nullptr;
			const double value = std::strtod(text.c_str(), &end);
			if (end == text.c_str() || !std::isfinite(value))
				return false;
			out = value;
			return true;
		}

		[[nodiscard]] int decimalPlaces(const std::string& text) noexcept
		{
			const std::size_t dot = text.find('.');
			if (dot == std::string::npos)
				return 0;
			std::size_t end = dot + 1;
			while (end < text.size() && std::isdigit(static_cast<unsigned char>(text[end])))
				++end;
			return static_cast<int>(end - dot - 1);
		}

		[[nodiscard]] int gregorianMidnightMjd(int year, int month, int day) noexcept
		{
			const int a = (14 - month) / 12;
			const int y = year + 4800 - a;
			const int m = month + 12 * a - 3;
			const long long jdn = day + (153LL * m + 2) / 5 + 365LL * y + y / 4 - y / 100 + y / 400 - 32045;
			return static_cast<int>(jdn - 2400001);
		}

		[[nodiscard]] bool parseObsTime(const std::string& text, std::string_view submissionFormat,
		                                obs::ParsedObservation& data) noexcept
		{
			if (text.size() < 19)
				return false;
			const auto num = [&text](std::size_t pos, std::size_t len) {
				return std::atoi(text.substr(pos, len).c_str());
			};
			if (text[4] != '-' || text[7] != '-' || text[10] != 'T' || text[13] != ':' || text[16] != ':')
				return false;

			const int year = num(0, 4);
			const int month = num(5, 2);
			const int day = num(8, 2);
			const int hour = num(11, 2);
			const int minute = num(14, 2);

			double seconds = 0.0;
			std::string secText = text.substr(17);
			if (!secText.empty() && (secText.back() == 'Z' || secText.back() == 'z'))
				secText.pop_back();
			if (!toDouble(secText, seconds))
				return false;

			{
				const std::size_t dot = secText.find('.');
				const int decimals = dot == std::string::npos ? 0 : static_cast<int>(secText.size() - dot - 1);
				data.timeAccuracyDays = std::pow(10.0, -decimals) / 86400.0;
			}

			data.utcMjdDay = gregorianMidnightMjd(year, month, day);
			data.utcSecondsOfDay = hour * 3600.0 + minute * 60.0 + seconds;
			data.hasSplitUtc = true;
			data.jd = 2400000.5 + (static_cast<double>(data.utcMjdDay) + data.utcSecondsOfDay / 86400.0);
			data.timeScale = submissionFormat == "PRE" || submissionFormat == "M47" || submissionFormat == "M92"
			                     ? time::legacyObservationScale(data.jd)
			                     : time::Scale::UTC;
			return true;
		}

		using FieldMap = std::unordered_map<std::string, std::string>;

		[[nodiscard]] bool decodeProgram(std::string_view text, char& program) noexcept
		{
			constexpr std::string_view base62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
			constexpr std::string_view codes =
			    R"(0123456789!"#$%&'()*+,-./[\]^_`{|}~:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz)";
			if (text.size() != 2)
				return false;
			const auto high = base62.find(text[0]);
			const auto low = base62.find(text[1]);
			if (high == base62.npos || low == base62.npos)
				return false;
			const auto index = high * 62 + low;
			program = index < codes.size() ? codes[index] : ' ';
			return true;
		}

		[[nodiscard]] const std::string* find(const FieldMap& fields, const char* key)
		{
			const auto it = fields.find(key);
			return it == fields.end() || it->second.empty() ? nullptr : &it->second;
		}

		[[nodiscard]] bool isInterstellarPermanentId(std::string_view id) noexcept
		{
			return id.size() >= 2 && id.front() >= '1' && id.front() <= '9' && id.back() == 'I' &&
			       std::all_of(id.begin(), id.end() - 1, [](char c) { return c >= '0' && c <= '9'; });
		}

		[[nodiscard]] bool buildRecord(const FieldMap& fields, obs::ParsedObservation& data, std::string& error)
		{
			const auto* permanentId = find(fields, "permID");
			if (permanentId)
				data.permanentDesignation = *permanentId;
			for (const char* key : {"permID", "provID"}) {
				if (const auto* designation = find(fields, key)) {
					data.classificationDesignation = *designation;
					break;
				}
			}
			const auto* format = find(fields, "subFmt");
			const std::string_view submissionFormat = format ? std::string_view(*format) : std::string_view{};
			data.isInterstellar = permanentId && isInterstellarPermanentId(*permanentId);
			if (fields.contains("trx") || fields.contains("rcv")) {
				error = "ADES radar requires obsTime, trx, rcv, com, frq and a measurement with positive RMS";
				const auto* time = find(fields, "obsTime");
				const auto* tx = find(fields, "trx");
				const auto* rx = find(fields, "rcv");
				const auto* com = find(fields, "com");
				const auto* frequency = find(fields, "frq");
				if (!time || !tx || !rx || !com || !frequency || !parseObsTime(*time, submissionFormat, data))
					return false;
				obs::RadarObservation radar;
				radar.transmitter = *tx;
				if (*com != "0" && *com != "1")
					return false;
				radar.returnType = *com == "1" ? obs::RadarReturn::CenterOfMass : obs::RadarReturn::Surface;
				if (!toDouble(*frequency, radar.frequencyHz))
					return false;
				radar.frequencyHz *= 1e6;
				const auto read = [&](const char* valueKey, const char* sigmaKey, double sigmaScale,
				                      std::optional<obs::RadarMeasurement>& channel) {
					const auto* value = find(fields, valueKey);
					const auto* sigma = find(fields, sigmaKey);
					if (!value && !sigma)
						return true;
					obs::RadarMeasurement measurement;
					if (!value || !sigma || !toDouble(*value, measurement.value) ||
					    !toDouble(*sigma, measurement.sigma))
						return false;
					measurement.sigma *= sigmaScale;
					channel = measurement;
					return true;
				};
				if (!read("delay", "rmsDelay", 1e-6, radar.delaySeconds) ||
				    !read("doppler", "rmsDoppler", 1.0, radar.dopplerHz) || !radar.valid())
					return false;
				for (const char* key : {"trkSub", "provID", "permID"})
					if (const auto* name = find(fields, key)) {
						data.objectName = *name;
						break;
					}
				if (data.objectName.empty())
					return false;
				data.obsCode = *rx;
				data.sourceFormat = "ADES";
				data.radar = std::move(radar);
				error.clear();
				return true;
			}
			const std::string* obsTime = find(fields, "obsTime");
			const std::string* raText = find(fields, "ra");
			const std::string* decText = find(fields, "dec");
			const std::string* raStarText = find(fields, "raStar");
			const std::string* decStarText = find(fields, "decStar");
			const bool occultation = raStarText || decStarText;
			if (occultation && (raText || decText)) {
				error = "ADES record mixes ra/dec with occultation raStar/decStar";
				return false;
			}
			if (occultation ? !(obsTime && raStarText && decStarText) : !(obsTime && raText && decText)) {
				error = occultation ? "ADES occultation missing obsTime/raStar/decStar"
				                    : "ADES record missing obsTime/ra/dec";
				return false;
			}
			if (!parseObsTime(*obsTime, submissionFormat, data)) {
				error = "ADES obsTime not ISO 8601 UTC: " + *obsTime;
				return false;
			}
			double raDeg = 0.0;
			double decDeg = 0.0;
			std::optional<double> offsetRaAccuracy;
			std::optional<double> offsetDecAccuracy;
			if (occultation) {
				const std::string* deltaRaText = find(fields, "deltaRA");
				const std::string* deltaDecText = find(fields, "deltaDec");
				if (!deltaRaText || !deltaDecText) {
					error = "ADES occultation requires deltaRA/deltaDec (dist/pa offsets are not supported)";
					return false;
				}
				double starRa = 0.0, starDec = 0.0, deltaRa = 0.0, deltaDec = 0.0;
				if (!toDouble(*raStarText, starRa) || !toDouble(*decStarText, starDec) ||
				    !toDouble(*deltaRaText, deltaRa) || !toDouble(*deltaDecText, deltaDec)) {
					error = "ADES occultation raStar/decStar/deltaRA/deltaDec not numeric";
					return false;
				}
				decDeg = starDec + deltaDec / 3600.0;
				raDeg = starRa + deltaRa / 3600.0 / std::cos(decDeg * Constants::DEG_TO_RAD);
				raDeg -= 360.0 * std::floor(raDeg / 360.0);
				offsetRaAccuracy = std::pow(10.0, -decimalPlaces(*deltaRaText));
				offsetDecAccuracy = std::pow(10.0, -decimalPlaces(*deltaDecText));
			} else if (!toDouble(*raText, raDeg) || !toDouble(*decText, decDeg)) {
				error = "ADES ra/dec not decimal degrees";
				return false;
			}
			data.raRad = raDeg * Constants::DEG_TO_RAD;
			data.decRad = decDeg * Constants::DEG_TO_RAD;

			for (const char* key : {"trkSub", "provID", "permID", "artSat"}) {
				if (const std::string* v = find(fields, key)) {
					data.objectName = *v;
					break;
				}
			}

			if (const std::string* v = find(fields, "stn"))
				data.obsCode = *v;
			if (const std::string* v = find(fields, "mode"))
				data.technique = techniqueFromMode(*v);
			if (const std::string* v = find(fields, "notes"))
				data.notes = *v;
			if (const std::string* v = find(fields, "prog")) {
				data.adesProgramCode = *v;
				if (!decodeProgram(*v, data.programCode)) {
					error = "ADES prog must be a two-character base-62 code: " + *v;
					return false;
				}
			}
			if (data.technique == obs::Technique::Unknown && data.notes.find('n') != std::string::npos)
				data.technique = obs::Technique::NormalPlace;
			if (occultation)
				data.technique = obs::Technique::Occultation;
			if (const std::string* v = find(fields, "deprecated"); v && *v == "X")
				data.technique = obs::Technique::ReplacedCcd;
			if (const std::string* v = find(fields, "astCat")) {
				data.astCat = *v;
				if (const auto* catalog = obs::findStarCatalog(*v))
					data.catalogCode = catalog->mpcCode;
			}
			if (const std::string* v = find(fields, "mag"))
				data.magnitudeText = *v;
			if (const std::string* v = find(fields, "band"))
				data.magBand = photometricBand(*v);
			if (const std::string* v = find(fields, "nucMag")) {
				if (*v == "1" || *v == "true" || *v == "True")
					data.nuclearMagnitude = true;
			}
			data.sourceFormat = "ADES";
			if (const std::string* v = find(fields, "subFmt")) {
				data.sourceFormat += " " + *v;
				if (*v == "M92" && (data.magBand == 'N' || data.magBand == 'T')) {
					const bool nuclear = data.magBand == 'N';
					const auto* flag = find(fields, "nucMag");
					const bool explicitlyNonNuclear = flag && (*flag == "0" || *flag == "false" || *flag == "False");
					if ((data.nuclearMagnitude && *data.nuclearMagnitude != nuclear) ||
					    (nuclear && explicitlyNonNuclear))
						data.nuclearMagnitude.reset();
					else
						data.nuclearMagnitude = nuclear;
				}
			}

			double value = 0.0;
			if (const std::string* v = find(fields, "rmsMag");
			    v && toDouble(*v, value) && std::isfinite(value) && value > 0.0)
				data.rmsMag = value;
			if (const std::string* v = find(fields, "rmsRA"); v && toDouble(*v, value)) {
				data.rmsRaArcsec = value;
				data.hasObserverSigma = true;
			}
			if (const std::string* v = find(fields, "rmsDec"); v && toDouble(*v, value)) {
				data.rmsDecArcsec = value;
				data.hasObserverSigma = true;
			}
			if (const std::string* v = find(fields, "rmsCorr"); v && toDouble(*v, value))
				data.rmsCorr = value;
			if (const std::string* v = find(fields, "rmsTime"); v && toDouble(*v, value)) {
				data.rmsTimeSec = value;
				data.hasRmsTime = true;
			}
			if (const std::string* v = find(fields, "uncTime"); v && toDouble(*v, value)) {
				data.uncTimeSec = value;
				data.hasUncTime = true;
			}

			if (const std::string* v = find(fields, "precRA"); v && toDouble(*v, value)) {
				const double cosDec = std::cos(data.decRad);
				data.raAccuracyArcsec = value * 15.0 * std::max(std::abs(cosDec), 1.0e-12);
			}
			if (const std::string* v = find(fields, "precDec"); v && toDouble(*v, value))
				data.decAccuracyArcsec = value;
			if (const std::string* v = find(fields, "precTime"); v && toDouble(*v, value)) {
				data.timeAccuracyDays = value * 1.0e-6;
				data.hasTimeAccuracy = true;
			}
			if (occultation) {
				const bool precGiven = find(fields, "precRA") || find(fields, "precDec");
				data.raAccuracyArcsec =
				    precGiven ? std::max(data.raAccuracyArcsec, *offsetRaAccuracy) : *offsetRaAccuracy;
				data.decAccuracyArcsec =
				    precGiven ? std::max(data.decAccuracyArcsec, *offsetDecAccuracy) : *offsetDecAccuracy;
				if (data.hasObserverSigma && data.raAccuracyArcsec > data.rmsRaArcsec &&
				    data.decAccuracyArcsec > data.rmsDecArcsec) {
					data.rmsRaArcsec = data.raAccuracyArcsec;
					data.rmsDecArcsec = data.decAccuracyArcsec;
					data.rmsCorr = 0.0;
				}
			}

			if (const std::string* sys = find(fields, "sys")) {
				const std::string tag = lower(*sys);
				double p1 = 0.0, p2 = 0.0, p3 = 0.0;
				const std::string* v1 = find(fields, "pos1");
				const std::string* v2 = find(fields, "pos2");
				const std::string* v3 = find(fields, "pos3");
				const auto* center = find(fields, "ctr");
				if (!v1 || !v2 || !v3 || !toDouble(*v1, p1) || !toDouble(*v2, p2) || !toDouble(*v3, p3) ||
				    !std::isfinite(p1) || !std::isfinite(p2) || !std::isfinite(p3) || (center && *center != "399")) {
					error = "Invalid or non-geocentric ADES observer position";
					return false;
				}
				if (tag == "wgs84") {
					data.rovingObserver = obs::RovingObserver{p1, p2, p3};
					if (!data.rovingObserver->valid()) {
						error = "Invalid ADES WGS84 observer position";
						return false;
					}
				} else if (tag == "icrf_km" || tag == "icrf_au") {
					data.isSatellite = true;
					data.offsetUnits = tag == "icrf_au" ? obs::OffsetUnits::Au : obs::OffsetUnits::Km;
					data.offsetX = p1;
					data.offsetY = p2;
					data.offsetZ = p3;
				} else {
					error = "Unsupported ADES observer coordinate system: " + *sys;
					return false;
				}
			}

			return true;
		}

		[[nodiscard]] std::string decodeEntities(std::string text)
		{
			static const std::array<std::pair<const char*, char>, 5> entities = {{
			    {"&lt;", '<'},
			    {"&gt;", '>'},
			    {"&quot;", '"'},
			    {"&apos;", '\''},
			    {"&amp;", '&'},
			}};
			for (const auto& [name, ch] : entities) {
				std::size_t pos = 0;
				const std::size_t len = std::char_traits<char>::length(name);
				while ((pos = text.find(name, pos)) != std::string::npos) {
					text.replace(pos, len, 1, ch);
					++pos;
				}
			}
			return text;
		}

		[[nodiscard]] std::vector<obs::ParsedObservation> parseXml(const std::string& text, std::string& error)
		{
			std::vector<obs::ParsedObservation> out;
			std::size_t pos = 0;
			auto optical = text.find("<optical>", pos);
			auto radar = text.find("<radar>", pos);
			auto occultation = text.find("<occultation>", pos);
			while (true) {
				pos = std::min({optical, radar, occultation});
				if (pos == std::string::npos)
					break;
				const std::string element = pos == radar ? "radar" : pos == occultation ? "occultation" : "optical";
				const std::string openTag = "<" + element + ">";
				const std::string closeTag = "</" + element + ">";
				const std::size_t bodyBegin = pos + openTag.size();
				const std::size_t bodyEnd = text.find(closeTag, bodyBegin);
				if (bodyEnd == std::string::npos) {
					error = "ADES XML: unterminated " + openTag + " block";
					return {};
				}

				const std::string_view body(text.data() + bodyBegin, bodyEnd - bodyBegin);
				FieldMap fields;
				std::size_t cursor = 0;
				while (cursor < body.size()) {
					const std::size_t tagBegin = body.find('<', cursor);
					if (tagBegin == std::string::npos)
						break;
					const std::size_t tagEnd = body.find('>', tagBegin);
					if (tagEnd == std::string::npos)
						break;

					std::string name(body.substr(tagBegin + 1, tagEnd - tagBegin - 1));
					if (!name.empty() && name.back() == '/') {
						cursor = tagEnd + 1;
						continue;
					}
					if (!name.empty() && name.front() == '/') {
						cursor = tagEnd + 1;
						continue;
					}

					const std::string endMarker = "</" + name + ">";
					const std::size_t valueEnd = body.find(endMarker, tagEnd + 1);
					if (valueEnd == std::string::npos) {
						cursor = tagEnd + 1;
						continue;
					}
					fields[name] = decodeEntities(trim(body.substr(tagEnd + 1, valueEnd - tagEnd - 1)));
					cursor = valueEnd + endMarker.size();
				}

				obs::ParsedObservation data;
				if (!buildRecord(fields, data, error))
					return {};
				out.push_back(std::move(data));
				pos = bodyEnd + closeTag.size();
				if (optical < pos)
					optical = text.find("<optical>", pos);
				if (radar < pos)
					radar = text.find("<radar>", pos);
				if (occultation < pos)
					occultation = text.find("<occultation>", pos);
			}

			if (out.empty())
				error = "ADES XML: no optical, occultation or radar records found";
			return out;
		}

		[[nodiscard]] std::vector<std::string> splitPipe(const std::string& line)
		{
			std::vector<std::string> out;
			std::size_t begin = 0;
			while (true) {
				const std::size_t bar = line.find('|', begin);
				if (bar == std::string::npos) {
					out.push_back(trim(std::string_view(line).substr(begin)));
					break;
				}
				out.push_back(trim(std::string_view(line).substr(begin, bar - begin)));
				begin = bar + 1;
			}
			return out;
		}

		[[nodiscard]] bool looksLikeKeywordRecord(const std::vector<std::string>& cells)
		{
			const bool hasTime = std::find(cells.begin(), cells.end(), "obsTime") != cells.end();
			const bool hasStn = std::find(cells.begin(), cells.end(), "stn") != cells.end();
			const bool hasReceiver = std::find(cells.begin(), cells.end(), "rcv") != cells.end();
			return hasTime && (hasStn || hasReceiver);
		}

		[[nodiscard]] std::vector<obs::ParsedObservation> parsePsv(const std::string& text, std::string& error)
		{
			std::vector<obs::ParsedObservation> out;
			std::vector<std::string> columns;
			std::size_t begin = 0;

			while (begin <= text.size()) {
				std::size_t end = text.find('\n', begin);
				if (end == std::string::npos)
					end = text.size();
				std::string line = trim(std::string_view(text).substr(begin, end - begin));
				begin = end + 1;

				if (line.empty() || line.front() == '#')
					continue;
				if (line.find('|') == std::string::npos)
					continue;

				std::vector<std::string> cells = splitPipe(line);
				if (looksLikeKeywordRecord(cells)) {
					columns = std::move(cells);
					continue;
				}
				if (columns.empty())
					continue;

				FieldMap fields;
				const std::size_t n = std::min(columns.size(), cells.size());
				for (std::size_t i = 0; i < n; ++i) {
					if (!cells[i].empty())
						fields[columns[i]] = cells[i];
				}

				obs::ParsedObservation data;
				if (!buildRecord(fields, data, error))
					return {};
				out.push_back(std::move(data));
			}

			if (out.empty())
				error = "ADES PSV: no data records found";
			return out;
		}

	}  // namespace

	const char* AdesReader::formatName() const noexcept { return "ADES"; }

	bool AdesReader::canRead(std::string_view text) const noexcept
	{
		const std::size_t adesRoot = text.find("<ades");
		if (adesRoot != std::string_view::npos && adesRoot < 4096)
			return true;

		const std::size_t scanEnd = std::min<std::size_t>(text.size(), 65536);
		std::size_t begin = 0;
		while (begin < scanEnd) {
			std::size_t end = text.find('\n', begin);
			if (end == std::string_view::npos || end > scanEnd)
				end = scanEnd;
			const std::string line = trim(text.substr(begin, end - begin));
			begin = end + 1;
			if (line.empty() || line.front() == '#')
				continue;
			if (line.find('|') != std::string::npos && looksLikeKeywordRecord(splitPipe(line))) {
				return true;
			}
		}
		return false;
	}

	std::vector<obs::ParsedObservation> AdesReader::read(const std::string& text, std::string& error) const
	{
		error.clear();
		const std::size_t adesRoot = text.find("<ades");
		if (adesRoot != std::string::npos && adesRoot < 4096)
			return parseXml(text, error);
		return parsePsv(text, error);
	}

}  // namespace orbital_core::ades
