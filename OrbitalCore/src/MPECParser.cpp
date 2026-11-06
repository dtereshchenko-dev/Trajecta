#include "MPECParser.h"
#include <string>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <string_view>
#include <vector>
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "logger/Logger.h"

namespace {

	[[nodiscard]] std::string_view trimLeft(std::string_view text) noexcept
	{
		while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front()))) {
			text.remove_prefix(1);
		}
		return text;
	}

	[[nodiscard]] bool startsWith(std::string_view text, std::string_view prefix) noexcept
	{
		return text.size() >= prefix.size() && text.substr(0, prefix.size()) == prefix;
	}

	[[nodiscard]] bool isInterstellarDirective(std::string_view line) noexcept
	{
		line = trimLeft(line);
		return startsWith(line, "#interstellar") || startsWith(line, "COM interstellar");
	}

	[[nodiscard]] bool isInterstellarPackedId(std::string_view packedId) noexcept
	{
		return packedId.size() >= 5 && packedId[0] == '0' && packedId[1] == '0' && packedId[4] == 'I';
	}

	[[nodiscard]] std::vector<std::string> splitFields(std::string_view text)
	{
		std::istringstream in{std::string(text)};
		std::vector<std::string> tokens;
		std::string token;
		while (in >> token)
			tokens.push_back(token);
		return tokens;
	}

	[[nodiscard]] std::string formatRa(int hours, int minutes, double seconds)
	{
		std::ostringstream out;
		out << std::setfill('0') << std::setw(2) << hours << ' ' << std::setw(2) << minutes << ' ' << std::setfill(' ')
		    << std::setprecision(12) << seconds;
		return out.str();
	}

	[[nodiscard]] std::string formatDec(bool negative, int degrees, int minutes, double seconds)
	{
		std::ostringstream out;
		out << (negative ? '-' : '+') << std::setfill('0') << std::setw(2) << degrees << ' ' << std::setw(2) << minutes
		    << ' ' << std::setfill(' ') << std::setprecision(12) << seconds;
		return out.str();
	}

	[[nodiscard]] int fractionalDigits(std::string_view token) noexcept
	{
		const std::size_t dot = token.find('.');
		if (dot == std::string_view::npos)
			return 0;
		int digits = 0;
		for (std::size_t i = dot + 1; i < token.size(); ++i) {
			if (!std::isdigit(static_cast<unsigned char>(token[i])))
				break;
			++digits;
		}
		return digits;
	}

	[[nodiscard]] double decimalResolution(std::string_view token) noexcept
	{
		const int digits = fractionalDigits(token);
		return digits > 0 ? std::pow(10.0, -digits) : 1.0;
	}

	[[nodiscard]] int gregorianMidnightMjd(int year, int month, int day) noexcept
	{
		return static_cast<int>(
		    std::llround(Convertors::convertGregorianToJD(year, month, static_cast<double>(day)) - 2400000.5));
	}

	[[nodiscard]] bool parseMpcDate(int year, int month, std::string_view dayToken, MPECInitData& data) noexcept
	{
		const std::size_t point = dayToken.find('.');

		int day = 0;
		if (!Convertors::tryParseNumber(dayToken.substr(0, point), day))
			return false;

		double fractionalDay = 0.0;
		if (point != std::string_view::npos && !Convertors::tryParseNumber(dayToken.substr(point), fractionalDay)) {
			return false;
		}

		data.utcMjdDay = gregorianMidnightMjd(year, month, day);
		data.utcSecondsOfDay = fractionalDay * 86400.0;
		data.hasSplitUtc = true;
		data.JD = 2400000.5 + (static_cast<double>(data.utcMjdDay) + data.utcSecondsOfDay / 86400.0);
		return true;
	}

	[[nodiscard]] bool parseMpc80RaField(std::string_view field, std::string& outRa)
	{
		const auto tokens = splitFields(field);
		if (tokens.size() != 2 && tokens.size() != 3)
			return false;

		int hours = 0;
		if (!Convertors::tryParseNumber(tokens[0], hours) || hours < 0 || hours > 23)
			return false;

		int minutes = 0;
		double seconds = 0.0;
		if (tokens.size() == 3) {
			if (!Convertors::tryParseNumber(tokens[1], minutes) || !Convertors::tryParseNumber(tokens[2], seconds)) {
				return false;
			}
		}

		else {
			double decimalMinutes = 0.0;
			if (!Convertors::tryParseNumber(tokens[1], decimalMinutes))
				return false;

			minutes = static_cast<int>(std::floor(decimalMinutes));
			seconds = (decimalMinutes - static_cast<double>(minutes)) * 60.0;
		}

		if (minutes < 0 || minutes >= 60 || seconds < 0.0 || seconds >= 60.0)
			return false;

		outRa = formatRa(hours, minutes, seconds);
		return true;
	}

	[[nodiscard]] double mpc80RaAccuracyArcsec(std::string_view field)
	{
		const auto tokens = splitFields(field);
		if (tokens.size() == 3)
			return 15.0 * decimalResolution(tokens[2]);
		if (tokens.size() == 2)
			return 15.0 * 60.0 * decimalResolution(tokens[1]);
		return 1.5e-2;
	}

	[[nodiscard]] bool parseMpc80DecField(std::string_view field, std::string& outDec)
	{
		const auto tokens = splitFields(field);
		if (tokens.size() != 2 && tokens.size() != 3)
			return false;

		const bool negative = !tokens[0].empty() && tokens[0].front() == '-';

		int signedDegrees = 0;
		if (!Convertors::tryParseNumber(tokens[0], signedDegrees))
			return false;

		const int degrees = std::abs(signedDegrees);

		int minutes = 0;
		double seconds = 0.0;
		if (tokens.size() == 3) {
			if (!Convertors::tryParseNumber(tokens[1], minutes) || !Convertors::tryParseNumber(tokens[2], seconds)) {
				return false;
			}
		}

		else {
			double decimalMinutes = 0.0;
			if (!Convertors::tryParseNumber(tokens[1], decimalMinutes))
				return false;

			minutes = static_cast<int>(std::floor(decimalMinutes));
			seconds = (decimalMinutes - static_cast<double>(minutes)) * 60.0;
		}

		if (degrees > 90 || minutes < 0 || minutes >= 60 || seconds < 0.0 || seconds >= 60.0) {
			return false;
		}

		outDec = formatDec(negative, degrees, minutes, seconds);
		return true;
	}

	[[nodiscard]] double mpc80DecAccuracyArcsec(std::string_view field)
	{
		const auto tokens = splitFields(field);
		if (tokens.size() == 3)
			return decimalResolution(tokens[2]);
		if (tokens.size() == 2)
			return 60.0 * decimalResolution(tokens[1]);
		return 1.0e-2;
	}

}  // namespace

MPECParser::MPECParser(const std::vector<std::string>& lines) : m_lines{lines} {}

std::vector<MPECInitData> MPECParser::parse() noexcept
{
	if (m_lines.empty()) {
		m_error = Error::InvalidData;
		return {};
	}

	std::vector<MPECInitData> data;
	std::size_t satCount = 0;
	std::size_t orphanContinuations = 0;
	std::size_t malformedContinuations = 0;
	bool interstellarObject = false;

	for (std::size_t i = 0; i < m_lines.size(); ++i) {
		const std::string& line = m_lines[i];

		if (isInterstellarDirective(line)) {
			interstellarObject = true;
			continue;
		}

		if (line.size() > 14 && (line[14] == 'R' || line[14] == 'Q')) {
			MPECInitData radar;
			if (i + 1 >= m_lines.size() || !parseRadarPair(line, m_lines[i + 1], radar)) {
				m_error = Error::InvalidFormat;
				return {};
			}
			data.push_back(std::move(radar));
			++i;
			continue;
		}
		if (line.size() > 14 && (line[14] == 'r' || line[14] == 'q')) {
			m_error = Error::InvalidFormat;
			return {};
		}

		if (isSatelliteContinuation(line) || isRovingContinuation(line)) {
			++orphanContinuations;
			continue;
		}

		MPECInitData obj = parseLine(line);
		if (m_error != Error::NoError) {
			LOG_WARN << "Skipping malformed line " << (i + 1) << ": " << errorString();
			m_error = Error::NoError;
			continue;
		}

		if (obj.isInterstellar)
			interstellarObject = true;

		if (line.size() >= 80 && (obj.note2 == 'S' || obj.note2 == 'V')) {
			const bool satellite = obj.note2 == 'S';
			bool valid = i + 1 < m_lines.size();
			if (valid) {
				const auto& next = m_lines[i + 1];
				valid = next.size() >= 80 && (satellite ? isSatelliteContinuation(next) : isRovingContinuation(next)) &&
				        next.substr(0, 12) == line.substr(0, 12) && next.substr(77, 3) == line.substr(77, 3);
				if (valid) {
					int year = 0, month = 0;
					std::string day;
					std::istringstream date(next.substr(15, 17));
					MPECInitData stamp;
					valid = static_cast<bool>(date >> year >> month >> day) && parseMpcDate(year, month, day, stamp) &&
					        stamp.utcMjdDay == obj.utcMjdDay &&
					        std::abs(stamp.utcSecondsOfDay - obj.utcSecondsOfDay) < 1.0e-6;
				}
				if (valid)
					valid = satellite ? parseSatelliteContinuation(next, obj) : parseRovingContinuation(next, obj);
			}
			if (!valid) {
				++malformedContinuations;
				LOG_WARN << "Skipping observation at line " << (i + 1)
				         << ": missing or mismatched observer continuation";
				continue;
			}
			obj.isSatellite = satellite;
			++satCount;
			++i;
		}

		data.push_back(obj);
	}

	if (interstellarObject) {
		for (auto& obs : data)
			obs.isInterstellar = true;
		LOG_INFO << "MPECParser: interstellar object flag enabled";
	}

	if (satCount > 0 || orphanContinuations > 0 || malformedContinuations > 0) {
		LOG_INFO << "MPECParser: parsed " << data.size() << " obs records  ("
		         << "ground=" << (data.size() - satCount) << ", satellite=" << satCount
		         << ", orphan-s-lines=" << orphanContinuations << ", malformed-continuations=" << malformedContinuations
		         << ")";
	}

	m_error = Error::NoError;
	return data;
}

bool MPECParser::isMPC80Col(const std::string& line) noexcept
{
	if (line.size() > 14 && (line[14] == 's' || line[14] == 'v'))
		return false;
	return line.size() >= 80 && (line[44] == '+' || line[44] == '-');
}

bool MPECParser::parseRadarPair(const std::string& first, const std::string& second, MPECInitData& data)
{
	if (first.size() < 80 || second.size() < 80 || second[14] != static_cast<char>(std::tolower(first[14])) ||
	    first.substr(0, 12) != second.substr(0, 12) || first.substr(15, 17) != second.substr(15, 17) ||
	    first.substr(68, 3) != second.substr(68, 3) || first.substr(77, 3) != second.substr(77, 3) ||
	    (second[32] != 'C' && second[32] != 'S'))
		return false;
	const auto implied = [](std::string_view field, std::size_t decimal, double& value) {
		if (field.find_first_not_of(' ') == std::string_view::npos)
			return false;
		std::string number;
		for (std::size_t i = 0; i < field.size(); ++i) {
			if (i == decimal)
				number += '.';
			const char c = field[i];
			if (c == '+' || c == '-') {
				if (i != 0)
					return false;
				number += c;
			} else if (c == ' ' || (c >= '0' && c <= '9')) {
				number += c == ' ' ? '0' : c;
			} else
				return false;
		}
		return Convertors::tryParseNumber(number, value) && std::isfinite(value);
	};
	int year = 0, month = 0;
	std::string day;
	std::istringstream date(first.substr(15, 17));
	double decimalDay = 0.0;
	if (!(date >> year >> month >> day) || year < 1600 || month < 1 || month > 12 ||
	    !Convertors::tryParseNumber(day, decimalDay) || decimalDay < 1.0 || decimalDay >= 32.0 ||
	    !parseMpcDate(year, month, day, data))
		return false;
	data.utcSecondsOfDay = std::round(data.utcSecondsOfDay);
	if (data.utcSecondsOfDay >= 86400.0) {
		++data.utcMjdDay;
		data.utcSecondsOfDay = 0.0;
	}
	data.JD = 2400000.5 + (static_cast<double>(data.utcMjdDay) + data.utcSecondsOfDay / 86400.0);
	data.objectName = first.substr(0, 5);
	data.classificationDesignation = first.substr(0, 12);
	if (data.objectName.find_first_not_of(' ') == std::string::npos)
		data.objectName = first.substr(5, 7);
	data.objectName.erase(std::remove(data.objectName.begin(), data.objectName.end(), ' '), data.objectName.end());
	if (data.objectName.empty())
		return false;
	data.code = first.substr(77, 3);
	data.note2 = 'R';
	data.timeAccuracyDays = 1.0 / 86400.0;
	orbital_core::obs::RadarObservation radar;
	radar.transmitter = first.substr(68, 3);
	radar.returnType =
	    second[32] == 'C' ? orbital_core::obs::RadarReturn::CenterOfMass : orbital_core::obs::RadarReturn::Surface;
	double frequency = 0.0;
	if (!implied(first.substr(62, 6) + second.substr(62, 6), 5, frequency))
		return false;
	radar.frequencyHz = frequency * 1.0e6;
	if (first.substr(32, 15).find_first_not_of(' ') != std::string::npos) {
		orbital_core::obs::RadarMeasurement delay;
		if (!implied(first.substr(32, 15), 5, delay.value) || !implied(second.substr(33, 14), 10, delay.sigma))
			return false;
		delay.sigma *= 1.0e-6;
		radar.delaySeconds = delay;
	}
	if (first.substr(47, 15).find_first_not_of(' ') != std::string::npos) {
		orbital_core::obs::RadarMeasurement doppler;
		if (!implied(first.substr(47, 15), 11, doppler.value) || !implied(second.substr(47, 15), 11, doppler.sigma))
			return false;
		radar.dopplerHz = doppler;
	}
	if (!radar.valid())
		return false;
	data.radar = std::move(radar);
	return true;
}

bool MPECParser::isSatelliteContinuation(const std::string& line) noexcept
{
	return line.size() > 14 && line[14] == 's';
}

bool MPECParser::isRovingContinuation(const std::string& line) noexcept { return line.size() > 14 && line[14] == 'v'; }

bool MPECParser::parseSatelliteContinuation(const std::string& line, MPECInitData& into) noexcept
{
	if (line.size() < 80 || (line[32] != '1' && line[32] != '2'))
		return false;

	const std::string payload = line.substr(33, std::min<std::size_t>(line.size() - 33, 44));

	std::istringstream ss(payload);
	std::vector<std::string> tokens;
	std::string tok;
	while (ss >> tok)
		tokens.push_back(tok);

	auto parseSigned = [](const std::string& s) -> std::pair<bool, double> {
		double value = 0.0;
		std::size_t consumed = 0;
		const bool ok = Convertors::tryParseNumber(s, value, &consumed) && consumed == s.size() && std::isfinite(value);
		return {ok, ok ? value : 0.0};
	};

	std::vector<double> values;
	values.reserve(3);
	for (std::size_t i = 0; i < tokens.size() && values.size() < 3; ++i) {
		if ((tokens[i] == "+" || tokens[i] == "-") && i + 1 < tokens.size()) {
			const auto [ok, v] = parseSigned(tokens[++i]);
			if (!ok)
				return false;
			values.push_back(tokens[i - 1] == "-" ? -v : v);
		} else {
			const auto [ok, v] = parseSigned(tokens[i]);
			if (!ok)
				return false;
			values.push_back(v);
		}
	}

	if (values.size() < 3)
		return false;

	const double x = values[0];
	const double y = values[1];
	const double z = values[2];

	into.satOffsetUnits = line[32];
	into.satOffsetX = x;
	into.satOffsetY = y;
	into.satOffsetZ = z;
	return true;
}

bool MPECParser::parseRovingContinuation(const std::string& line, MPECInitData& into) noexcept
{
	if (line.size() < 80)
		return false;

	const std::string payload = line.substr(33, std::min<std::size_t>(line.size() - 33, 44));
	std::istringstream ss(payload);
	std::vector<std::string> tokens;
	std::string tok;
	while (ss >> tok)
		tokens.push_back(tok);

	if (tokens.size() < 3)
		return false;

	double lonDeg = 0.0;
	double latDeg = 0.0;
	double altMeters = 0.0;
	if (!Convertors::tryParseNumber(tokens[0], lonDeg) || !Convertors::tryParseNumber(tokens[1], latDeg) ||
	    !Convertors::tryParseNumber(tokens[2], altMeters)) {
		return false;
	}

	const orbital_core::obs::RovingObserver position{lonDeg, latDeg, altMeters};
	if (!position.valid())
		return false;
	into.rovingObserver = position;
	return true;
}

MPECInitData MPECParser::parseLine(const std::string& line)
{
	if (isMPC80Col(line))
		return parseLine80Col(line);
	return parseLineTokenised(line);
}

MPECInitData MPECParser::parseLine80Col(const std::string& line)
{
	MPECInitData data;

	auto trim = [](std::string s) {
		const auto a = s.find_first_not_of(' ');
		if (a == std::string::npos)
			return std::string{};
		const auto b = s.find_last_not_of(' ');
		return s.substr(a, b - a + 1);
	};

	data.objectName = trim(line.substr(0, 5));
	data.classificationDesignation = line.substr(0, 12);
	if (data.objectName.size() < 3)
		data.objectName = trim(line.substr(5, 7));
	if (data.objectName.size() < 3) {
		m_error = Error::InvalidObjectName;
		return {};
	}
	data.isInterstellar = isInterstellarPackedId(data.objectName);

	data.note1 = line[13];
	data.note2 = line[14];

	{
		const std::string dateField = line.substr(15, 17);
		std::istringstream iss(dateField);
		int year = 0, month = 0;
		std::string dayToken;
		iss >> year >> month >> dayToken;
		if (!iss || year < 1800 || year > 2200 || month < 1 || month > 12) {
			m_error = Error::InvalidDate;
			return {};
		}
		data.timeAccuracyDays = decimalResolution(dayToken);
		if (!parseMpcDate(year, month, dayToken, data)) {
			m_error = Error::InvalidDate;
			return {};
		}
	}

	{
		const std::string raField = line.substr(32, 12);
		if (!parseMpc80RaField(raField, data.ra)) {
			m_error = Error::InvalidRAData;
			return {};
		}
		data.raAccuracyArcsec = mpc80RaAccuracyArcsec(raField);
	}

	{
		const std::string decField = line.substr(44, 12);
		if (!parseMpc80DecField(decField, data.dec)) {
			m_error = Error::InvalidDecData;
			return {};
		}
		data.decAccuracyArcsec = mpc80DecAccuracyArcsec(decField);
	}

	{
		const std::string magField = trim(line.substr(65, 5));
		data.magnitude = magField;
		const char band = line[70];
		data.magBand = (band != ' ' && band != '\0') ? std::string{band} : std::string{};
	}

	{
		const char cat = line[71];
		data.catalog = (cat != ' ' && cat != '\0') ? std::string{cat} : std::string{};
		data.filter = data.catalog;
	}

	data.code = trim(line.substr(77, 3));
	if (data.code.size() != 3) {
		m_error = Error::InvalidCode;
		return {};
	}

	return data;
}

MPECInitData MPECParser::parseLineTokenised(const std::string& line)
{
	std::istringstream lineStream(line);

	std::vector<std::string> tokens;
	std::string token;
	while (lineStream >> token) {
		tokens.push_back(token);
	}

	if (tokens.size() < MIN_TOKENS_COUNT) {
		m_error = Error::InsufficientData;
		return {};
	}

	if (tokens[0].empty() || tokens[0].size() < 3) {
		m_error = Error::InvalidObjectName;
		return {};
	}

	if (!isNumber(tokens[4]) || !isNumber(tokens[5]) || !isDouble(tokens[6])) {
		m_error = Error::InvalidRAData;
		return {};
	}

	if (tokens[7][0] != DEC_SIGN && tokens[7][0] != DEC_SIGN_NEGATIVE && !isNumber(tokens[7].substr(1))) {
		m_error = Error::InvalidDecData;
		return {};
	}

	if (!isNumber(tokens[8]) || !isDouble(tokens[9])) {
		m_error = Error::InvalidDecData;
		return {};
	}

	if (tokens.size() <= 10 || !isDouble(tokens[10])) {
		m_error = Error::InvalidMagnitude;
		return {};
	}

	MPECInitData data;

	data.objectName = tokens[0];
	data.classificationDesignation = data.objectName;
	data.isInterstellar = isInterstellarPackedId(data.objectName);
	if (tokens[1].size() < 4) {
		m_error = Error::InvalidDate;
		return {};
	}
	const std::size_t prefixSize = tokens[1].size() - 4;
	if (prefixSize >= 1)
		data.note2 = tokens[1][prefixSize - 1];
	if (prefixSize >= 2)
		data.note1 = tokens[1][prefixSize - 2];
	const std::string year = tokens[1].substr(prefixSize);
	const std::string& month = tokens[2];
	const std::string& day = tokens[3];
	data.timeAccuracyDays = decimalResolution(day);

	if (!parseMpcDate(Convertors::stringToInt(year), Convertors::stringToInt(month), day, data)) {
		m_error = Error::InvalidDate;
		return {};
	}

	data.ra = tokens[4] + RA_DELIMITER + tokens[5] + RA_DELIMITER + tokens[6];
	data.dec = tokens[7] + DEC_DELIMITER + tokens[8] + DEC_DELIMITER + tokens[9];

	std::size_t lastIndex = 10;
	if (tokens.size() > 10 && (isdigit(tokens[10][0]) || tokens[10][0] == '-')) {
		data.magnitude = tokens[10];
		lastIndex = 11;
	}

	else
		data.magnitude = "";

	if (tokens.size() > lastIndex && tokens[lastIndex].size() == 1 && isalpha(tokens[lastIndex][0])) {
		data.magBand = tokens[lastIndex][0];
		data.filter = data.magBand;
		lastIndex++;
	}

	else
		data.filter = "";

	data.code = tokens.back().substr(tokens.back().size() - 3);

	return data;
}

bool MPECParser::isDouble(const std::string& s) noexcept
{
	double value = 0.0;
	return Convertors::tryParseNumber(s, value);
}

bool MPECParser::isNumber(const std::string& s) noexcept
{
	return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}

MPECParser::Error MPECParser::error() const noexcept { return m_error; }

std::string MPECParser::errorString() const noexcept
{
	switch (m_error) {
		case Error::InvalidData:
			return "Invalid data";

		case Error::InsufficientData:
			return "Insufficient data";

		case Error::InvalidRAData:
			return "Invalid RA data";

		case Error::InvalidDecData:
			return "Invalid Dec data";

		case Error::InvalidMagnitude:
			return "Invalid magnitude";

		case Error::InvalidObjectName:
			return "Invalid object name";

		case Error::InvalidDate:
			return "Invalid date conversion";

		case Error::InvalidCode:
			return "Invalid observatory code";

		case Error::InvalidFormat:
			return "Invalid data format";

		default:
			return "No error";
	}
}
