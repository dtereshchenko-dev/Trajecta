#include "orbital_core/obs/ObservatoryUpdates.h"
#include "obs/ObservatoryCatalog.h"
#include "HttpClient.h"

#include <nlohmann/json.hpp>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace orbital_core::obs {
	namespace {
		std::string trim(std::string_view value)
		{
			const auto first = value.find_first_not_of(" \t\r\n");
			if (first == std::string_view::npos)
				return {};
			return std::string(value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1));
		}

		bool validCode(std::string_view code)
		{
			return code.size() == 3 && std::all_of(code.begin(), code.end(), [](char c) {
				       return (c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z');
			       });
		}

		double number(const std::string& text)
		{
			double value = 0.0;
			const char* first = text.data();
			const char* last = first + text.size();
			if (first != last && *first == '+')
				++first;
			const auto [end, error] = std::from_chars(first, last, value);
			if (error != std::errc{} || end != last || !std::isfinite(value))
				throw std::runtime_error("Invalid observatory coordinate");
			return value;
		}

		double jsonNumber(const nlohmann::json& object, const char* field)
		{
			const auto& value = object.at(field);
			const double result = value.is_string() ? number(value.get<std::string>()) : value.get<double>();
			if (!std::isfinite(result))
				throw std::runtime_error("Invalid radar coordinate");
			return result;
		}

		bool sameCatalog(std::string_view left, std::string_view right)
		{
			const auto rows = [](std::string_view text) {
				std::map<std::string, std::string> result;
				std::istringstream input{std::string(text)};
				std::string line;
				while (std::getline(input, line)) {
					std::istringstream row(line);
					std::string code, token;
					if (!(row >> code) || code == "Code" || code[0] == '!')
						continue;
					std::ostringstream normalized;
					normalized.imbue(std::locale::classic());
					while (row >> token) {
						double value = 0;
						const char* first = token.data();
						const char* last = first + token.size();
						if (first != last && *first == '+')
							++first;
						const auto [end, error] = std::from_chars(first, last, value);
						if (value == 0)
							value = 0;
						if (error == std::errc{} && end == last && std::isfinite(value))
							normalized << 'N' << std::hexfloat << value;
						else
							normalized << 'S' << token;
						normalized << '\x1f';
					}
					if (!result.emplace(code, normalized.str()).second)
						throw std::runtime_error("Duplicate observatory code");
				}
				return result;
			};
			try {
				return rows(left) == rows(right);
			} catch (const std::runtime_error&) {
				return false;
			}
		}
	}  // namespace

	void validateOpticalObservatoryCatalog(std::string_view contents)
	{
		std::istringstream input{std::string(contents)};
		std::string line;
		if (!std::getline(input, line))
			throw std::runtime_error("Invalid optical catalog header");
		std::istringstream header(line);
		std::string codeHeader, lonHeader, cosHeader, sinHeader, nameHeader, extraHeader;
		if (!(header >> codeHeader >> lonHeader >> cosHeader >> sinHeader >> nameHeader) || codeHeader != "Code" ||
		    lonHeader != "Long." || cosHeader != "cos" || sinHeader != "sin" || nameHeader != "Name" ||
		    (header >> extraHeader))
			throw std::runtime_error("Invalid optical catalog header");
		std::set<std::string> codes;
		while (std::getline(input, line)) {
			if (trim(line).empty())
				continue;
			std::istringstream row(line);
			std::string code, token;
			if (!(row >> code >> token) || !validCode(code) || !codes.insert(code).second)
				throw std::runtime_error("Invalid optical catalog row");
			if ((token[0] >= '0' && token[0] <= '9') || token[0] == '-' || token[0] == '+' || token[0] == '.') {
				std::string xy, z, name;
				if (!(row >> xy >> z) || !std::getline(row >> std::ws, name) || trim(name).empty())
					throw std::runtime_error("Incomplete optical catalog row");
				const double lon = number(token), x = number(xy), zz = number(z);
				if (lon < 0 || lon >= 360 || x < 0 || std::hypot(x, zz) > 2)
					throw std::runtime_error("Out-of-range optical catalog coordinates");
			}
		}
		if (codes.size() < 1000 || !codes.contains("500"))
			throw std::runtime_error("Incomplete optical catalog");
	}

	ObservatoryUpdate prepareOpticalObservatoryUpdate(std::string_view source, std::string_view current)
	{
		const auto begin = source.find("<pre>");
		const auto end = source.find("</pre>");
		if (begin == std::string_view::npos || end == std::string_view::npos || end <= begin + 5)
			throw std::runtime_error("Invalid MPC observatory response");
		std::map<std::string, std::string> rows;
		std::istringstream previous{std::string(current)};
		std::string line;
		std::getline(previous, line);
		while (std::getline(previous, line)) {
			if (line.size() >= 3 && validCode(std::string_view(line).substr(0, 3)))
				rows.emplace(line.substr(0, 3), trim(line));
		}
		std::istringstream input{std::string(source.substr(begin + 5, end - begin - 5))};
		std::set<std::string> received;
		while (std::getline(input, line)) {
			if (line.size() < 30 || !validCode(std::string_view(line).substr(0, 3)))
				continue;
			const auto code = line.substr(0, 3);
			if (!received.insert(code).second)
				throw std::runtime_error("Duplicate MPC observatory code");
			const auto lon = trim(std::string_view(line).substr(4, 9));
			const auto rhoXY = trim(std::string_view(line).substr(13, 8));
			const auto rhoZ = trim(std::string_view(line).substr(21, 9));
			const auto name = trim(std::string_view(line).substr(30));
			if (name.empty())
				throw std::runtime_error("MPC observatory without a name");
			std::string row = code + " ";
			if (!lon.empty() || !rhoXY.empty() || !rhoZ.empty()) {
				if (lon.empty() || rhoXY.empty() || rhoZ.empty())
					throw std::runtime_error("Incomplete MPC observatory coordinates");
				const double l = number(lon), x = number(rhoXY), z = number(rhoZ);
				if (l < 0 || l >= 360 || x < 0 || std::hypot(x, z) > 2)
					throw std::runtime_error("Out-of-range MPC observatory coordinates");
				row += lon + " " + rhoXY + " " + rhoZ + " ";
			}
			rows[code] = row + name;
		}
		if (received.size() < 1000 || !received.contains("500"))
			throw std::runtime_error("Incomplete MPC observatory catalog");
		std::ostringstream output;
		output << "Code Long. cos sin Name\n";
		for (const auto& [code, row] : rows)
			output << row << '\n';
		const auto result = output.str();
		validateOpticalObservatoryCatalog(result);
		return {result, received.size(), 0, current.empty() || !sameCatalog(current, result)};
	}

	ObservatoryUpdate prepareRadarObservatoryUpdate(std::string_view source, std::string_view current)
	{
		RadarCatalog catalog;
		if (!current.empty()) {
			std::istringstream validation{std::string(current)};
			catalog = readRadarCatalog(validation);
		}
		const auto document = nlohmann::json::parse(source.begin(), source.end());
		if (!document.is_object() || !document.contains("signature") ||
		    document.at("signature").value("version", std::string{}) != "1.1" || !document.contains("coords") ||
		    !document.at("coords").is_object())
			throw std::runtime_error("Invalid JPL radar response");
		const auto& coords = document.at("coords");
		if (coords.size() < 10)
			throw std::runtime_error("Incomplete JPL radar catalog");
		std::set<std::string> used;
		std::size_t withoutMpc = 0;
		std::ostringstream output;
		output.imbue(std::locale::classic());
		output << std::setprecision(17);
		if (current.empty())
			output << "! RADCODE format=2\n";
		std::istringstream input{std::string(current)};
		for (std::string line; std::getline(input, line);) {
			if (line.empty() || line[0] == '!') {
				output << line << '\n';
				continue;
			}
			std::istringstream row(line);
			std::string code, jpl, mpc, lon, xy, z, name;
			if (!(row >> code >> jpl >> mpc >> lon >> xy >> z))
				throw std::runtime_error("Invalid existing radar catalog");
			std::getline(row >> std::ws, name);
			if (jpl != "-" && coords.contains(jpl)) {
				const auto& position = coords.at(jpl);
				const double l = jsonNumber(position, "longitude"), x = jsonNumber(position, "d_xy"),
				             zz = jsonNumber(position, "d_z");
				if (l < 0 || l >= 360 || x < 0 || std::hypot(x, zz) <= 0)
					throw std::runtime_error("Out-of-range JPL radar coordinate");
				output << code << ' ' << jpl << ' ' << mpc << ' ' << l << ' ' << x << ' ' << zz << ' ' << name << '\n';
				used.insert(jpl);
				if (mpc == "-")
					++withoutMpc;
			} else {
				output << line << '\n';
			}
		}
		if (used.empty() && !catalog.stations.empty())
			throw std::runtime_error("JPL radar catalog did not match any known stations");
		std::size_t added = 0;
		for (auto it = coords.begin(); it != coords.end(); ++it) {
			if (used.contains(it.key()))
				continue;
			const auto& code = it.key();
			if (code.size() < 2 || code[0] != '-' || code[1] < '1' || code[1] > '9' ||
			    !std::all_of(code.begin() + 2, code.end(), [](char c) { return c >= '0' && c <= '9'; }))
				throw std::runtime_error("Invalid JPL radar station identifier");
			const auto& position = it.value();
			const auto name = position.at("name").get<std::string>();
			const double l = jsonNumber(position, "longitude"), x = jsonNumber(position, "d_xy"),
			             z = jsonNumber(position, "d_z");
			if (trim(name).empty() || name.find_first_of("\r\n") != std::string::npos ||
			    name.find('\0') != std::string::npos || l < 0 || l >= 360 || x < 0 ||
			    !std::isfinite(std::hypot(x, z)) || std::hypot(x, z) <= 0)
				throw std::runtime_error("Invalid JPL radar station");
			output << code << ' ' << code << " - " << l << ' ' << x << ' ' << z << ' ' << name << '\n';
			++added;
		}
		const auto result = output.str();
		std::istringstream updated(result);
		if (readRadarCatalog(updated).stations.size() != catalog.stations.size() + added)
			throw std::runtime_error("Radar catalog validation failed");
		return {result, used.size() + added, withoutMpc + added, current.empty() || !sameCatalog(current, result)};
	}

	ObservatoryUpdate downloadOpticalObservatoryUpdate(std::string_view current, std::stop_token cancel)
	{
		const auto response = trajecta::http::get("https://www.minorplanetcenter.net/iau/lists/ObsCodes.html", 30,
		                                          4 * 1024 * 1024, cancel);
		return prepareOpticalObservatoryUpdate(response.body, current);
	}

	ObservatoryUpdate downloadRadarObservatoryUpdate(std::string_view current, std::stop_token cancel)
	{
		const auto response = trajecta::http::get("https://ssd-api.jpl.nasa.gov/sb_radar.api?coords=1&c-coords=1", 60,
		                                          8 * 1024 * 1024, cancel);
		return prepareRadarObservatoryUpdate(response.body, current);
	}
}  // namespace orbital_core::obs
