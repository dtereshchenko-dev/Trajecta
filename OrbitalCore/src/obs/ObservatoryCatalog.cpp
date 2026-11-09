#include "obs/ObservatoryCatalog.h"

#include "FileReader.h"
#include "logger/Logger.h"
#include "orbital_core/Types.h"
#include "orbital_core/RuntimeData.h"
#include "orbital_core/obs/ObservatoryUpdates.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace orbital_core::obs {
	namespace {

		std::string readCatalogFile(const std::filesystem::path& path)
		{
			FileReader reader;
			reader.read(path);
			if (reader.error() != FileReader::Error::NoError)
				throw std::runtime_error(reader.errorString() + ": " + path.string());
			return reader.data();
		}

		ObservatoryCatalog loadOpticalCatalog(const std::filesystem::path& path)
		{
			std::istringstream stream(readCatalogFile(path));
			ObservatoryCatalog catalog;
			std::string line;
			std::getline(stream, line);
			while (std::getline(stream, line)) {
				if (line.empty() || line.front() == '!')
					continue;
				std::istringstream row(line);
				row.imbue(std::locale::classic());
				ObservatoryRecord record;
				if (!(row >> record.code))
					continue;
				row >> record.longitude >> record.cosPhi >> record.sinPhi;
				row.clear();
				std::getline(row, record.name);
				record.name.erase(record.name.begin(), std::ranges::find_if(record.name, [](unsigned char ch) {
					                  return !std::isspace(ch);
				                  }));
				const auto code = record.code;
				catalog.emplace(code, std::move(record));
			}
			if (catalog.empty())
				throw std::runtime_error(path.filename().string() + " is empty or malformed");
			LOG_INFO << "ObservatoryCatalog: loaded " << catalog.size() << " stations from " << path.string();
			return catalog;
		}

		template <typename Catalog>
		struct LoadedCatalog {
			Catalog data;
			std::filesystem::path path;
		};

		const LoadedCatalog<ObservatoryCatalog>& opticalCatalogData()
		{
			static const LoadedCatalog<ObservatoryCatalog> data = [] {
				const auto directory = observatoryDataDirectory();
				if (directory.empty())
					throw std::runtime_error("Observatory data directory is not configured");
				const auto path = directory / "OBSCODE.dat";
				validateOpticalObservatoryCatalog(readCatalogFile(path));
				return LoadedCatalog<ObservatoryCatalog>{loadOpticalCatalog(path), path};
			}();
			return data;
		}

		const LoadedCatalog<RadarCatalog>& radarCatalogData()
		{
			static const LoadedCatalog<RadarCatalog> data = [] {
				const auto directory = observatoryDataDirectory();
				if (directory.empty())
					throw std::runtime_error("Observatory data directory is not configured");
				const auto path = directory / "RADCODE.dat";
				std::istringstream stream(readCatalogFile(path));
				auto catalog = readRadarCatalog(stream);
				LOG_INFO << "ObservatoryCatalog: loaded " << catalog.stations.size() << " stations from "
				         << path.string();
				return LoadedCatalog<RadarCatalog>{std::move(catalog), path};
			}();
			return data;
		}

	}  // namespace

	RadarCatalog readRadarCatalog(std::istream& stream)
	{
		RadarCatalog data;
		bool versionSeen = false;
		std::size_t lineNumber = 0;
		const auto fail = [&](const std::string& reason) {
			throw std::runtime_error("RADCODE.dat line " + std::to_string(lineNumber) + ": " + reason);
		};
		const auto stationCode = [](std::string_view code) {
			const auto digit = [](char ch) { return ch >= '0' && ch <= '9'; };
			return code.size() == 3 && (digit(code[0]) || (code[0] >= 'A' && code[0] <= 'Z')) && digit(code[1]) &&
			       digit(code[2]);
		};
		const auto jplCode = [](std::string_view code) {
			return code.size() >= 2 && code.front() == '-' && code[1] >= '1' && code[1] <= '9' &&
			       std::ranges::all_of(code.substr(2), [](char ch) { return ch >= '0' && ch <= '9'; });
		};
		const auto coordinate = [&](const std::string& token) {
			double value = 0;
			const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), value);
			if (error != std::errc{} || end != token.data() + token.size() || !std::isfinite(value))
				fail("invalid coordinate " + token);
			return value;
		};
		for (std::string line; std::getline(stream, line);) {
			++lineNumber;
			const auto start = line.find_first_not_of(" \t\r");
			if (start == std::string::npos)
				continue;
			line = line.substr(start, line.find_last_not_of(" \t\r") - start + 1);
			if (line.starts_with("! RADCODE")) {
				if (versionSeen || line != "! RADCODE format=2")
					fail("unsupported or duplicate format header");
				versionSeen = true;
				continue;
			}
			if (line.starts_with("! ALIAS"))
				fail("aliases must be in the station's JPL/MPC columns");
			if (line.front() == '!')
				continue;
			if (!versionSeen)
				fail("expected ! RADCODE format=2 before station rows");
			std::istringstream row(line);
			row.imbue(std::locale::classic());
			ObservatoryRecord record;
			std::string jpl, mpc, longitude, rhoXY, rhoZ;
			if (!(row >> record.code >> jpl >> mpc >> longitude >> rhoXY >> rhoZ))
				fail("malformed station row");
			if ((!stationCode(record.code) && !jplCode(record.code)) || (jpl != "-" && !jplCode(jpl)) ||
			    (mpc != "-" && !stationCode(mpc)))
				fail("invalid station identifier");
			record.longitude = coordinate(longitude);
			record.cosPhi = coordinate(rhoXY);
			record.sinPhi = coordinate(rhoZ);
			const double radius = std::hypot(record.cosPhi, record.sinPhi);
			if (record.longitude < 0 || record.longitude >= 360 || record.cosPhi < 0 || !std::isfinite(radius) ||
			    radius <= 0)
				fail("invalid station coordinates");
			std::getline(row >> std::ws, record.name);
			if (record.name.empty())
				fail("missing station name");
			record.cosPhi *= 1e-10 / Constants::EARTH_RADIUS_IN_AU;
			record.sinPhi *= 1e-10 / Constants::EARTH_RADIUS_IN_AU;
			const auto code = record.code;
			if (!data.stations.emplace(code, std::move(record)).second)
				fail("duplicate station " + code);
			for (const auto& alias : {jpl, mpc}) {
				if (alias == "-" || alias == code)
					continue;
				if (!data.aliases.emplace(alias, code).second)
					fail("duplicate identifier " + alias);
			}
		}
		if (!versionSeen || data.stations.empty())
			fail("missing format header or empty catalog");
		if (stream.bad())
			fail("catalog read failed");
		for (const auto& [alias, target] : data.aliases)
			if (data.stations.contains(alias))
				fail("ambiguous identifier " + alias + " for " + target);
		return data;
	}

	const ObservatoryCatalog& observatoryCatalog() { return opticalCatalogData().data; }

	std::filesystem::path activeOpticalCatalogPath() { return opticalCatalogData().path; }

	const ObservatoryCatalog& radarObservatoryCatalog() { return radarCatalogData().data.stations; }

	std::filesystem::path activeRadarCatalogPath() { return radarCatalogData().path; }

	const ObservatoryRecord& radarObservatory(std::string_view code)
	{
		const auto& data = radarCatalogData().data;
		std::string canonical(code);
		if (const auto alias = data.aliases.find(canonical); alias != data.aliases.end())
			canonical = alias->second;
		const auto station = data.stations.find(canonical);
		if (station == data.stations.end())
			throw std::invalid_argument("Radar station missing from RADCODE.dat: " + std::string(code));
		return station->second;
	}

}  // namespace orbital_core::obs
