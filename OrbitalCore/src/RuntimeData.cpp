#include "orbital_core/RuntimeData.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"
#include "AstrometryBias.h"
#include "earth_orientation/EopData.h"
#include "earth_orientation/EarthPck.h"
#include "obs/ObservatoryCatalog.h"
#include "obs/Vfcc17Rules.h"
#include "EnvironmentVariable.h"

#include <fstream>
#include <iterator>
#include <utility>

namespace orbital_core {
	namespace {
		std::filesystem::path defaultObservatoryDirectory()
		{
#ifdef _WIN32
			if (const auto local = environmentVariable(L"LOCALAPPDATA"); local && !local->empty())
				return std::filesystem::path(*local) / "Trajecta" / "Trajecta" / "observatories";
#elif defined(__APPLE__)
			if (const auto home = environmentVariable("HOME"); home && !home->empty())
				return std::filesystem::path(*home) / "Library" / "Application Support" / "Trajecta" / "Trajecta" /
				       "observatories";
#else
			if (const auto data = environmentVariable("XDG_DATA_HOME"); data && !data->empty())
				return std::filesystem::path(*data) / "Trajecta" / "Trajecta" / "observatories";
			if (const auto home = environmentVariable("HOME"); home && !home->empty())
				return std::filesystem::path(*home) / ".local" / "share" / "Trajecta" / "Trajecta" / "observatories";
#endif
			return {};
		}

		std::filesystem::path observatoryDirectory = defaultObservatoryDirectory();
		std::filesystem::path ephemerisDirectory =
		    observatoryDirectory.empty() ? std::filesystem::path{} : observatoryDirectory.parent_path() / "ephemerides";
		std::filesystem::path pckDirectory = observatoryDirectory.empty()
		                                         ? std::filesystem::path{}
		                                         : observatoryDirectory.parent_path() / "earth_orientation";
	}  // namespace

	void setObservatoryDataDirectory(std::filesystem::path directory) { observatoryDirectory = std::move(directory); }

	std::filesystem::path observatoryDataDirectory() { return observatoryDirectory; }

	void setEarthPckDataDirectory(std::filesystem::path directory) { pckDirectory = std::move(directory); }

	std::filesystem::path earthPckDataDirectory() { return pckDirectory; }

	std::filesystem::path earthC04DataDirectory()
	{
		if (!pckDirectory.empty())
			return pckDirectory;
		return getExecutableDir() / "earth_orientation";
	}

	std::filesystem::path earthPckFilePath(const std::string& filename)
	{
		const auto directory = earthPckDataDirectory();
		const auto updated = directory / filename;
		if (!directory.empty())
			return updated;
		return getExecutableDir() / "earth_orientation" / filename;
	}

	void setEphemerisDataDirectory(std::filesystem::path directory) { ephemerisDirectory = std::move(directory); }

	std::filesystem::path ephemerisDataDirectory() { return ephemerisDirectory; }

	bool planetaryEphemerisAvailable() noexcept
	{
		const auto* source = dynamic_cast<const DE440Source*>(getDefaultEphemerisSource());
		return source && source->isAvailable();
	}

	bool asteroidEphemerisAvailable() noexcept
	{
		const auto* source = getDefaultAsteroidSource();
		return source && source->isAvailable();
	}

	std::vector<RuntimeDataEntry> inspectRuntimeData(const std::filesystem::path& downloadedAsteroidFile)
	{
		const auto base = getExecutableDir();
		std::vector<RuntimeDataEntry> entries;
		const auto inspect = [&](std::string id, std::filesystem::path path, auto&& fn) {
			RuntimeDataEntry entry;
			entry.id = std::move(id);
			entry.path = std::move(path);
			try {
				fn(entry);
				if (!entry.valid && entry.error.empty())
					entry.error = "Invalid or unavailable data";
			} catch (const std::exception& error) {
				entry.error = error.what();
			}
			entries.push_back(std::move(entry));
		};
		inspect("planetary", ephemerisDataDirectory() / "linux_p1550p2650.440", [](auto& e) {
			const auto* source = dynamic_cast<const DE440Source*>(getDefaultEphemerisSource());
			e.active = planetaryEphemerisAvailable();
			e.valid = e.active;
			if (e.active) {
				const auto path = source->sourcePath();
				e.path = std::filesystem::path(std::u8string(path.begin(), path.end()));
			} else
				e.error = "Planetary ephemeris not connected";
		});
		inspect("asteroids", base / "sb441-n16.bsp", [&](auto& e) {
			const auto* source = getDefaultAsteroidSource();
			e.active = asteroidEphemerisAvailable();
			e.valid = e.active;
			if (e.active) {
				const auto path = source->sourcePath();
				e.path = std::filesystem::path(std::u8string(path.begin(), path.end()));
				e.records = source->targetIds().size();
			} else if (!downloadedAsteroidFile.empty() && std::filesystem::is_regular_file(downloadedAsteroidFile)) {
				e.path = downloadedAsteroidFile;
				const auto utf8 = e.path.u8string();
				const SpkReader downloaded(std::string(utf8.begin(), utf8.end()));
				e.records = downloaded.targetIds().size();
				e.valid = downloaded.isAvailable() && e.records == 16;
				if (!e.valid)
					e.error = "Downloaded asteroid ephemeris is invalid";
			} else
				e.error = "Asteroid ephemeris not connected";
		});
		for (const auto& [id, file] :
		     {std::pair{"earth_pck_historic", "earth_historic.bpc"}, std::pair{"earth_pck_latest", "earth_latest.bpc"},
		      std::pair{"earth_pck_predict", "earth_predict.bpc"}}) {
			inspect(id, earthPckFilePath(file), [file](auto& e) {
				e.active = std::filesystem::is_regular_file(e.path) && earthPckLoaded(file);
				e.valid = e.active;
				if (!e.valid)
					e.error = "Earth PCK not loaded";
			});
		}
		inspect("c04", earthC04DataDirectory() / "eopc04.62-now", [&](auto& e) {
			earth_orientation_detail::C04Data data(earthC04DataDirectory(), EarthOrientationProfile::IAU1980);
			e.valid = data.isAvailable();
			e.records = data.size();
			if (!e.valid && !std::filesystem::is_regular_file(earthC04DataDirectory() / "eopc04.cnc"))
				e.error = "C04 consistency file is missing: eopc04.cnc";
		});
		inspect("delta_t", earthC04DataDirectory() / "ET-UT.dat", [](auto& e) {
			earth_orientation_detail::DeltaTData data(e.path);
			e.records = data.size();
			e.valid = e.records > 0;
		});
		inspect("weights", base / "weights" / "vfcc17.rules", [](auto& e) {
			std::ifstream input(e.path, std::ios::binary);
			if (!input) {
				e.error = "Cannot read file";
				return;
			}
			const std::string text{std::istreambuf_iterator<char>(input), {}};
			obs::Vfcc17Rules rules;
			e.valid = rules.load(text, e.error);
			e.records = rules.size();
		});
		for (const auto& [id, scheme] : {std::pair{"debias_fcct14", Calculations::AstrometryBiasScheme::Fcct14},
		                                 std::pair{"debias_efcc18", Calculations::AstrometryBiasScheme::Efcc18}}) {
			inspect(id, Calculations::astrometryBiasPathForScheme(base, scheme), [](auto& e) {
				AstrometryBias data;
				data.load(e.path);
				e.valid = data.loaded();
				e.records = data.rowCount();
				e.catalogs = data.catalogCount();
			});
		}
		inspect("observatories", observatoryDataDirectory() / "OBSCODE.dat", [](auto& e) {
			e.records = obs::observatoryCatalog().size();
			e.path = obs::activeOpticalCatalogPath();
			e.valid = e.records > 0;
			e.active = e.valid;
		});
		inspect("radar", observatoryDataDirectory() / "RADCODE.dat", [](auto& e) {
			e.records = obs::radarObservatoryCatalog().size();
			e.path = obs::activeRadarCatalogPath();
			e.valid = e.records > 0;
			e.active = e.valid;
		});
		return entries;
	}
}  // namespace orbital_core
