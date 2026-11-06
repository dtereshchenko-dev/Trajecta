#include "orbital_core/RuntimeData.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"
#include "AstrometryBias.h"
#include "earth_orientation/EopData.h"
#include "obs/ObservatoryCatalog.h"
#include "obs/Vfcc17Rules.h"

#include <fstream>
#include <iterator>
#include <utility>

namespace orbital_core {

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
		inspect("planetary", base / "linux_p1550p2650.440", [](auto& e) {
			const auto* source = dynamic_cast<const DE440Source*>(getDefaultEphemerisSource());
			e.active = planetaryEphemerisAvailable();
			e.valid = e.active;
			if (e.active)
				e.path = source->sourcePath();
			else
				e.error = "Planetary ephemeris not connected";
		});
		inspect("asteroids", base / "sb441-n16.bsp", [&](auto& e) {
			const auto* source = getDefaultAsteroidSource();
			e.active = asteroidEphemerisAvailable();
			e.valid = e.active;
			if (e.active) {
				e.path = std::filesystem::u8path(source->sourcePath());
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
		inspect("finals", base / "earth_orientation" / "finals2000A.all", [](auto& e) {
			earth_orientation_detail::FinalsData data(e.path);
			e.valid = data.isAvailable();
			e.records = data.size();
			e.firstMjdUtc = data.firstMjdUtc();
			e.lastMjdUtc = data.lastMjdUtc();
		});
		inspect("c04", base / "earth_orientation" / "eopc04.62-now", [&](auto& e) {
			earth_orientation_detail::C04Data data(base / "earth_orientation", EarthOrientationProfile::IAU1980);
			e.valid = data.isAvailable();
			e.records = data.size();
			if (!e.valid && !std::filesystem::is_regular_file(base / "earth_orientation" / "eopc04.cnc"))
				e.error = "C04 consistency file is missing: eopc04.cnc";
		});
		inspect("delta_t", base / "earth_orientation" / "ET-UT.dat", [](auto& e) {
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
		inspect("observatories", base / "OBSCODE.dat", [](auto& e) {
			e.records = obs::observatoryCatalog().size();
			e.valid = e.records > 0;
			e.active = e.valid;
		});
		inspect("radar", base / "RADCODE.dat", [](auto& e) {
			e.records = obs::radarObservatoryCatalog().size();
			e.valid = e.records > 0;
			e.active = e.valid;
		});
		return entries;
	}
}  // namespace orbital_core
