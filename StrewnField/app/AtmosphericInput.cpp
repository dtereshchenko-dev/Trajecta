#include "strewn_core/AtmosphericInput.h"
#include "JobFiles.h"
#include "AutoEnvironment.h"
#include "Convergence.h"
#include "VectorMath.h"
#include "../output/FieldOutputInternal.h"
#include <set>
#include <chrono>
#include <atomic>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace strewn {
	using namespace jobfiles;

	JobSpec prepareAtmosphericJob(const std::filesystem::path& mcArchive)
	{
		const auto directory = std::filesystem::absolute(mcArchive).lexically_normal();
		const auto manifest = readJson(directory / "manifest.json");
		const auto request = readJson(directory / "request.json");
		const auto data = readJson(directory / "entries.json");
		if (manifest.at("schema_version") != 1 || manifest.at("status") != "complete" ||
		    !manifest.at("entry_states_ready").get<bool>() || data.at("schema_version") != 1 ||
		    data.at("frame") != "ECEF" || data.at("altitude_km") != 100)
			throw std::invalid_argument("Atmospheric MC requires complete Earth entries at 100 km");
		std::vector<EntryState> entries;
		std::set<std::uint64_t> ids;
		for (const auto& row : data.at("entries")) {
			const auto entry = readEntry(row, entries.size());
			if (!detail::finite(entry.positionEcefM))
				throw std::invalid_argument("Non-finite entry position");
			const auto position = geodetic(entry.positionEcefM);
			const auto up = enuToEcef({0, 0, 1}, position.latitudeRad, position.longitudeRad);
			const double vertical =
			    up[0] * entry.velocityEcefMps[0] + up[1] * entry.velocityEcefMps[1] + up[2] * entry.velocityEcefMps[2];
			if (!ids.insert(entry.id).second || !detail::finite(entry.velocityEcefMps) ||
			    !std::isfinite(entry.epochMjdUtc) || entry.weight != 1.0 || !std::isfinite(position.heightWgs84M) ||
			    std::abs(position.heightWgs84M - 100000) > 1 || !std::isfinite(vertical) || vertical >= 0)
				throw std::invalid_argument("Invalid, duplicate or non-inbound entry state");
			entries.push_back(entry);
		}
		if (entries.empty() || entries.size() != request.at("monte_carlo").at("clones").get<std::size_t>())
			throw std::invalid_argument("Incomplete MC entry archive");
		const auto& physical = request.at("physical");
		if (!physical.at("absolute_magnitude").is_number())
			throw std::invalid_argument("Fit has no absolute magnitude for photometric size");
		JobSpec job;
		job.input = EntryFile{directory / "entries.json"};
		job.environment = AutomaticEnvironment{};
		job.environmentDirectory = directory.parent_path().parent_path() / "StrewnField";
		job.outputDirectory = job.environmentDirectory / "runs" / (directory.filename().string() + "-atmosphere");
		job.physics.representativesPerBin = 3;
		job.monteCarlo.runs = 1000;
		job.monteCarlo.seed = request.at("monte_carlo").at("seed").get<std::uint64_t>();
		job.monteCarlo.parameters.size = PhotometricSize{physical.at("absolute_magnitude").get<double>()};
		applyMaterialPreset(job, MaterialPreset::OrdinaryChondrite);
		job.monteCarlo.convergence = atmosphericConvergenceGrid(entries);
		job.maps = MapOptions{};
		job.maps->title = request.at("object").get<std::string>();
		job.maps->basemap = Basemap::EsriImagery;
		validate(job.monteCarlo);
		return job;
	}

	ConvergenceConfig atmosphericConvergenceGrid(std::span<const EntryState> entries, ConvergenceConfig convergence)
	{
		const auto bounds = sources::environmentBounds(AutomaticEnvironment{}, entries);
		convergence.centerDeg = {(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2};
		convergence.halfSizeKm = {
		    (bounds[1] - bounds[0]) * 111.195 / 2,
		    (bounds[3] - bounds[2]) * 111.195 * std::cos(convergence.centerDeg[0] * detail::deg) / 2};
		while ((2000 * convergence.halfSizeKm[0] / convergence.cellM) *
		           (2000 * convergence.halfSizeKm[1] / convergence.cellM) >
		       2000000)
			convergence.cellM *= 2;
		detail::validateConvergence(convergence);
		return convergence;
	}

	void saveJobConfiguration(const std::filesystem::path& path, const JobSpec& job)
	{
		validate(job.physics, job.integration);
		validate(job.monteCarlo);
		if (job.maps)
			validate(*job.maps);
		const auto base = std::filesystem::absolute(path).parent_path();
		const auto relative = [&](const std::filesystem::path& value) {
			const auto rel = std::filesystem::absolute(value).lexically_relative(base);
			return pathText(rel.empty() ? std::filesystem::absolute(value) : rel);
		};
		Json input;
		if (const auto* file = std::get_if<EntryFile>(&job.input))
			input = {{"kind", "state_samples"}, {"path", relative(file->path)}};
		else if (const auto* entry = std::get_if<EntryState>(&job.input))
			input = {{"kind", "fixed_state"}, {"state", entryJson(*entry)}};
		else
			throw std::invalid_argument("Archive entry samples before saving their configuration");
		Json environment;
		if (const auto* automatic = std::get_if<AutomaticEnvironment>(&job.environment)) {
			sources::validateAutomatic(*automatic);
			environment = {{"mode", "automatic"},
			               {"margin_km", automatic->marginKm},
			               {"terrain_step_deg", automatic->terrainStepDeg},
			               {"timeout_seconds", automatic->timeoutSeconds},
			               {"max_download_bytes", automatic->maxDownloadBytes},
			               {"offline", automatic->offline},
			               {"on_unavailable", automatic->idealizedFallback ? "idealized" : "error"}};
			if (automatic->boundsDeg)
				environment["bounds_deg"] = *automatic->boundsDeg;
			if (automatic->cycleMjdUtc)
				environment["cycle_mjd_utc"] = *automatic->cycleMjdUtc;
			if (!automatic->gfsFile.empty())
				environment["gfs_file"] = relative(automatic->gfsFile);
		} else if (const auto* idealized = std::get_if<IdealizedEnvironment>(&job.environment)) {
			environment = {{"mode", "idealized"}, {"ground_height_wgs84_m", idealized->groundHeightWgs84M}};
		} else if (const auto* file = std::get_if<EnvironmentFile>(&job.environment)) {
			environment = {{"mode", "file"}, {"path", relative(file->path)}};
		} else {
			const auto& remote = std::get<EnvironmentUrls>(job.environment);
			environment = {{"mode", "remote"},
			               {"weather_url", remote.weather},
			               {"terrain_url", remote.terrain},
			               {"timeout_seconds", remote.timeoutSeconds},
			               {"max_download_bytes", remote.maxDownloadBytes}};
		}
		Json value = {{"schema_version", 2},
		              {"input", input},
		              {"environment", environment},
		              {"physics", physicsJson(job.physics)},
		              {"integration", integrationJson(job.integration)},
		              {"monte_carlo", samplingJson(job.monteCarlo)},
		              {"maps", job.maps ? output::optionsJson(*job.maps) : Json(nullptr)},
		              {"output_directory", relative(job.outputDirectory)}};
		if (!job.environmentDirectory.empty())
			value["environment_directory"] = relative(job.environmentDirectory);
		static std::atomic<unsigned> sequence{0};
		auto temporary = path;
		temporary += ".tmp-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "-" +
		             std::to_string(sequence++);
		try {
			{
				std::ofstream stream(temporary, std::ios::binary);
				stream.exceptions(std::ios::failbit | std::ios::badbit);
				stream << value.dump(2) << '\n';
				stream.close();
			}
#ifdef _WIN32
			if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
				throw std::system_error(static_cast<int>(GetLastError()), std::system_category(),
				                        "Save atmospheric config");
#else
			std::filesystem::rename(temporary, path);
#endif
		} catch (...) {
			std::error_code ignored;
			std::filesystem::remove(temporary, ignored);
			throw;
		}
	}

}  // namespace strewn
