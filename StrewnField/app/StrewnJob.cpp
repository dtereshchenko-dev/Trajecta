#include "JobFiles.h"
#include "EnvironmentLoader.h"
#include "HttpClient.h"
#include "AutoEnvironment.h"
#include "../output/FieldOutputInternal.h"
#include <optional>
#include <stdexcept>

namespace strewn {
	using namespace jobfiles;

	JobResult StrewnJob::run(const JobSpec& spec, const RunControl& control) const
	{
		const auto directory = std::filesystem::absolute(spec.outputDirectory).lexically_normal();
		if (spec.outputDirectory.empty() ||
		    (std::filesystem::exists(directory) && !std::filesystem::is_empty(directory)))
			throw std::invalid_argument("Output directory must be new or empty");
		if (spec.maps)
			validate(*spec.maps);
		if (control.cancel.stop_requested()) {
			RunSummary summary;
			summary.requested = spec.monteCarlo.runs;
			summary.cancelled = true;
			return {summary, {}, directory};
		}
		std::vector<EntryState> entries;
		if (const auto* state = std::get_if<EntryState>(&spec.input)) {
			entries.push_back(*state);
		} else if (const auto* samples = std::get_if<EntrySamples>(&spec.input)) {
			entries = samples->values;
		} else {
			const auto data = readJson(std::get<EntryFile>(spec.input).path);
			const bool legacy = !data.contains("schema_version");
			if (!legacy && (data.at("schema_version") != 1 || data.at("frame") != "ECEF"))
				throw std::invalid_argument("Unsupported entry sample schema or frame");
			for (const auto& row : data.at("entries"))
				entries.push_back(readEntry(row, entries.size(), legacy));
		}
		Json environmentSnapshot;
		bool idealized = std::holds_alternative<IdealizedEnvironment>(spec.environment);
		auto prepareEnvironment = [&] {
			if (idealized) {
				auto prepared =
				    Environment::idealized(std::get<IdealizedEnvironment>(spec.environment).groundHeightWgs84M);
				environmentSnapshot = environmentJson(prepared);
				environmentSnapshot["source"] = "US76 table / zero wind / constant WGS84 surface";
				return prepared;
			}
			if (const auto* remote = std::get_if<EnvironmentUrls>(&spec.environment)) {
				const auto cache = spec.environmentDirectory.empty() ? directory / "inputs" : spec.environmentDirectory;
				environmentSnapshot = downloadEnvironment(*remote, cache, control.cancel);
			} else if (const auto* automatic = std::get_if<AutomaticEnvironment>(&spec.environment)) {
				const auto cache = spec.environmentDirectory.empty() ? directory / "inputs" : spec.environmentDirectory;
				environmentSnapshot = sources::automaticEnvironment(*automatic, entries, cache, control.cancel);
				idealized = environmentSnapshot.value("idealized_fallback", false);
			} else {
				environmentSnapshot = readJson(std::get<EnvironmentFile>(spec.environment).path);
			}
			return readEnvironment(environmentSnapshot);
		};
		std::optional<Environment> prepared;
		try {
			prepared.emplace(prepareEnvironment());
		} catch (const trajecta::http::Cancelled&) {
			RunSummary summary;
			summary.requested = spec.monteCarlo.runs;
			summary.cancelled = true;
			return {summary, {}, directory};
		}
		auto& environment = *prepared;
		Simulator simulator(environment, spec.physics, spec.integration);
		std::filesystem::create_directories(directory);
		Json inputSnapshot = {{"schema_version", 1}, {"frame", "ECEF"}, {"entries", Json::array()}};
		for (const auto& entry : entries)
			inputSnapshot["entries"].push_back(entryJson(entry));
		writeJson(directory / "entries.json", inputSnapshot);
		writeJson(directory / "environment.json", environmentSnapshot);
		writeJson(directory / "job.json", {{"schema_version", 2},
		                                   {"input", {{"kind", "state_samples"}, {"path", "entries.json"}}},
		                                   {"environment", {{"mode", "file"}, {"path", "environment.json"}}},
		                                   {"physics", physicsJson(spec.physics)},
		                                   {"integration", integrationJson(spec.integration)},
		                                   {"monte_carlo", samplingJson(spec.monteCarlo)},
		                                   {"maps", spec.maps ? output::optionsJson(*spec.maps) : Json(nullptr)},
		                                   {"output_directory", "replay"}});
		Json manifest = {
		    {"schema_version", 1},
		    {"engine_version", "0.2.0"},
		    {"status", "running"},
		    {"rng", "mt19937_64 / nested SplitMix64 v1 / top-53-bit uniform / streams 0,1,2"},
		    {"conditioning", "provided atmospheric entry states; not an Earth impact probability"},
		    {"weights",
		     "entry weights used in categorical sampling; run_weight=1 is unnormalized; multiplicity is separate"},
		    {"storage", "float64 round-trip CSV and JSONL; no restart checkpoints in version 0.1"},
		    {"idealized_environment", idealized},
		    {"shape", shapeName(spec.physics.fragmentShape)},
		    {"shape_warning", spec.physics.fragmentShape == Shape::Sphere
		                          ? ""
		                          : "Experimental fixed area, shared spherical Cd(M), no attitude model"}};
		writeJson(directory / "manifest.json", manifest);
		JobResult result;
		result.directory = directory;
		try {
			ResultArchive archive(directory);
			auto& summary = result.simulation;
			summary =
			    simulator.run(entries, spec.monteCarlo, control, [&](const RunResult& run) { archive.append(run); });
			writeJson(directory / "summary.json", summaryJson(summary));
			manifest["status"] =
			    summary.cancelled ? "cancelled" : (summary.failed || summary.truncated ? "partial" : "complete");
			writeJson(directory / "manifest.json", manifest);
		} catch (const std::exception& error) {
			manifest["status"] = "failed";
			manifest["error"] = error.what();
			writeJson(directory / "manifest.json", manifest);
			throw;
		}
		if (spec.maps) {
			if (result.simulation.cancelled) {
				result.artifacts.status = ArtifactStatus::Cancelled;
				result.artifacts.error = "Maps not generated: calculation cancelled";
			} else if (result.simulation.failed || result.simulation.truncated) {
				result.artifacts.error = "Maps not generated: incomplete ensemble (" +
				                         std::to_string(result.simulation.failed) + " failed, " +
				                         std::to_string(result.simulation.truncated) +
				                         " particle-limited realizations). Numerical results saved.";
				if (result.simulation.truncated)
					result.artifacts.error += " Disable (0) or increase physics.max_particles (currently " +
					                          std::to_string(spec.physics.maxParticles) +
					                          ") and rerun with the same inputs and seed.";
				if (result.simulation.wallTimeouts)
					result.artifacts.error += " " + std::to_string(result.simulation.wallTimeouts) +
					                          " realizations exceeded monte_carlo.run_timeout_s (" +
					                          std::to_string(spec.monteCarlo.runTimeoutS) + " s).";
			} else
				try {
					const auto cache = spec.environmentDirectory.empty() ? directory / "inputs" / "imagery"
					                                                     : spec.environmentDirectory / "imagery";
					result.artifacts = renderSavedRun(directory, directory / "maps", *spec.maps, control.cancel, cache);
				} catch (const std::exception& error) {
					result.artifacts.status =
					    control.cancel.stop_requested() ? ArtifactStatus::Cancelled : ArtifactStatus::Failed;
					result.artifacts.error = error.what();
				}
			manifest["artifacts"] = output::reportJson(result.artifacts);
			writeJson(directory / "manifest.json", manifest);
		}
		return result;
	}

}  // namespace strewn
