#pragma once

#include "strewn_core/Simulator.h"
#include "strewn_core/FieldOutput.h"
#include <filesystem>
#include <optional>

namespace strewn {

	struct EntryFile {
		std::filesystem::path path;
	};

	struct EntrySamples {
		std::vector<EntryState> values;
	};

	struct IdealizedEnvironment {
		double groundHeightWgs84M{};
	};

	struct EnvironmentFile {
		std::filesystem::path path;
	};

	struct EnvironmentUrls {
		std::string weather;
		std::string terrain;
		unsigned timeoutSeconds{60};
		std::size_t maxDownloadBytes{64 * 1024 * 1024};
	};

	struct AutomaticEnvironment {
		double marginKm{25.0};
		double terrainStepDeg{0.02};
		// South, north, west, east; absent derives coverage from entry trajectories.
		std::optional<std::array<double, 4>> boundsDeg;
		std::optional<double> cycleMjdUtc;
		std::filesystem::path gfsFile;
		unsigned timeoutSeconds{120};
		std::size_t maxDownloadBytes{256 * 1024 * 1024};
		bool offline{};
		bool idealizedFallback{};
	};

	struct JobSpec {
		std::variant<EntryState, EntryFile, EntrySamples> input;
		std::variant<IdealizedEnvironment, EnvironmentFile, EnvironmentUrls, AutomaticEnvironment> environment;
		PhysicsConfig physics;
		IntegrationOptions integration;
		MonteCarloConfig monteCarlo;
		std::optional<MapOptions> maps;
		std::filesystem::path outputDirectory;
		// Shared per-asteroid weather/terrain cache; empty uses outputDirectory/inputs.
		std::filesystem::path environmentDirectory;
	};

	struct JobResult {
		RunSummary simulation;
		ArtifactReport artifacts;
		std::filesystem::path directory;
	};

	STREWN_CORE_EXPORT JobSpec loadJob(const std::filesystem::path& path);
	// The producer validates its input section; atmospheric options retain the
	// ordinary schema. Prepared samples are owned by the returned job.
	STREWN_CORE_EXPORT JobSpec loadJob(const std::filesystem::path& path, EntrySamples prepared);

	class STREWN_CORE_EXPORT StrewnJob {
	   public:
		JobResult run(const JobSpec& spec, const RunControl& control = {}) const;
	};

}  // namespace strewn
