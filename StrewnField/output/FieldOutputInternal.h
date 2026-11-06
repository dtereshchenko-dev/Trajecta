#pragma once

#include "strewn_core/FieldOutput.h"
#include <nlohmann/json.hpp>

namespace strewn::output {

	using Json = nlohmann::json;

	struct ImpactDiagnostic {
		std::uint64_t runId{}, fragmentId{};
		double timeS{}, speedMps{};
	};

	struct ArchivedField {
		std::vector<FieldPoint> points;
		std::vector<ImpactDiagnostic> impactDetails;
		Json summary, job;
		std::string impactFingerprint;
	};

	MapOptions readOptions(const Json& value, const std::filesystem::path& base);
	Json optionsJson(const MapOptions& value);
	Json reportJson(const ArtifactReport& report);
	const char* statusName(ArtifactStatus status);
	ArchivedField readArchive(const std::filesystem::path& directory, std::stop_token cancel);
	Json exportFields(const std::filesystem::path& directory, const DensityFields& fields, std::stop_token cancel);
	Json embedBasemap(const Json& bounds, const MapOptions& options, const std::filesystem::path& cache,
	                  ArtifactReport& report, std::stop_token cancel);
	void writeViewer(const std::filesystem::path& directory, const std::string& title, Json data);

}  // namespace strewn::output
