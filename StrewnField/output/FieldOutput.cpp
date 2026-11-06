#include "FieldOutputInternal.h"
#include "ReportData.h"
#include "ReportPreview.h"
#include "../app/JobFiles.h"
#include <chrono>
#include <algorithm>

namespace strewn {

	MapOptions loadMapOptions(const std::filesystem::path& path)
	{
		return output::readOptions(jobfiles::readJson(path), std::filesystem::absolute(path).parent_path());
	}

	ArtifactReport renderSavedRun(const std::filesystem::path& runDirectory,
	                              const std::filesystem::path& outputDirectory, const MapOptions& options,
	                              std::stop_token cancel, const std::filesystem::path& cacheDirectory)
	{
		validate(options);
		const auto input = std::filesystem::weakly_canonical(runDirectory);
		const auto output = std::filesystem::weakly_canonical(outputDirectory);
		if (outputDirectory.empty() || input == output ||
		    (std::filesystem::exists(output) && !std::filesystem::is_empty(output)))
			throw std::invalid_argument("Artifact directory must be new or empty and distinct from the run archive");
		ArtifactReport report;
		const auto start = std::chrono::steady_clock::now();
		if (cancel.stop_requested()) {
			report.status = ArtifactStatus::Cancelled;
			return report;
		}
		std::filesystem::create_directories(output);
		output::Json provenance = {{"schema_version", 1},
		                           {"options", output::optionsJson(options)},
		                           {"run_directory", jobfiles::pathText(input)}};
		try {
			const auto archive = output::readArchive(input, cancel);
			provenance["impacts_fnv1a64"] = archive.impactFingerprint;
			provenance["realizations"] = archive.summary.at("completed");
			provenance["tracked_impacts"] = archive.points.size();
			output::Json data = {{"schema_version", 2},
			                     {"title", options.title},
			                     {"summary", archive.summary},
			                     {"job", archive.job},
			                     {"impacts_fnv1a64", archive.impactFingerprint},
			                     {"runs", output::Json::array()},
			                     {"impacts", output::Json::array()},
			                     {"bounds", nullptr},
			                     {"grid", nullptr},
			                     {"contours", {{"type", "FeatureCollection"}, {"features", output::Json::array()}}}};
			std::map<std::uint64_t, std::vector<std::size_t>> pointRuns;
			for (std::size_t i = 0; i < archive.points.size(); ++i)
				pointRuns[archive.impactDetails[i].runId].push_back(i);
			const auto pointLimit = output::previewLimit(archive.summary.at("completed"));
			std::vector<std::size_t> preview;
			for (auto& [id, indices] : pointRuns) {
				std::sort(indices.begin(), indices.end(), [&](auto a, auto b) {
					return archive.points[a].massKg != archive.points[b].massKg
					           ? archive.points[a].massKg > archive.points[b].massKg
					           : a < b;
				});
				const auto count = std::min(indices.size(), pointLimit);
				for (std::size_t j = 0; j < count; ++j)
					preview.push_back(indices[count == 1 ? 0 : j * (indices.size() - 1) / (count - 1)]);
			}
			std::sort(preview.begin(), preview.end());
			for (const auto i : preview) {
				const auto& p = archive.points[i];
				const auto& d = archive.impactDetails[i];
				data["impacts"].push_back({p.latitudeDeg, p.longitudeDeg, p.massKg, p.multiplicity, p.runWeight,
				                           d.runId, d.fragmentId, d.timeS, d.speedMps});
			}
			data["impact_preview"] = {
			    {"total", archive.points.size()},
			    {"shown", preview.size()},
			    {"method", "mass-rank-stratified points per realization; not used for statistics or KDE"}};
			data["impact_histograms"] = {
			    {"mass", output::impactHistogram(
			                 archive, [&](auto i) { return archive.points[i].massKg; }, true)},
			    {"speed", output::impactHistogram(
			                  archive, [&](auto i) { return archive.impactDetails[i].speedMps; }, false)},
			    {"time",
			     output::impactHistogram(archive, [&](auto i) { return archive.impactDetails[i].timeS; }, false)}};
			data["impact_columns"] = {"latitude_deg", "longitude_deg", "mass_kg", "multiplicity", "run_weight",
			                          "run_id",       "fragment_id",   "time_s",  "speed_mps"};
			if (!archive.points.empty()) {
				const auto fields = estimateField(archive.points, options, cancel);
				auto resolved = options;
				resolved.cellM = fields.cellM;
				resolved.autoCellSize = false;
				provenance["resolved_options"] = output::optionsJson(resolved);
				if (fields.cellM != options.cellM)
					report.warnings.push_back("Map cell size increased from " + std::to_string(options.cellM) + " to " +
					                          std::to_string(fields.cellM) +
					                          " m to fit grid/FFT limits; kernel bandwidths unchanged");
				const auto contours = output::exportFields(output, fields, cancel);
				for (const auto* file : {"density.tif", "density.json", "contours.geojson"})
					report.files.push_back(output / file);
				if (fields.count.fixedKernelFallback || fields.mass.fixedKernelFallback)
					report.warnings.push_back(
					    "Degenerate sample covariance: fixed kernel used for the affected density");
				data["contours"] = contours;
				data["grid"] = {{"encoding", "base64-deflate-byte-shuffle-float32-le"},
				                {"rows", fields.rows},
				                {"columns", fields.columns},
				                {"cell_m", fields.cellM},
				                {"latitude", fields.firstLatitudeDeg},
				                {"longitude", fields.firstLongitudeDeg},
				                {"dlat", fields.latitudeStepDeg},
				                {"dlon", fields.longitudeStepDeg},
				                {"count", output::encodeFloats(fields.count.probability)},
				                {"mass", output::encodeFloats(fields.mass.probability)},
				                {"metadata", jobfiles::readJson(output / "density.json")}};
				double west = 180, east = -180, south = 90, north = -90;
				const double countCut =
				    *std::max_element(fields.count.probability.begin(), fields.count.probability.end()) * 2e-4;
				const double massCut =
				    *std::max_element(fields.mass.probability.begin(), fields.mass.probability.end()) * 2e-4;
				for (std::size_t y = 0; y < fields.rows; ++y)
					for (std::size_t x = 0; x < fields.columns; ++x) {
						const auto i = y * fields.columns + x;
						if (fields.count.probability[i] < countCut && fields.mass.probability[i] < massCut)
							continue;
						west = std::min(west, fields.firstLongitudeDeg + (x - .5) * fields.longitudeStepDeg);
						east = std::max(east, fields.firstLongitudeDeg + (x + .5) * fields.longitudeStepDeg);
						south = std::min(south, fields.firstLatitudeDeg + (y - .5) * fields.latitudeStepDeg);
						north = std::max(north, fields.firstLatitudeDeg + (y + .5) * fields.latitudeStepDeg);
					}
				for (const auto& p : archive.points) {
					west = std::min(west, p.longitudeDeg);
					east = std::max(east, p.longitudeDeg);
					south = std::min(south, p.latitudeDeg);
					north = std::max(north, p.latitudeDeg);
				}
				for (const auto& feature : contours["features"])
					for (const auto& point : feature["geometry"]["coordinates"]) {
						west = std::min(west, point[0].get<double>());
						east = std::max(east, point[0].get<double>());
						south = std::min(south, point[1].get<double>());
						north = std::max(north, point[1].get<double>());
					}
				const double dx = std::max((east - west) * .05, fields.longitudeStepDeg * 4),
				             dy = std::max((north - south) * .05, fields.latitudeStepDeg * 4);
				data["bounds"] = {{south - dy, west - dx}, {north + dy, east + dx}};
			}
			try {
				data.update(output::buildReportData(input, archive, report, cancel));
			} catch (const std::exception& error) {
				if (cancel.stop_requested())
					throw;
				report.warnings.push_back(std::string("Ensemble diagnostics unavailable: ") + error.what());
			}
			data["tiles"] = output::embedBasemap(
			    data["bounds"], options, cacheDirectory.empty() ? output / "imagery" : cacheDirectory, report, cancel);
			data["warnings"] = report.warnings;
			data["options"] = output::optionsJson(options);
			if (cancel.stop_requested())
				throw std::runtime_error("Report cancelled");
			output::writeViewer(output, options.title, std::move(data));
			report.files.push_back(output / "field.html");
			report.status = archive.points.empty() ? ArtifactStatus::Empty : ArtifactStatus::Complete;
		} catch (const std::exception& error) {
			report.status = cancel.stop_requested() ? ArtifactStatus::Cancelled : ArtifactStatus::Failed;
			report.error = error.what();
		}
		report.elapsedS = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
		report.files.push_back(output / "artifacts.json");
		provenance["report"] = output::reportJson(report);
		jobfiles::writeJson(output / "artifacts.json", provenance);
		return report;
	}

}  // namespace strewn
