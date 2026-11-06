#include "ReportData.h"
#include "EntrySummary.h"
#include "ReportPreview.h"
#include "../app/JobFiles.h"
#include "strewn_core/Simulator.h"
#include "VectorMath.h"
#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <thread>
#include <cpl_conv.h>

namespace strewn::output {
	namespace {
		using Sample = std::array<double, 8>;

		std::vector<Json> readLines(const std::filesystem::path& path, std::stop_token cancel)
		{
			std::ifstream stream(path);
			if (!stream)
				throw std::runtime_error("Cannot read report archive: " + jobfiles::pathText(path));
			std::vector<Json> result;
			std::string line;
			while (std::getline(stream, line)) {
				if (cancel.stop_requested())
					throw std::runtime_error("Report cancelled");
				if (!line.empty())
					result.push_back(Json::parse(line));
			}
			if (stream.bad())
				throw std::runtime_error("Cannot finish reading report archive");
			return result;
		}

		Sample sample(const FlightState& state, double time, const GeodeticPosition& origin,
		              const Environment& environment, double windScale)
		{
			const auto pos = geodetic(state.positionEcefM);
			const auto wind = enuToEcef(environment.air(pos).windEnuMps, pos.latitudeRad, pos.longitudeRad);
			Vector3 relative;
			for (std::size_t i = 0; i < 3; ++i)
				relative[i] = state.velocityEcefMps[i] - windScale * wind[i];
			const double angle = 2 * std::asin(std::sqrt(std::clamp(
			                             std::pow(std::sin((pos.latitudeRad - origin.latitudeRad) / 2), 2) +
			                                 std::cos(origin.latitudeRad) * std::cos(pos.latitudeRad) *
			                                     std::pow(std::sin((pos.longitudeRad - origin.longitudeRad) / 2), 2),
			                             0., 1.)));
			return {time,
			        pos.heightWgs84M / 1000,
			        state.massKg,
			        detail::norm(state.velocityEcefMps) / 1000,
			        6367 * angle,
			        pos.latitudeRad / detail::deg,
			        pos.longitudeRad / detail::deg,
			        detail::norm(relative) / 1000};
		}

		struct Track {
			Json metadata;
			std::vector<Sample> samples;
			Sample last{};
			std::size_t received{}, stride{1};
		};

		std::vector<Sample> compact(std::vector<Sample> rows)
		{
			if (rows.size() <= 96)
				return rows;
			auto normalized = rows;
			for (auto& row : normalized) {
				row[0] = std::log1p(row[0]);
				row[2] = std::log(std::max(1e-30, row[2]));
			}
			for (std::size_t i = 1; i < normalized.size(); ++i)
				normalized[i][6] = normalized[i - 1][6] + std::remainder(rows[i][6] - rows[i - 1][6], 360.);
			for (std::size_t c = 0; c < Sample{}.size(); ++c) {
				double lo = normalized.front()[c], hi = lo;
				for (const auto& row : normalized) {
					lo = std::min(lo, row[c]);
					hi = std::max(hi, row[c]);
				}
				for (auto& row : normalized)
					row[c] = (row[c] - lo) / std::max(1e-12, hi - lo);
			}
			std::vector<std::size_t> indices{0, rows.size() - 1};
			while (indices.size() < 96) {
				double largest = -1;
				std::size_t chosen = 0, insert = 0;
				for (std::size_t s = 1; s < indices.size(); ++s) {
					const auto a = indices[s - 1], b = indices[s];
					for (auto i = a + 1; i < b; ++i) {
						const double t = (normalized[i][0] - normalized[a][0]) /
						                 std::max(1e-15, normalized[b][0] - normalized[a][0]);
						double error = 0;
						for (std::size_t c = 1; c < Sample{}.size(); ++c)
							error = std::max(
							    error, std::abs(normalized[i][c] - std::lerp(normalized[a][c], normalized[b][c], t)));
						if (error > largest) {
							largest = error;
							chosen = i;
							insert = s;
						}
					}
				}
				indices.insert(indices.begin() + insert, chosen);
			}
			std::vector<Sample> result;
			for (auto i : indices)
				result.push_back(rows[i]);
			return result;
		}

		Json replayRun(const Json& run, const Json& breakups, const EntryState& entry, const JobSpec& spec,
		               const Environment& environment, const ArchivedField& archive,
		               const std::map<std::uint64_t, std::size_t>& expected, std::stop_token cancel,
		               std::size_t traceLimit)
		{
			RunParameters p;
			const auto& j = run.at("parameters");
			p.massKg = j.at("mass_kg");
			p.densityKgM3 = j.at("density_kg_m3");
			p.strengthPa = j.at("strength_pa");
			p.dustFraction = j.at("dust_fraction");
			p.largestFraction = j.at("largest_fraction");
			p.beta = j.at("beta");
			p.windScale = j.at("wind_scale");
			p.dragScale = j.at("drag_scale");
			p.ablationKgJ = j.at("ablation_kg_j");
			const auto origin = geodetic(entry.positionEcefM);
			std::map<std::uint64_t, Track> tracks;
			const auto replay = Simulator(environment, spec.physics, spec.integration)
			                        .runOne(
			                            entry, p, spec.monteCarlo.seed, run.at("run_id"), cancel,
			                            [&](const Fragment& fragment, double time, const FlightState& state) {
				                            auto& track = tracks[fragment.id];
				                            if (track.metadata.is_null())
					                            track.metadata = {{"id", fragment.id},
					                                              {"parent_id", fragment.parentId},
					                                              {"generation", fragment.generation},
					                                              {"multiplicity", fragment.multiplicity},
					                                              {"strength_pa", fragment.strengthPa}};
				                            if (track.received && time == track.last[0])
					                            return;
				                            track.last = sample(state, time, origin, environment, p.windScale);
				                            if (track.received++ % track.stride != 0)
					                            return;
				                            track.samples.push_back(track.last);
				                            if (track.samples.size() == 512) {
					                            for (std::size_t i = 1; i < 256; ++i)
						                            track.samples[i] = track.samples[i * 2];
					                            track.samples.resize(256);
					                            track.stride *= 2;
				                            }
			                            },
			                            spec.monteCarlo.runTimeoutS);
			if (replay.status != RunStatus::Complete || expected.size() != replay.impacts.size())
				throw std::runtime_error("Replay does not reproduce archived impact count/status");
			double maxDelta = 0;
			for (const auto& impact : replay.impacts) {
				if (!expected.contains(impact.fragment.id))
					throw std::runtime_error("Replay fragment mismatch");
				const auto i = expected.at(impact.fragment.id);
				const auto& original = archive.points[i];
				const double delta =
				    111195 * std::hypot(impact.latitudeDeg - original.latitudeDeg,
				                        std::remainder(impact.longitudeDeg - original.longitudeDeg, 360.) *
				                            std::cos(origin.latitudeRad));
				maxDelta = std::max(maxDelta, delta);
				if (delta > .01 ||
				    std::abs(impact.fragment.state.massKg - original.massKg) > 1e-10 * std::max(1., original.massKg) ||
				    std::abs(impact.fragment.timeS - archive.impactDetails[i].timeS) > 1e-6 ||
				    impact.fragment.multiplicity != original.multiplicity)
					throw std::runtime_error("Replay differs from archived impacts");
				tracks.at(impact.fragment.id).metadata["end_event"] = "ground";
			}
			const auto& mass = run.at("mass_balance_kg");
			const std::array<std::pair<const char*, double>, 7> ledger{{{"initial", replay.mass.initialKg},
			                                                            {"landed", replay.mass.landedKg},
			                                                            {"ablated", replay.mass.ablatedKg},
			                                                            {"dust", replay.mass.dustKg},
			                                                            {"residual", replay.mass.residualKg},
			                                                            {"below_cutoff", replay.mass.belowCutoffKg},
			                                                            {"unresolved", replay.mass.unresolvedKg}}};
			for (const auto& [key, value] : ledger)
				if (std::abs(value - mass.at(key).get<double>()) > 1e-9 * std::max(1., replay.mass.initialKg))
					throw std::runtime_error("Replay mass ledger mismatch");
			if (replay.breakups.size() != breakups.size())
				throw std::runtime_error("Replay breakup count mismatch");
			Json events = Json::array();
			for (std::size_t i = 0; i < replay.breakups.size(); ++i) {
				const auto& event = replay.breakups[i];
				const auto& saved = breakups[i];
				if (event.parent.id != saved.at("fragment_id") ||
				    std::abs(event.parent.timeS - saved.at("time_s").get<double>()) > 1e-6 ||
				    std::abs(event.heightWgs84M - saved.at("height_wgs84_m").get<double>()) > .01)
					throw std::runtime_error("Replay breakup event mismatch");
				tracks.at(event.parent.id).metadata["end_event"] = "breakup";
				const auto state = sample(event.parent.state, event.parent.timeS, origin, environment, p.windScale);
				const auto air = environment.air(geodetic(event.parent.state.positionEcefM));
				const double pressure = air.densityKgM3 * std::pow(state[7] * 1000, 2) *
				                        (spec.physics.ramPressure == RamPressure::RhoV2 ? 1. : .5);
				events.push_back({{"fragment_id", event.parent.id}, {"state", state}, {"pressure_pa", pressure}});
			}
			std::vector<PreviewTrack> candidates;
			for (const auto& [id, track] : tracks)
				candidates.push_back({id, track.metadata.at("parent_id"), track.last[2]});
			const auto chosen = previewTracks(candidates, traceLimit);
			events.erase(
			    std::remove_if(events.begin(), events.end(),
			                   [&](const auto& event) {
				                   return !chosen.contains(event.at("fragment_id").template get<std::uint64_t>());
			                   }),
			    events.end());
			Json traces = Json::array();
			std::vector<double> values;
			for (auto& [id, track] : tracks) {
				if (!chosen.contains(id))
					continue;
				if (track.samples.back() != track.last)
					track.samples.push_back(track.last);
				const auto rows = compact(std::move(track.samples));
				track.metadata["offset"] = values.size();
				track.metadata["count"] = rows.size();
				if (!track.metadata.contains("end_event"))
					track.metadata["end_event"] = "mass_cutoff";
				track.metadata["first"] = rows.front();
				track.metadata["last"] = rows.back();
				for (const auto& row : rows)
					values.insert(values.end(), row.begin(), row.end());
				traces.push_back(std::move(track.metadata));
			}
			Json wind = Json::array();
			const double top = std::min(
			    50000., environment.atmosphereGrid().firstHeightM +
			                (environment.atmosphereGrid().heightCount - 1) * environment.atmosphereGrid().heightStepM);
			for (double h = 0; h <= top; h += 250) {
				auto pos = origin;
				pos.heightWgs84M = h;
				const auto air = environment.air(pos);
				wind.push_back({h / 1000, air.windEnuMps[0] * p.windScale, air.windEnuMps[1] * p.windScale,
				                air.windEnuMps[2] * p.windScale});
			}
			return {{"verified", true},
			        {"max_impact_delta_m", maxDelta},
			        {"traces", std::move(traces)},
			        {"trace_selection",
			         {{"total", tracks.size()},
			          {"shown", chosen.size()},
			          {"limit", traceLimit},
			          {"sampled", chosen.size() != tracks.size()},
			          {"method", "terminal mass ranks with complete ancestry"}}},
			        {"samples", encodeFloats(values)},
			        {"breakups", std::move(events)},
			        {"wind", std::move(wind)}};
		}
	}  // namespace

	std::string base64(std::string_view bytes)
	{
		constexpr char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
		std::string result;
		result.reserve((bytes.size() + 2) / 3 * 4);
		for (std::size_t i = 0; i < bytes.size(); i += 3) {
			const unsigned value = (static_cast<unsigned char>(bytes[i]) << 16) |
			                       (i + 1 < bytes.size() ? static_cast<unsigned char>(bytes[i + 1]) << 8 : 0) |
			                       (i + 2 < bytes.size() ? static_cast<unsigned char>(bytes[i + 2]) : 0);
			result += alphabet[value >> 18];
			result += alphabet[(value >> 12) & 63];
			result += i + 1 < bytes.size() ? alphabet[(value >> 6) & 63] : '=';
			result += i + 2 < bytes.size() ? alphabet[value & 63] : '=';
		}
		return result;
	}

	std::string encodeFloats(std::span<const double> values)
	{
		std::string bytes(values.size() * 4, '\0');
		for (std::size_t i = 0; i < values.size(); ++i) {
			const auto bits = std::bit_cast<std::uint32_t>(static_cast<float>(values[i]));
			for (unsigned byte = 0; byte < 4; ++byte)
				bytes[byte * values.size() + i] = static_cast<char>((bits >> (byte * 8)) & 255);
		}
		return encodeText(bytes);
	}

	std::string encodeText(std::string_view bytes)
	{
		std::size_t size = 0;
		std::unique_ptr<void, decltype(&CPLFree)> compressed(
		    CPLZLibDeflate(bytes.data(), bytes.size(), 6, nullptr, 0, &size), CPLFree);
		if (!compressed)
			throw std::runtime_error("Cannot compress report samples");
		return base64(std::string_view(static_cast<const char*>(compressed.get()), size));
	}

	Json buildReportData(const std::filesystem::path& input, const ArchivedField& archive, ArtifactReport& report,
	                     std::stop_token cancel)
	{
		auto runs = readLines(input / "runs.jsonl", cancel);
		if (runs.size() != archive.summary.at("completed").get<std::size_t>())
			throw std::runtime_error("Realization archive count mismatch");
		std::map<std::uint64_t, std::size_t> indices;
		for (std::size_t i = 0; i < runs.size(); ++i) {
			const auto id = runs[i].at("run_id").get<std::uint64_t>();
			if (!indices.emplace(id, i).second || runs[i].at("status") != "complete")
				throw std::runtime_error("Invalid archived realization");
			runs[i]["breakups"] = Json::array();
			runs[i]["ground_count"] = 0.;
			runs[i]["tracked_count"] = 0;
			runs[i]["largest_mass_kg"] = 0.;
			runs[i]["first_breakup_km"] = nullptr;
		}
		for (auto& event : readLines(input / "breakups.jsonl", cancel)) {
			auto& run = runs.at(indices.at(event.at("run_id")));
			if (run["breakups"].empty() ||
			    event.at("time_s").get<double>() < run["breakups"][0]["time_s"].get<double>())
				run["first_breakup_km"] = event.at("height_wgs84_m").get<double>() / 1000;
			run["breakups"].push_back(std::move(event));
		}
		std::vector<std::map<std::uint64_t, std::size_t>> expected(runs.size());
		for (std::size_t i = 0; i < archive.points.size(); ++i) {
			const auto index = indices.at(archive.impactDetails[i].runId);
			auto& run = runs[index];
			const auto& point = archive.points[i];
			if (!expected[index].emplace(archive.impactDetails[i].fragmentId, i).second)
				throw std::runtime_error("Duplicate impact fragment ID");
			run["ground_count"] = run["ground_count"].get<double>() + point.multiplicity;
			run["tracked_count"] = run["tracked_count"].get<std::size_t>() + 1;
			run["largest_mass_kg"] = std::max(run["largest_mass_kg"].get<double>(), point.massKg);
		}
		const auto spec = loadJob(input / "job.json");
		const auto environment = jobfiles::readEnvironment(jobfiles::readJson(input / "environment.json"));
		std::map<std::uint64_t, EntryState> entries;
		std::vector<EntryState> sourceEntries;
		const auto entryFile = jobfiles::readJson(input / "entries.json");
		for (const auto& row : entryFile.at("entries")) {
			const auto entry = jobfiles::readEntry(row, 0);
			entries.emplace(entry.id, entry);
			sourceEntries.push_back(entry);
		}
		const auto entrySummary = summarizeEntries(sourceEntries);
		std::atomic_size_t next{};
		const auto workers =
		    std::min(runs.size(), std::max(std::size_t{1}, static_cast<std::size_t>(spec.monteCarlo.workers)));
		std::vector<std::jthread> threads;
		for (std::size_t worker = 0; worker < workers; ++worker)
			threads.emplace_back([&] {
				for (auto i = next.fetch_add(1); i < runs.size() && !cancel.stop_requested(); i = next.fetch_add(1)) {
					auto& run = runs[i];
					try {
						const auto& entry = entries.at(run.at("entry_id"));
						run["entry"] = jobfiles::entryJson(entry);
						run["flight"] = replayRun(run, run.at("breakups"), entry, spec, environment, archive,
						                          expected[i], cancel, previewLimit(runs.size()));
					} catch (const std::exception& error) {
						run["flight"] = {{"verified", false}, {"error", error.what()}};
					}
				}
			});
		threads.clear();
		if (cancel.stop_requested())
			throw std::runtime_error("Report cancelled");
		std::size_t failed = 0;
		for (const auto& run : runs)
			if (!run["flight"].value("verified", false))
				++failed;
		if (failed)
			report.warnings.push_back(std::to_string(failed) + " realization replays unavailable; see per-run errors");
		return {
		    {"runs", std::move(runs)},
		    {"entry_summary", entrySummary},
		    {"trace_encoding", "base64-deflate-byte-shuffle-float32-le"},
		    {"trace_columns",
		     {"time_s", "height_km", "mass_kg", "ecef_speed_km_s", "ground_range_km", "latitude_deg", "longitude_deg",
		      "air_speed_km_s"}},
		    {"trace_sampling",
		     "Visualization only: bounded terminal-mass-rank track selection with complete ancestry; up to 96 samples "
		     "per displayed track, birth/end retained. Ensemble statistics and KDE use the full archive. Displayed "
		     "track "
		     "counts are not population estimates."},
		    {"range_definition",
		     "Great-circle displacement from entry subpoint, radius 6367 km; not accumulated path length"}};
	}
}  // namespace strewn::output
