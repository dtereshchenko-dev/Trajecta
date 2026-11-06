#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include "orbital_core/TimeScale.h"
#include "IntervalMcAssets.h"
#include "hazard/McBoundaryScan.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <cmath>

namespace orbital_core::hazard {
	namespace {
		using Json = nlohmann::json;

		template <std::size_t N>
		std::string asset(const unsigned char (&value)[N])
		{
			return {reinterpret_cast<const char*>(value), N};
		}

		void replace(std::string& text, const std::string& key, const std::string& value)
		{
			std::size_t position = 0;
			while ((position = text.find(key, position)) != std::string::npos) {
				text.replace(position, key.size(), value);
				position += value.size();
			}
		}

		Json utc(const Json& epoch)
		{
			if (!epoch.is_number())
				return nullptr;
			const double value = time::fromMjdTdb(epoch.get<double>(), time::Scale::UTC);
			return std::isfinite(value) ? Json(value) : Json(nullptr);
		}

		void addUtc(Json& visit)
		{
			for (const auto* key : {"start", "end"})
				visit[std::string(key) + "_utc"] = utc(visit.at(key));
		}
	}  // namespace

	void renderIntervalMonteCarloReport(const std::filesystem::path& directory)
	{
		std::ifstream input(directory / "result.json", std::ios::binary);
		if (!input)
			throw std::runtime_error("Cannot read MC numerical archive");
		Json data = Json::parse(input);
		const bool legacy = data.at("schema") == "orbital-interval-mc-1";
		const bool tdbArchive = data.at("schema") == "orbital-interval-mc-4";
		if (!legacy && data.at("schema") != "orbital-interval-mc-2" && data.at("schema") != "orbital-interval-mc-3" &&
		    !tdbArchive)
			throw std::runtime_error("Unsupported MC archive schema");
		if (data.at("time_scale") != (tdbArchive ? "MJD TDB" : "MJD TT"))
			throw std::runtime_error("Unsupported MC archive time scale");
		const auto sourceScale = tdbArchive ? time::Scale::TDB : time::Scale::TT;
		const auto convertEpoch = [sourceScale](Json& epoch) {
			if (epoch.is_number())
				epoch = time::toMjdTdb(epoch.get<double>(), sourceScale);
		};
		for (const auto* key : {"start", "end"})
			convertEpoch(data.at(key));
		const auto convertClone = [&](Json& clone) {
			if (clone.contains("epoch"))
				convertEpoch(clone.at("epoch"));
			if (legacy) {
				for (auto& visit : clone.at("encounters")) {
					for (const auto* key : {"start", "end", "closest"})
						if (visit.contains(key))
							convertEpoch(visit.at(key));
					if (visit.contains("entry") && visit.at("entry").is_object())
						convertEpoch(visit.at("entry").at("epoch"));
				}
			} else if (clone.at("entry").is_object())
				convertEpoch(clone.at("entry").at("epoch"));
		};
		convertClone(data.at("nominal"));
		for (auto& clone : data.at("clones"))
			convertClone(clone);
		data["source_time_scale"] = data.at("time_scale");
		data["time_scale"] = "MJD TDB";
		const int target = data.at("target").get<int>();
		if (target < 0 || target >= static_cast<int>(Planet::Count_))
			throw std::runtime_error("Invalid MC target");
		data["target_name"] = longLabel(static_cast<Planet>(target));
		for (const auto* key : {"start", "end"})
			data[std::string(key) + "_utc"] = utc(data.at(key));
		Json entries = Json::array();
		IntervalMcResult summary;
		auto collect = [&](const Json& clone, bool nominal) {
			McClone trial;
			trial.index = clone.at("index");
			trial.status = static_cast<McCloneStatus>(clone.at("status").get<int>());
			trial.entryBeforeInterval = clone.value("entry_before_interval", false);
			Json entry = nullptr;
			if (trial.status == McCloneStatus::Complete) {
				if (legacy) {
					for (const auto& visit : clone.at("encounters"))
						if (visit.contains("entry")) {
							entry = visit.at("entry");
							break;
						}
				} else
					entry = clone.at("entry");
			}
			if (!entry.is_null()) {
				entry["clone"] = clone.at("index");
				entry["epoch_utc"] = utc(entry.at("epoch"));
				trial.entry = McBoundaryEntry{.mjdTdb = entry.at("epoch").get<double>()};
				entries.push_back(std::move(entry));
			}
			if (!nominal)
				summary.clones.push_back(std::move(trial));
		};
		collect(data.at("nominal"), true);
		for (const auto& clone : data.at("clones"))
			collect(clone, false);
		mc::summarize(summary);
		data["events"] = Json::array();
		for (const auto& event : summary.events) {
			Json value = {{"start", event.startMjdTdb},
			              {"end", event.endMjdTdb},
			              {"entries", event.entryClones},
			              {"probability", event.probability},
			              {"probability_low", event.probabilityLow},
			              {"probability_high", event.probabilityHigh}};
			addUtc(value);
			data["events"].push_back(std::move(value));
		}
		for (auto& entry : entries) {
			entry["event"] = nullptr;
			const auto epoch = entry.at("epoch").get<double>();
			for (std::size_t i = 0; i < summary.events.size(); ++i)
				if (epoch >= summary.events[i].startMjdTdb && epoch <= summary.events[i].endMjdTdb) {
					entry["event"] = i;
					break;
				}
		}
		data.erase("clones");
		data.erase("reporting_radius_au");
		data["event_grouping"] = "consecutive crossing gaps; presentation only";
		data["event_max_gap_days"] = mc::eventMaxGapDays;
		data["entries"] = std::move(entries);
		std::string html = asset(report_assets::html);
		replace(html, "__LEAFLET_CSS__", std::string(report_assets::leafletCss));
		replace(html, "__REPORT_CSS__", asset(report_assets::style));
		replace(html, "__LEAFLET_JS__", asset(report_assets::leaflet));
		replace(html, "__GRATICULE_JS__", asset(report_assets::graticule));
		replace(html, "__VIEWER_JS__", asset(report_assets::script));
		replace(html, "__LICENSES__", std::string(report_assets::license));
		auto payload = data.dump();
		replace(payload, "<", "\\u003c");
		replace(html, "__DATA__", payload);
		std::ofstream file(directory / "report.html", std::ios::binary);
		file.exceptions(std::ios::badbit | std::ios::failbit);
		file << html;
		file.close();
	}
}  // namespace orbital_core::hazard
