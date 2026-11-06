#include "FieldOutputInternal.h"
#include "ViewerAssets.h"
#include "ReportData.h"
#include <fstream>

namespace strewn::output {
	namespace {
		std::string escape(const std::string& text)
		{
			std::string result;
			for (char c : text) {
				switch (c) {
					case '&':
						result += "&amp;";
						break;
					case '<':
						result += "&lt;";
						break;
					case '>':
						result += "&gt;";
						break;
					case '"':
						result += "&quot;";
						break;
					default:
						result += c;
				}
			}
			return result;
		}

	}  // namespace

	void writeViewer(const std::filesystem::path& directory, const std::string& title, Json data)
	{
		std::vector<Json> flights;
		for (auto& run : data.at("runs")) {
			auto flight = std::move(run.at("flight"));
			run["breakup_count"] = run.at("breakups").size();
			double weighted = 0;
			for (const auto& event : run.at("breakups"))
				weighted += event.at("multiplicity").get<double>();
			run["breakup_weighted_count"] = weighted;
			flight["breakup_rows"] = std::move(run.at("breakups"));
			run.erase("breakups");
			run["flight"] = {{"verified", flight.value("verified", false)},
			                 {"element", "flight-data-" + std::to_string(flights.size())}};
			flights.push_back(std::move(flight));
		}
		const auto replace = [](std::string& text, const std::string& token, const std::string& value) {
			std::size_t pos = 0;
			while ((pos = text.find(token, pos)) != std::string::npos) {
				text.replace(pos, token.size(), value);
				pos += value.size();
			}
		};
		std::string html(reinterpret_cast<const char*>(assets::html), sizeof(assets::html)), json = data.dump();
		replace(json, "<", "\\u003c");
		replace(html, "__LEAFLET_CSS__", std::string(assets::css));
		replace(html, "__LEAFLET_JS__",
		        "/*\n" + std::string(assets::license) + "\n*/\n" +
		            std::string(reinterpret_cast<const char*>(assets::leaflet), sizeof(assets::leaflet)));
		replace(html, "__VIEWER_JS__",
		        std::string(reinterpret_cast<const char*>(assets::script), sizeof(assets::script)));
		replace(html, "__PLOTLY_JS__",
		        std::string(reinterpret_cast<const char*>(assets::plotly), sizeof(assets::plotly)));
		replace(html, "__REPORT_JS__", std::string(assets::reportScript));
		replace(html, "__PRINT_JS__", std::string(assets::printScript));
		replace(html, "__GRATICULE_JS__", std::string(assets::graticuleScript));
		replace(html, "__TRAJECTORY_JS__", std::string(assets::trajectoryScript));
		replace(html, "__MAP_FLIGHT_JS__", std::string(assets::mapFlightScript));
		replace(html, "__TITLE__", escape(title));
		replace(html, "__DATA__", json);
		const std::string marker = "<!-- STREWN_FLIGHT_PAYLOADS -->";
		const auto split = html.find(marker);
		if (split == std::string::npos)
			throw std::runtime_error("Missing trajectory payload placeholder");
		std::ofstream file(directory / "field.html", std::ios::binary);
		file.write(html.data(), static_cast<std::streamsize>(split));
		for (std::size_t i = 0; i < flights.size(); ++i) {
			auto payload = flights[i].dump();
			replace(payload, "<", "\\u003c");
			if (payload.size() > 65536)
				payload = Json{{"encoding", "base64-deflate-json"}, {"payload", encodeText(payload)}}.dump();
			file << "<script id=\"flight-data-" << i << "\" type=\"application/json\">" << payload << "</script>\n";
		}
		const auto tail = split + marker.size();
		file.write(html.data() + tail, static_cast<std::streamsize>(html.size() - tail));
		if (!file)
			throw std::runtime_error("Cannot write HTML field viewer");
	}
}  // namespace strewn::output
