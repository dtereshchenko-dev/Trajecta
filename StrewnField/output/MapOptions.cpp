#include "FieldOutputInternal.h"
#include "../app/JobFiles.h"

namespace strewn::output {

	MapOptions readOptions(const Json& v, const std::filesystem::path& base)
	{
		using namespace jobfiles;
		keys(v, {"title", "method", "cell_m", "fixed_sigma_m", "adaptive_alpha", "adaptive_limits", "adaptive_classes",
		         "center_deg", "half_size_km", "image_width", "basemap", "offline", "font_file", "auto_cell_size"});
		MapOptions o;
		o.title = v.value("title", o.title);
		const auto method = v.value("method", "adaptive");
		if (method == "fixed")
			o.method = KdeMethod::Fixed;
		else if (method != "adaptive")
			throw std::invalid_argument("Map method must be adaptive or fixed");
		o.cellM = v.value("cell_m", o.cellM);
		o.autoCellSize = v.value("auto_cell_size", o.autoCellSize);
		o.fixedSigmaM = v.value("fixed_sigma_m", o.fixedSigmaM);
		o.adaptiveAlpha = v.value("adaptive_alpha", o.adaptiveAlpha);
		if (v.contains("adaptive_limits"))
			o.adaptiveLimits = numbers<2>(v.at("adaptive_limits"));
		const auto classes = count(v, "adaptive_classes", o.adaptiveClasses),
		           width = count(v, "image_width", o.imageWidth);
		if (classes > 32 || width > 4096)
			throw std::invalid_argument("Map dimension/class limit exceeded");
		o.adaptiveClasses = static_cast<unsigned>(classes);
		o.imageWidth = static_cast<unsigned>(width);
		if (v.contains("center_deg") && !v.at("center_deg").is_null())
			o.centerDeg = numbers<2>(v.at("center_deg"));
		if (v.contains("half_size_km"))
			o.halfSizeKm = numbers<2>(v.at("half_size_km"));
		const auto basemap = v.value("basemap", "none");
		if (basemap == "esri_imagery")
			o.basemap = Basemap::EsriImagery;
		else if (basemap != "none")
			throw std::invalid_argument("Map basemap must be none or esri_imagery");
		o.offline = v.value("offline", false);
		if (v.contains("font_file")) {
			const auto text = v.at("font_file").get<std::string>();
			o.fontFile =
			    std::filesystem::absolute(base / std::filesystem::path(std::u8string(text.begin(), text.end())))
			        .lexically_normal();
		}
		validate(o);
		return o;
	}

	Json optionsJson(const MapOptions& o)
	{
		Json value = {{"title", o.title},
		              {"method", o.method == KdeMethod::Fixed ? "fixed" : "adaptive"},
		              {"cell_m", o.cellM},
		              {"auto_cell_size", o.autoCellSize},
		              {"fixed_sigma_m", o.fixedSigmaM},
		              {"adaptive_alpha", o.adaptiveAlpha},
		              {"adaptive_limits", o.adaptiveLimits},
		              {"adaptive_classes", o.adaptiveClasses},
		              {"center_deg", o.centerDeg ? Json(*o.centerDeg) : Json(nullptr)},
		              {"half_size_km", o.halfSizeKm},
		              {"basemap", o.basemap == Basemap::None ? "none" : "esri_imagery"},
		              {"offline", o.offline}};
		return value;
	}

	const char* statusName(ArtifactStatus s)
	{
		switch (s) {
			case ArtifactStatus::Skipped:
				return "skipped";
			case ArtifactStatus::Complete:
				return "complete";
			case ArtifactStatus::Empty:
				return "empty";
			case ArtifactStatus::Failed:
				return "failed";
			case ArtifactStatus::Cancelled:
				return "cancelled";
		}
		throw std::invalid_argument("Unknown artifact status");
	}

	Json reportJson(const ArtifactReport& report)
	{
		Json files = Json::array();
		for (const auto& file : report.files)
			files.push_back(jobfiles::pathText(file));
		return {{"status", statusName(report.status)},
		        {"files", files},
		        {"warnings", report.warnings},
		        {"error", report.error},
		        {"elapsed_s", report.elapsedS}};
	}

}  // namespace strewn::output
