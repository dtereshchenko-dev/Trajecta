#include "ReportData.h"
#include "../app/SourceCache.h"
#include <algorithm>
#include <cmath>
#include <numbers>

namespace strewn::output {

	Json embedBasemap(const Json& bounds, const MapOptions& options, const std::filesystem::path& cache,
	                  ArtifactReport& report, std::stop_token cancel)
	{
		Json result = Json::array();
		if (options.basemap == Basemap::None || bounds.is_null())
			return result;
		const double west = bounds[0][1], east = bounds[1][1], south = bounds[0][0], north = bounds[1][0];
		const auto tileY = [](double latitude, double scale) {
			return (1 -
			        std::asinh(std::tan(std::clamp(latitude, -85., 85.) * std::numbers::pi / 180)) / std::numbers::pi) *
			       scale / 2;
		};
		int zoom = std::clamp(static_cast<int>(std::ceil(std::log2(1600 * 360.0 / (256 * (east - west))))), 1, 15);
		int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
		do {
			const double scale = std::ldexp(1.0, zoom);
			x0 = static_cast<int>(std::floor((west + 180) / 360 * scale));
			x1 = static_cast<int>(std::floor((east + 180) / 360 * scale));
			y0 = static_cast<int>(std::floor(tileY(north, scale)));
			y1 = static_cast<int>(std::floor(tileY(south, scale)));
			if ((x1 - x0 + 1) * (y1 - y0 + 1) <= 64)
				break;
			--zoom;
		} while (zoom > 0);
		AutomaticEnvironment network;
		network.offline = options.offline;
		network.timeoutSeconds = 15;
		network.maxDownloadBytes = 2 * 1024 * 1024;
		try {
			for (int y = y0; y <= y1; ++y)
				for (int x = x0; x <= x1; ++x) {
					if (cancel.stop_requested())
						throw trajecta::http::Cancelled{};
					const int wrapped = (x % (1 << zoom) + (1 << zoom)) % (1 << zoom);
					const auto url =
					    "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/" +
					    std::to_string(zoom) + "/" + std::to_string(y) + "/" + std::to_string(wrapped) +
					    "?blankTile=false";
					const auto bytes =
					    sources::readBytes(sources::fetch(url, cache, network, cancel), network.maxDownloadBytes);
					const bool jpeg = bytes.size() > 2 && static_cast<unsigned char>(bytes[0]) == 255 &&
					                  static_cast<unsigned char>(bytes[1]) == 216;
					const bool png = bytes.size() > 8 && bytes.substr(0, 8) == std::string("\x89PNG\r\n\x1a\n", 8);
					if (!jpeg && !png)
						throw std::runtime_error("Invalid imagery tile signature");
					result.push_back({{"x", x},
					                  {"y", y},
					                  {"z", zoom},
					                  {"url", std::string(jpeg ? "data:image/jpeg;base64," : "data:image/png;base64,") +
					                              base64(bytes)}});
				}
		} catch (const sources::Unavailable& error) {
			report.warnings.push_back(
			    std::string("Satellite imagery partially/unavailable; geographic layers remain usable: ") +
			    error.what());
		}
		return result;
	}
}  // namespace strewn::output
