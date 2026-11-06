#include "EnvironmentLoader.h"
#include "HttpClient.h"
#include <chrono>
#include <iomanip>
#include <sstream>

namespace strewn::jobfiles {
	namespace {

		std::string cacheKey(const std::string& url)
		{
			std::uint64_t hash = 14695981039346656037ULL;
			for (const unsigned char byte : url)
				hash = (hash ^ byte) * 1099511628211ULL;
			std::ostringstream out;
			out << std::hex << std::setfill('0') << std::setw(16) << hash;
			return out.str();
		}

		class CacheLock {
		   public:
			explicit CacheLock(std::filesystem::path path) : m_path(std::move(path))
			{
				if (!std::filesystem::create_directory(m_path))
					throw std::runtime_error("Environment cache is busy: " + m_path.string());
			}

			~CacheLock()
			{
				std::error_code ignored;
				std::filesystem::remove(m_path, ignored);
			}

			CacheLock(const CacheLock&) = delete;
			CacheLock& operator=(const CacheLock&) = delete;

		   private:
			std::filesystem::path m_path;
		};

		void validatePart(const Json& document, const char* part)
		{
			auto combined = environmentJson(Environment::idealized(0));
			combined["schema_version"] = document.at("schema_version");
			combined["height_datum"] = document.at("height_datum");
			combined[part] = document.at(part);
			(void)readEnvironment(combined);
		}

		Json acquire(const std::string& url, const char* part, const std::filesystem::path& directory,
		             const EnvironmentUrls& options, std::stop_token cancel)
		{
			if (url.size() > 8192 || url.find('\0') != std::string::npos ||
			    (!url.starts_with("http://") && !url.starts_with("https://")))
				throw std::invalid_argument("Environment URL must use HTTP or HTTPS and contain at most 8192 bytes");
			if (cancel.stop_requested())
				throw trajecta::http::Cancelled{};
			std::filesystem::create_directories(directory);
			const auto key = cacheKey(url);
			CacheLock lock(directory / (key + ".lock"));
			const auto file = directory / (key + ".json");
			if (std::filesystem::exists(file)) {
				if (std::filesystem::file_size(file) > options.maxDownloadBytes * 2ULL + 32768)
					throw std::runtime_error("Cached environment exceeds download limits");
				auto cached = readJson(file);
				if (cached.at("source_url") != url)
					throw std::runtime_error("Environment cache key collision");
				validatePart(cached.at("document"), part);
				return cached;
			}
			const auto response = trajecta::http::get(url, options.timeoutSeconds, options.maxDownloadBytes, cancel);
			const auto document = Json::parse(response.body);
			validatePart(document, part);
			if (cancel.stop_requested())
				throw trajecta::http::Cancelled{};
			Json cached = {{"source_url", url},
			               {"effective_url", response.effectiveUrl},
			               {"downloaded_unix_seconds", std::chrono::duration_cast<std::chrono::seconds>(
			                                               std::chrono::system_clock::now().time_since_epoch())
			                                               .count()},
			               {"document", document}};
			const auto temporary = directory / (key + ".part");
			try {
				const auto text = cached.dump();
				if (text.size() > options.maxDownloadBytes * 2ULL + 32768)
					throw std::runtime_error("Serialized environment cache exceeds download limits");
				{
					std::ofstream out(temporary, std::ios::binary);
					out.exceptions(std::ios::badbit | std::ios::failbit);
					out << text;
				}
				std::filesystem::rename(temporary, file);
			} catch (...) {
				std::error_code ignored;
				std::filesystem::remove(temporary, ignored);
				throw;
			}
			return cached;
		}

		Json metadata(Json source)
		{
			source["document"].erase("atmosphere");
			source["document"].erase("terrain");
			return source;
		}

	}  // namespace

	Json downloadEnvironment(const EnvironmentUrls& sources, const std::filesystem::path& directory,
	                         std::stop_token cancel)
	{
		const auto weather = acquire(sources.weather, "atmosphere", directory / "weather", sources, cancel);
		const auto terrain = acquire(sources.terrain, "terrain", directory / "terrain", sources, cancel);
		return {{"schema_version", 1},
		        {"height_datum", "WGS84"},
		        {"atmosphere", weather.at("document").at("atmosphere")},
		        {"terrain", terrain.at("document").at("terrain")},
		        {"download_sources", {{"weather", metadata(weather)}, {"terrain", metadata(terrain)}}}};
	}

}  // namespace strewn::jobfiles
