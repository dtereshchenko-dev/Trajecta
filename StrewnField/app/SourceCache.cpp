#include "SourceCache.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace strewn::sources {

	std::string fingerprint(std::string_view text)
	{
		std::uint64_t hash = 14695981039346656037ULL;
		for (const unsigned char byte : text)
			hash = (hash ^ byte) * 1099511628211ULL;
		std::ostringstream out;
		out << std::hex << std::setfill('0') << std::setw(16) << hash;
		return out.str();
	}

	std::string readBytes(const std::filesystem::path& path, std::size_t limit)
	{
		const auto size = std::filesystem::file_size(path);
		if (size > limit)
			throw std::runtime_error("Cached source exceeds max_download_bytes");
		std::string data(static_cast<std::size_t>(size), '\0');
		std::ifstream in(path, std::ios::binary);
		in.exceptions(std::ios::failbit | std::ios::badbit);
		in.read(data.data(), static_cast<std::streamsize>(data.size()));
		return data;
	}

	std::filesystem::path fetch(const std::string& url, const std::filesystem::path& directory,
	                            const AutomaticEnvironment& options, std::stop_token cancel,
	                            std::optional<trajecta::http::ByteRange> range)
	{
		if (cancel.stop_requested())
			throw trajecta::http::Cancelled{};
		const auto identity =
		    url + (range ? "#bytes=" + std::to_string(range->first) + "-" + std::to_string(range->last) : "");
		const auto key = fingerprint(identity);
		std::filesystem::create_directories(directory);
		const auto lock = directory / (key + ".lock");
		if (!std::filesystem::create_directory(lock))
			throw std::runtime_error("Source cache is busy: " + lock.string());

		struct Unlock {
			std::filesystem::path path;

			~Unlock()
			{
				std::error_code error;
				std::filesystem::remove(path, error);
			}
		} unlock{lock};

		const auto file = directory / (key + ".bin"), meta = directory / (key + ".json");
		if (std::filesystem::exists(meta)) {
			const auto record = Json::parse(readBytes(meta, 32768));
			const auto bytes = readBytes(file, options.maxDownloadBytes);
			if (record.at("identity") != identity || record.at("bytes") != bytes.size() ||
			    record.at("fnv1a64") != fingerprint(bytes))
				throw std::runtime_error("Source cache integrity check failed");
			return file;
		}
		if (options.offline)
			throw Unavailable("Offline cache miss: " + url);
		trajecta::http::Response response;
		try {
			response = trajecta::http::get(url, options.timeoutSeconds, options.maxDownloadBytes, cancel, range);
		} catch (const trajecta::http::Error& error) {
			if (error.unavailable())
				throw Unavailable(url + ": " + error.what());
			throw;
		}
		const Json record = {{"identity", identity},
		                     {"source_url", url},
		                     {"effective_url", response.effectiveUrl},
		                     {"bytes", response.body.size()},
		                     {"fnv1a64", fingerprint(response.body)},
		                     {"downloaded_unix_seconds", std::chrono::duration_cast<std::chrono::seconds>(
		                                                     std::chrono::system_clock::now().time_since_epoch())
		                                                     .count()}};
		const auto part = directory / (key + ".part"), metaPart = directory / (key + ".json.part");
		try {
			{
				std::ofstream out(part, std::ios::binary);
				out.exceptions(std::ios::failbit | std::ios::badbit);
				out.write(response.body.data(), static_cast<std::streamsize>(response.body.size()));
			}
			{
				std::ofstream out(metaPart);
				out.exceptions(std::ios::failbit | std::ios::badbit);
				out << record.dump(2);
			}
			if (cancel.stop_requested())
				throw trajecta::http::Cancelled{};
			std::filesystem::remove(file);
			std::filesystem::rename(part, file);
			std::filesystem::rename(metaPart, meta);
		} catch (...) {
			std::error_code error;
			std::filesystem::remove(part, error);
			std::filesystem::remove(metaPart, error);
			throw;
		}
		return file;
	}

}  // namespace strewn::sources
