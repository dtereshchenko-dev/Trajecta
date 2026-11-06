#pragma once

#include "HttpClient.h"
#include "strewn_core/StrewnJob.h"
#include <nlohmann/json.hpp>

namespace strewn::sources {

	using Json = nlohmann::json;

	class Unavailable : public std::runtime_error {
	   public:
		using std::runtime_error::runtime_error;
	};

	std::string fingerprint(std::string_view text);
	std::string readBytes(const std::filesystem::path& path, std::size_t limit);
	std::filesystem::path fetch(const std::string& url, const std::filesystem::path& directory,
	                            const AutomaticEnvironment& options, std::stop_token cancel,
	                            std::optional<trajecta::http::ByteRange> range = {});

}  // namespace strewn::sources
