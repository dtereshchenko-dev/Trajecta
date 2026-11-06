#pragma once

#include <cstddef>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <optional>
#include <cstdint>

namespace trajecta::http {

	struct Response {
		std::string body;
		std::string effectiveUrl;
		std::string location;
	};

	class Cancelled : public std::runtime_error {
	   public:
		Cancelled() : std::runtime_error("Download cancelled") {}
	};

	class Error : public std::runtime_error {
	   public:
		Error(std::string message, long status, bool unavailable)
		    : std::runtime_error(std::move(message)), m_status(status), m_unavailable(unavailable)
		{
		}

		long status() const { return m_status; }

		bool unavailable() const { return m_unavailable; }

	   private:
		long m_status;
		bool m_unavailable;
	};

	struct ByteRange {
		std::uint64_t first;
		std::uint64_t last;
	};

	Response get(const std::string& url, unsigned timeoutSeconds, std::size_t maxBytes, std::stop_token cancel,
	             std::optional<ByteRange> range = {}, unsigned lowSpeedSeconds = 15);

	Response postFormForRedirect(const std::string& url, const std::string& form, unsigned timeoutSeconds,
	                             std::size_t maxBytes, std::stop_token cancel);

}  // namespace trajecta::http
