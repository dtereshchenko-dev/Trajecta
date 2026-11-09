#include "HttpClient.h"
#include <curl/curl.h>
#include <array>
#include <exception>
#include <limits>
#include <memory>
#include <charconv>
#include <string_view>

namespace trajecta::http {
	namespace {

		class Runtime {
		   public:
			Runtime()
			{
				if (curl_global_init(CURL_GLOBAL_DEFAULT) != CURLE_OK)
					throw std::runtime_error("Cannot initialize libcurl");
			}

			~Runtime() { curl_global_cleanup(); }
		};

		struct Transfer {
			std::string body;
			std::size_t maxBytes;
			std::stop_token cancel;
			std::exception_ptr error;
			bool tooLarge{};
			std::optional<ByteRange> receivedRange;
			std::string location;
		};

		std::size_t header(char* data, std::size_t size, std::size_t count, void* opaque) noexcept
		{
			auto& transfer = *static_cast<Transfer*>(opaque);
			const std::string_view line(data, size * count);
			if (line.starts_with("HTTP/")) {
				transfer.receivedRange.reset();
				transfer.location.clear();
			}

			constexpr std::string_view locationName = "location:";
			bool isLocation = line.size() >= locationName.size();
			for (std::size_t i = 0; isLocation && i < locationName.size(); ++i) {
				const auto c = line[i];
				isLocation = (c >= 'A' && c <= 'Z' ? c + ('a' - 'A') : c) == locationName[i];
			}

			if (isLocation) {
				try {
					auto value = line.substr(locationName.size());
					const auto start = value.find_first_not_of(" \t");
					const auto end = value.find_last_not_of(" \t\r\n");
					transfer.location =
					    start == std::string_view::npos ? "" : std::string(value.substr(start, end - start + 1));
				}

				catch (...) {
					transfer.error = std::current_exception();
					return 0;
				}
			}

			constexpr std::string_view name = "content-range: bytes ";
			bool matches = line.size() >= name.size();
			for (std::size_t i = 0; matches && i < name.size(); ++i) {
				const auto c = line[i];
				matches = (c >= 'A' && c <= 'Z' ? c + ('a' - 'A') : c) == name[i];
			}

			if (matches) {
				ByteRange range{};
				const auto* end = line.data() + line.size();
				const auto first = std::from_chars(line.data() + name.size(), end, range.first);
				if (first.ec == std::errc{} && first.ptr < end && *first.ptr == '-') {
					const auto last = std::from_chars(first.ptr + 1, end, range.last);
					if (last.ec == std::errc{} && last.ptr < end && *last.ptr == '/')
						transfer.receivedRange = range;
				}
			}

			return size * count;
		}

		std::size_t receive(char* data, std::size_t size, std::size_t count, void* opaque) noexcept
		{
			auto& transfer = *static_cast<Transfer*>(opaque);
			if (transfer.cancel.stop_requested())
				return 0;

			if (size && count > (transfer.maxBytes - transfer.body.size()) / size) {
				transfer.tooLarge = true;
				return 0;
			}

			const auto bytes = size * count;
			try {
				transfer.body.append(data, bytes);
			}

			catch (...) {
				transfer.error = std::current_exception();
				return 0;
			}

			return bytes;
		}

		int progress(void* opaque, curl_off_t, curl_off_t, curl_off_t, curl_off_t) noexcept
		{
			return static_cast<Transfer*>(opaque)->cancel.stop_requested() ? 1 : 0;
		}

		template <class T>
		void option(CURL* handle, CURLoption name, T value)
		{
			if (curl_easy_setopt(handle, name, value) != CURLE_OK)
				throw std::runtime_error("Cannot configure libcurl transfer");
		}

	}  // namespace

	static Response request(const std::string& url, unsigned timeoutSeconds, std::size_t maxBytes,
	                        std::stop_token cancel, std::optional<ByteRange> range, unsigned lowSpeedSeconds,
	                        const std::string* form)
	{
		if (url.find('\0') != std::string::npos || (!url.starts_with("https://") && !url.starts_with("http://")))
			throw std::invalid_argument("Download URL must use HTTP or HTTPS");

		if (timeoutSeconds == 0 || timeoutSeconds > 3600 || lowSpeedSeconds == 0 || lowSpeedSeconds > 3600 ||
		    maxBytes == 0 || maxBytes > static_cast<std::size_t>((std::numeric_limits<curl_off_t>::max)()))
			throw std::invalid_argument("Invalid download limits");

		if (range && (range->last < range->first || range->last - range->first >= maxBytes))
			throw std::invalid_argument("Invalid download byte range");

		if (cancel.stop_requested())
			throw Cancelled{};

		static const Runtime runtime;
		std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> handle(curl_easy_init(), curl_easy_cleanup);
		if (!handle)
			throw std::runtime_error("Cannot create libcurl transfer");

		Transfer transfer{{}, maxBytes, cancel, {}, false, {}, {}};
		std::array<char, CURL_ERROR_SIZE> error{};
		option(handle.get(), CURLOPT_URL, url.c_str());
		option(handle.get(), CURLOPT_PROTOCOLS_STR, "http,https");
		option(handle.get(), CURLOPT_REDIR_PROTOCOLS_STR, url.starts_with("https://") ? "https" : "http,https");
		option(handle.get(), CURLOPT_FOLLOWLOCATION, form ? 0L : 1L);

		if (form) {
			option(handle.get(), CURLOPT_POST, 1L);
			option(handle.get(), CURLOPT_POSTFIELDS, form->c_str());
			option(handle.get(), CURLOPT_POSTFIELDSIZE_LARGE, static_cast<curl_off_t>(form->size()));
		}

		option(handle.get(), CURLOPT_MAXREDIRS, 5L);
		option(handle.get(), CURLOPT_NOSIGNAL, 1L);
		option(handle.get(), CURLOPT_CONNECTTIMEOUT, 10L);
		option(handle.get(), CURLOPT_TIMEOUT, static_cast<long>(timeoutSeconds));
		option(handle.get(), CURLOPT_LOW_SPEED_LIMIT, 1L);
		option(handle.get(), CURLOPT_LOW_SPEED_TIME, static_cast<long>(lowSpeedSeconds));
		option(handle.get(), CURLOPT_SSL_VERIFYPEER, 1L);
		option(handle.get(), CURLOPT_SSL_VERIFYHOST, 2L);
#ifdef _WIN32
		option(handle.get(), CURLOPT_SSL_OPTIONS, static_cast<long>(CURLSSLOPT_NATIVE_CA));
		option(handle.get(), CURLOPT_PROXY_SSL_OPTIONS, static_cast<long>(CURLSSLOPT_NATIVE_CA));
#endif
		option(handle.get(), CURLOPT_FAILONERROR, 1L);
		option(handle.get(), CURLOPT_MAXFILESIZE_LARGE, static_cast<curl_off_t>(maxBytes));
		option(handle.get(), CURLOPT_USERAGENT, "Trajecta/" TRAJECTA_VERSION);
		option(handle.get(), CURLOPT_ERRORBUFFER, error.data());
		option(handle.get(), CURLOPT_WRITEFUNCTION, &receive);
		option(handle.get(), CURLOPT_WRITEDATA, &transfer);
		option(handle.get(), CURLOPT_HEADERFUNCTION, &header);
		option(handle.get(), CURLOPT_HEADERDATA, &transfer);
		const auto rangeText = range ? std::to_string(range->first) + "-" + std::to_string(range->last) : std::string{};

		if (range)
			option(handle.get(), CURLOPT_RANGE, rangeText.c_str());

		option(handle.get(), CURLOPT_XFERINFOFUNCTION, &progress);
		option(handle.get(), CURLOPT_XFERINFODATA, &transfer);
		option(handle.get(), CURLOPT_NOPROGRESS, 0L);
		const auto status = curl_easy_perform(handle.get());
		if (cancel.stop_requested())
			throw Cancelled{};

		if (transfer.error)
			std::rethrow_exception(transfer.error);

		if (transfer.tooLarge || status == CURLE_FILESIZE_EXCEEDED)
			throw std::runtime_error("Download exceeds max_download_bytes");

		long code = 0;
		curl_easy_getinfo(handle.get(), CURLINFO_RESPONSE_CODE, &code);
		if (status != CURLE_OK) {
			const bool unavailable = code == 404 || code == 403 || code == 429 || code >= 500 ||
			                         status == CURLE_COULDNT_RESOLVE_HOST || status == CURLE_COULDNT_RESOLVE_PROXY ||
			                         status == CURLE_COULDNT_CONNECT || status == CURLE_OPERATION_TIMEDOUT;
			throw Error(std::string("Download failed: ") + (error[0] ? error.data() : curl_easy_strerror(status)), code,
			            unavailable);
		}

		char* effective = nullptr;
		if (code != (form ? 303 : (range ? 206 : 200)))
			throw std::runtime_error("Unexpected HTTP status " + std::to_string(code));

		if (range &&
		    (!transfer.receivedRange || transfer.receivedRange->first != range->first ||
		     transfer.receivedRange->last != range->last || transfer.body.size() != range->last - range->first + 1))
			throw std::runtime_error("Server returned an incorrect or incomplete byte range");

		if (curl_easy_getinfo(handle.get(), CURLINFO_EFFECTIVE_URL, &effective) != CURLE_OK || !effective)
			throw std::runtime_error("Cannot read the final download URL");

		return {std::move(transfer.body), effective, std::move(transfer.location)};
	}

	Response get(const std::string& url, unsigned timeoutSeconds, std::size_t maxBytes, std::stop_token cancel,
	             std::optional<ByteRange> range, unsigned lowSpeedSeconds)
	{
		return request(url, timeoutSeconds, maxBytes, cancel, range, lowSpeedSeconds, nullptr);
	}

	Response postFormForRedirect(const std::string& url, const std::string& form, unsigned timeoutSeconds,
	                             std::size_t maxBytes, std::stop_token cancel)
	{
		return request(url, timeoutSeconds, maxBytes, cancel, {}, timeoutSeconds, &form);
	}

}  // namespace trajecta::http
