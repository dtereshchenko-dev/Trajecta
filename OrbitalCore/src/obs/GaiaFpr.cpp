#include "orbital_core/obs/GaiaFpr.h"
#include "HttpClient.h"

#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <fstream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#ifdef _WIN32
#include <windows.h>
#endif

namespace orbital_core::obs {
	namespace {
		using Json = nlohmann::json;
		constexpr int kMaxRows = 50000;
		constexpr std::size_t kMaxBytes = 64 * 1024 * 1024;
		constexpr std::string_view kJobs = "https://gea.esac.esa.int/tap-server/tap/async/";
		constexpr double kTcbToTdb = 1.0 - 1.550519768e-8;
		constexpr std::array<std::string_view, 14> fields{"transit_id",
		                                                  "number_mp",
		                                                  "epoch_utc",
		                                                  "ra",
		                                                  "dec",
		                                                  "ra_error_systematic",
		                                                  "dec_error_systematic",
		                                                  "ra_dec_correlation_systematic",
		                                                  "ra_error_random",
		                                                  "dec_error_random",
		                                                  "ra_dec_correlation_random",
		                                                  "x_gaia_geocentric",
		                                                  "y_gaia_geocentric",
		                                                  "z_gaia_geocentric"};

		double number(const Json& value)
		{
			if (!value.is_number())
				throw std::runtime_error("Gaia FPR response contains a missing or invalid number");
			const double result = value.get<double>();
			if (!std::isfinite(result))
				throw std::runtime_error("Gaia FPR response contains a non-finite number");
			return result;
		}

		std::string transitId(const Json& value)
		{
			if (value.is_string())
				return value.get<std::string>();
			if (value.is_number_integer() || value.is_number_unsigned())
				return value.dump();
			throw std::runtime_error("Gaia FPR response has an invalid transit ID");
		}

		std::string readCache(const std::filesystem::path& file)
		{
			std::error_code error;
			const auto size = std::filesystem::file_size(file, error);
			if (error || size == 0 || size > kMaxBytes)
				throw std::runtime_error("Gaia FPR cache is unavailable");
			std::ifstream input(file, std::ios::binary);
			if (!input)
				throw std::runtime_error("Cannot open Gaia FPR cache");
			std::string data(static_cast<std::size_t>(size), '\0');
			if (!input.read(data.data(), static_cast<std::streamsize>(data.size())))
				throw std::runtime_error("Cannot read Gaia FPR cache");
			return data;
		}

		void writeCache(const std::filesystem::path& file, const std::string& data)
		{
			std::filesystem::create_directories(file.parent_path());
			static std::atomic<unsigned> serial{};
			const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
			auto temp = file;
			temp += "." + std::to_string(stamp) + "." + std::to_string(serial++) + ".tmp";
			try {
				std::ofstream out(temp, std::ios::binary | std::ios::trunc);
				if (!out.write(data.data(), static_cast<std::streamsize>(data.size())) || !out.flush())
					throw std::runtime_error("Cannot write Gaia FPR cache");
				out.close();
				if (!out)
					throw std::runtime_error("Cannot close Gaia FPR cache");
#ifdef _WIN32
				if (!MoveFileExW(temp.c_str(), file.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
					throw std::runtime_error("Cannot replace Gaia FPR cache");
#else
				std::filesystem::rename(temp, file);
#endif
			} catch (...) {
				std::error_code ignored;
				std::filesystem::remove(temp, ignored);
				throw;
			}
		}

		std::string fetchArchive(int number, std::stop_token cancel)
		{
			std::string query = "SELECT ";
			for (std::size_t i = 0; i < fields.size(); ++i) {
				if (i)
					query += ',';
				query += fields[i];
			}
			query += " FROM gaiafpr.sso_observation WHERE number_mp=" + std::to_string(number) + " ORDER BY epoch_utc";
			std::unique_ptr<char, decltype(&curl_free)> encoded(
			    curl_easy_escape(nullptr, query.c_str(), static_cast<int>(query.size())), &curl_free);
			if (!encoded)
				throw std::runtime_error("Cannot encode Gaia Archive query");
			const std::string form =
			    "REQUEST=doQuery&LANG=ADQL&FORMAT=json&PHASE=RUN&MAXREC=50000&QUERY=" + std::string(encoded.get());
			const auto response = trajecta::http::postFormForRedirect(std::string(kJobs.substr(0, kJobs.size() - 1)),
			                                                          form, 30, 65536, cancel);
			const auto& job = response.location;
			if (!job.starts_with(kJobs) || job.size() == kJobs.size() ||
			    job.find_first_not_of("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-",
			                          kJobs.size()) != std::string::npos)
				throw std::runtime_error("Invalid Gaia Archive job URL");
			const auto deadline = std::chrono::steady_clock::now() + std::chrono::minutes(5);
			while (std::chrono::steady_clock::now() < deadline) {
				if (cancel.stop_requested())
					throw trajecta::http::Cancelled{};
				auto phase = trajecta::http::get(job + "/phase", 20, 4096, cancel).body;
				const auto first = phase.find_first_not_of(" \t\r\n");
				const auto last = phase.find_last_not_of(" \t\r\n");
				phase = first == std::string::npos ? "" : phase.substr(first, last - first + 1);
				if (phase == "COMPLETED")
					return trajecta::http::get(job + "/results/result", 60, kMaxBytes, cancel, {}, 60).body;
				if (phase != "PENDING" && phase != "QUEUED" && phase != "EXECUTING")
					throw std::runtime_error("Gaia Archive job failed: " + phase);
				std::mutex mutex;
				std::condition_variable_any wake;
				std::unique_lock lock(mutex);
				wake.wait_for(lock, cancel, std::chrono::seconds(2), [] { return false; });
			}
			throw trajecta::http::Error("Gaia Archive job timed out", 0, true);
		}
	}  // namespace

	std::vector<GaiaFprObservation> parseGaiaFprResponse(std::string_view text, int minorPlanetNumber)
	{
		if (minorPlanetNumber <= 0)
			throw std::invalid_argument("Gaia FPR requires a numbered asteroid");
		const auto document = Json::parse(text.begin(), text.end());
		if (!document.is_object() || !document.contains("metadata") || !document["metadata"].is_array() ||
		    !document.contains("data") || !document["data"].is_array())
			throw std::runtime_error("Invalid Gaia FPR TAP response");
		const auto& metadata = document["metadata"];
		if (metadata.size() != fields.size())
			throw std::runtime_error("Unexpected Gaia FPR columns");
		for (std::size_t i = 0; i < fields.size(); ++i) {
			if (!metadata[i].is_object() || metadata[i].value("name", std::string{}) != fields[i])
				throw std::runtime_error("Unexpected Gaia FPR column order");
		}
		const auto& rows = document["data"];
		if (rows.empty())
			throw std::runtime_error("No Gaia FPR observations for this object");
		if (rows.size() >= kMaxRows)
			throw std::runtime_error("Gaia FPR result reached the download limit");
		std::unordered_map<std::string, std::size_t> transitCounts;
		for (const auto& row : rows) {
			if (!row.is_array() || row.size() != fields.size() || number(row[1]) != minorPlanetNumber)
				throw std::runtime_error("Gaia FPR result does not match the requested object");
			++transitCounts[transitId(row[0])];
		}
		std::vector<GaiaFprObservation> result;
		result.reserve(rows.size());
		for (const auto& row : rows) {
			const double sysRa = number(row[5]) / 1000.0, sysDec = number(row[6]) / 1000.0;
			const double randomRa = number(row[8]) / 1000.0, randomDec = number(row[9]) / 1000.0;
			const double cs = number(row[7]), cr = number(row[10]);
			if (sysRa < 0.0 || sysDec < 0.0 || randomRa < 0.0 || randomDec < 0.0 || std::abs(cs) > 1.0 ||
			    std::abs(cr) > 1.0)
				throw std::runtime_error("Invalid Gaia FPR uncertainty");
			// A conservative per-CCD approximation to transit-correlated systematics.
			const double count = static_cast<double>(transitCounts.at(transitId(row[0])));
			const double sigmaRa = std::sqrt(count * sysRa * sysRa + randomRa * randomRa);
			const double sigmaDec = std::sqrt(count * sysDec * sysDec + randomDec * randomDec);
			if (!(sigmaRa > 0.0) || !(sigmaDec > 0.0))
				throw std::runtime_error("Invalid Gaia FPR covariance");
			const double correlation = (count * cs * sysRa * sysDec + cr * randomRa * randomDec) / (sigmaRa * sigmaDec);
			const double epoch = number(row[2]) + 55197.0, ra = number(row[3]), dec = number(row[4]);
			if (!std::isfinite(epoch) || epoch < 0.0 || epoch > 1000000.0 || ra < 0.0 || ra >= 360.0 ||
			    std::abs(dec) > 90.0 || !std::isfinite(sigmaRa) || !std::isfinite(sigmaDec) ||
			    !std::isfinite(correlation) || std::abs(correlation) >= 1.0)
				throw std::runtime_error("Invalid Gaia FPR astrometry or covariance");
			result.push_back(
			    {epoch, ra, dec, sigmaRa, sigmaDec, correlation,
			     Vec3D(number(row[11]) * kTcbToTdb, number(row[12]) * kTcbToTdb, number(row[13]) * kTcbToTdb)});
		}
		return result;
	}

	GaiaFprDownload downloadGaiaFpr(int minorPlanetNumber, const std::filesystem::path& cacheDirectory,
	                                std::stop_token cancel)
	{
		if (minorPlanetNumber <= 0)
			throw std::invalid_argument("Gaia FPR requires a numbered asteroid");
		const auto cache = cacheDirectory.empty()
		                       ? std::filesystem::path{}
		                       : cacheDirectory / ("gaiafpr-" + std::to_string(minorPlanetNumber) + ".json");
		std::string data;
		try {
			data = fetchArchive(minorPlanetNumber, cancel);
		} catch (const trajecta::http::Error& error) {
			if (!error.unavailable() || cache.empty() || cancel.stop_requested())
				throw;
			try {
				return {parseGaiaFprResponse(readCache(cache), minorPlanetNumber), GaiaFprSource::Cache, error.what()};
			} catch (const std::exception&) {
				throw std::runtime_error(std::string(error.what()) + "; no valid Gaia FPR cache available");
			}
		}
		auto observations = parseGaiaFprResponse(data, minorPlanetNumber);
		if (cancel.stop_requested())
			throw trajecta::http::Cancelled{};
		std::string notice;
		if (!cache.empty()) {
			try {
				writeCache(cache, data);
			} catch (const std::exception& error) {
				notice = std::string("Gaia cache not saved: ") + error.what();
			}
		}
		return {std::move(observations), GaiaFprSource::Archive, std::move(notice)};
	}
}  // namespace orbital_core::obs
