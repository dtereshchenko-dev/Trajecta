#include "orbital_core/obs/SbdbPhysical.h"
#include "HttpClient.h"

#include <nlohmann/json.hpp>
#include <charconv>
#include <cmath>
#include <stdexcept>
#include <string>

namespace orbital_core::obs {
	namespace {
		std::optional<double> numericValue(const nlohmann::json& value)
		{
			if (value.is_number())
				return value.get<double>();
			if (!value.is_string())
				return std::nullopt;
			const auto& text = value.get_ref<const std::string&>();
			double number{};
			const auto result = std::from_chars(text.data(), text.data() + text.size(), number);
			if (result.ec != std::errc{} || result.ptr != text.data() + text.size())
				return std::nullopt;
			return number;
		}
	}  // namespace

	std::optional<double> parseSbdbRadiusKm(std::string_view text, int minorPlanetNumber)
	{
		if (minorPlanetNumber <= 0)
			throw std::invalid_argument("SBDB requires a numbered asteroid");
		const auto data = nlohmann::json::parse(text.begin(), text.end());
		if (!data.is_object() || !data.contains("object") || !data["object"].is_object() ||
		    data["object"].value("des", std::string{}) != std::to_string(minorPlanetNumber))
			throw std::runtime_error("SBDB response does not match the requested asteroid");
		if (!data.contains("phys_par") || data["phys_par"].is_null())
			return std::nullopt;
		if (!data["phys_par"].is_array())
			throw std::runtime_error("Invalid SBDB physical parameters");
		for (const auto& parameter : data["phys_par"]) {
			if (!parameter.is_object() || parameter.value("name", std::string{}) != "diameter")
				continue;
			if (parameter.value("units", std::string{}) != "km" || !parameter.contains("value"))
				return std::nullopt;
			const auto diameter = numericValue(parameter["value"]);
			if (!diameter || !std::isfinite(*diameter) || *diameter <= 0.0)
				return std::nullopt;
			return *diameter / 2.0;
		}
		return std::nullopt;
	}

	std::optional<double> downloadSbdbRadiusKm(int minorPlanetNumber, std::stop_token cancel)
	{
		if (minorPlanetNumber <= 0)
			throw std::invalid_argument("SBDB requires a numbered asteroid");
		const auto url =
		    "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=" + std::to_string(minorPlanetNumber) + "&phys-par=1";
		return parseSbdbRadiusKm(trajecta::http::get(url, 15, 2 * 1024 * 1024, cancel).body, minorPlanetNumber);
	}
}  // namespace orbital_core::obs
