#include "strewn_core/Environment.h"
#include "VectorMath.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace strewn {
	namespace {

		void validateAxis(const std::vector<double>& values)
		{
			if (values.empty())
				throw std::invalid_argument("Environment axis is empty");
			for (std::size_t i = 0; i < values.size(); ++i) {
				if (!std::isfinite(values[i]) || (i && values[i] <= values[i - 1]))
					throw std::invalid_argument("Environment axes must be finite and increasing");
				if (i > 1 && std::abs((values[i] - values[i - 1]) - (values[1] - values[0])) > 1e-9)
					throw std::invalid_argument("Environment grid must be regular");
			}
		}

		struct Interval {
			std::size_t low;
			std::size_t high;
			double fraction;
		};

		Interval locate(const std::vector<double>& axis, double value)
		{
			if (!std::isfinite(value))
				throw std::invalid_argument("Non-finite environment coordinate");
			if (axis.size() == 1)
				return {0, 0, 0.0};
			if (value < axis.front() - 1e-9 || value > axis.back() + 1e-9)
				throw std::out_of_range("Trajectory left the prepared environment grid");
			const double index =
			    std::clamp((value - axis[0]) / (axis[1] - axis[0]), 0.0, static_cast<double>(axis.size() - 1));
			const auto low = std::min(static_cast<std::size_t>(index), axis.size() - 2);
			return {low, low + 1, index - static_cast<double>(low)};
		}

		AirSample mix(const AirSample& a, const AirSample& b, double fraction)
		{
			AirSample result{a.densityKgM3 + fraction * (b.densityKgM3 - a.densityKgM3),
			                 a.temperatureK + fraction * (b.temperatureK - a.temperatureK),
			                 {}};
			for (std::size_t i = 0; i < 3; ++i)
				result.windEnuMps[i] = a.windEnuMps[i] + fraction * (b.windEnuMps[i] - a.windEnuMps[i]);
			return result;
		}

		std::size_t checkedSize(std::size_t a, std::size_t b, std::size_t c = 1)
		{
			if (!a || !b || !c || a > std::numeric_limits<std::size_t>::max() / b ||
			    a * b > std::numeric_limits<std::size_t>::max() / c)
				throw std::invalid_argument("Invalid environment dimensions");
			return a * b * c;
		}

	}  // namespace

	GeodeticPosition geodetic(const Vector3& r)
	{
		const auto [x, y, z] = r;
		const double p = std::hypot(x, y);
		if (!detail::finite(r) || detail::norm(r) < 1.0)
			throw std::invalid_argument("Invalid ECEF position");
		if (p < 1e-7)
			return {std::copysign(detail::pi / 2.0, z), 0.0, std::abs(z) - detail::earthA * (1.0 - detail::earthF)};
		double lat = std::atan2(z, p * (1.0 - detail::earthE2));
		for (int i = 0; i < 3; ++i) {
			const double s = std::sin(lat);
			const double n = detail::earthA / std::sqrt(1.0 - detail::earthE2 * s * s);
			lat = std::atan2(z, p * (1.0 - detail::earthE2 * n / (n + p / std::cos(lat) - n)));
		}
		const double s = std::sin(lat);
		const double n = detail::earthA / std::sqrt(1.0 - detail::earthE2 * s * s);
		return {lat, std::atan2(y, x), p / std::cos(lat) - n};
	}

	Vector3 ecefPosition(const GeodeticPosition& position)
	{
		const auto [lat, lon, h] = position;
		if (!std::isfinite(lat) || !std::isfinite(lon) || !std::isfinite(h) || std::abs(lat) > detail::pi / 2.0 ||
		    h < -10000.0)
			throw std::invalid_argument("Invalid WGS84 position");
		const double s = std::sin(lat), c = std::cos(lat);
		const double n = detail::earthA / std::sqrt(1.0 - detail::earthE2 * s * s);
		return {(n + h) * c * std::cos(lon), (n + h) * c * std::sin(lon), (n * (1.0 - detail::earthE2) + h) * s};
	}

	Vector3 enuToEcef(const Vector3& v, double lat, double lon)
	{
		const double sl = std::sin(lat), cl = std::cos(lat), so = std::sin(lon), co = std::cos(lon);
		return {-so * v[0] - sl * co * v[1] + cl * co * v[2], co * v[0] - sl * so * v[1] + cl * so * v[2],
		        cl * v[1] + sl * v[2]};
	}

	AirSample standardAtmosphere(double heightM)
	{
		if (!std::isfinite(heightM) || heightM < -1000.0 || heightM > 140000.0)
			throw std::out_of_range("Standard atmosphere supports heights -1000..140000 m");
		const double z = heightM / 1000.0;
		if (z < 86.0) {
			constexpr double layers[][4] = {{0.0, 288.15, -6.5, 101325.0}, {11.0, 216.65, 0.0, 22632.06},
			                                {20.0, 216.65, 1.0, 5474.889}, {32.0, 228.65, 2.8, 868.0187},
			                                {47.0, 270.65, 0.0, 110.9063}, {51.0, 270.65, -2.8, 66.93887},
			                                {71.0, 214.65, -2.0, 3.956420}};
			const double h = 6356.766 * z / (6356.766 + z);
			std::size_t i = 0;
			while (i + 1 < std::size(layers) && layers[i + 1][0] <= h)
				++i;
			const auto& layer = layers[i];
			const double temperature = layer[1] + layer[2] * (h - layer[0]);
			const double pressure =
			    layer[3] * (layer[2] == 0.0 ? std::exp(-34.163195 * (h - layer[0]) / layer[1])
			                                : std::pow(layer[1] / temperature, 34.163195 / layer[2]));
			return {pressure / (287.0531 * temperature), temperature, {}};
		}
		constexpr double heights[] = {86, 90, 95, 100, 105, 110, 115, 120, 130, 140};
		constexpr double densities[] = {6.958e-6, 3.416e-6, 1.393e-6, 5.604e-7, 2.325e-7,
		                                9.708e-8, 4.289e-8, 2.222e-8, 8.152e-9, 3.831e-9};
		constexpr double temperatures[] = {186.87, 186.87, 188.42, 195.08, 208.84, 240, 300, 360, 469.27, 559.63};
		std::size_t i = 0;
		while (i + 2 < std::size(heights) && heights[i + 1] < z)
			++i;
		const double f = (z - heights[i]) / (heights[i + 1] - heights[i]);
		return {std::exp(std::log(densities[i]) + f * (std::log(densities[i + 1]) - std::log(densities[i]))),
		        temperatures[i] + f * (temperatures[i + 1] - temperatures[i]),
		        {}};
	}

	Environment::Environment(AtmosphereGrid atmosphere, TerrainGrid terrain)
	    : m_atmosphere(std::move(atmosphere)), m_terrain(std::move(terrain))
	{
		validateAxis(m_atmosphere.latitudesDeg);
		validateAxis(m_atmosphere.longitudesDeg);
		validateAxis(m_terrain.latitudesDeg);
		validateAxis(m_terrain.longitudesDeg);
		if (!std::isfinite(m_atmosphere.firstHeightM) || !std::isfinite(m_atmosphere.heightStepM) ||
		    !std::isfinite(m_atmosphere.firstHeightM +
		                   m_atmosphere.heightStepM * static_cast<double>(m_atmosphere.heightCount)) ||
		    m_atmosphere.heightStepM <= 0.0 || m_atmosphere.heightCount < 2 ||
		    m_atmosphere.samples.size() != checkedSize(m_atmosphere.latitudesDeg.size(),
		                                               m_atmosphere.longitudesDeg.size(), m_atmosphere.heightCount) ||
		    m_terrain.heightsWgs84M.size() !=
		        checkedSize(m_terrain.latitudesDeg.size(), m_terrain.longitudesDeg.size()))
			throw std::invalid_argument("Environment grid dimensions do not match samples");
		for (const auto& sample : m_atmosphere.samples)
			if (!std::isfinite(sample.densityKgM3) || sample.densityKgM3 <= 0.0 ||
			    !std::isfinite(sample.temperatureK) || sample.temperatureK <= 0.0 || !detail::finite(sample.windEnuMps))
				throw std::invalid_argument("Invalid atmosphere sample");
		for (double height : m_terrain.heightsWgs84M)
			if (!std::isfinite(height))
				throw std::invalid_argument("Invalid terrain height");
		m_maxGroundHeightM = *std::max_element(m_terrain.heightsWgs84M.begin(), m_terrain.heightsWgs84M.end());
	}

	Environment Environment::idealized(double groundHeightWgs84M)
	{
		AtmosphereGrid grid;
		grid.heightCount = 2801;
		grid.samples.reserve(grid.heightCount);
		for (std::size_t i = 0; i < grid.heightCount; ++i)
			grid.samples.push_back(standardAtmosphere(static_cast<double>(i) * grid.heightStepM));
		TerrainGrid terrain;
		terrain.heightsWgs84M = {groundHeightWgs84M};
		return Environment(std::move(grid), std::move(terrain));
	}

	AirSample Environment::air(const GeodeticPosition& position) const
	{
		if (!std::isfinite(position.heightWgs84M) || !std::isfinite(position.latitudeRad) ||
		    !std::isfinite(position.longitudeRad))
			throw std::invalid_argument("Non-finite atmospheric position");
		const auto& grid = m_atmosphere;
		const double top = grid.firstHeightM + grid.heightStepM * static_cast<double>(grid.heightCount - 1);
		if (position.heightWgs84M > top + 0.01)
			throw std::out_of_range("Trajectory above the prepared atmosphere");
		const double index = std::clamp((position.heightWgs84M - grid.firstHeightM) / grid.heightStepM, 0.0,
		                                static_cast<double>(grid.heightCount - 1));
		const auto k = std::min(static_cast<std::size_t>(index), grid.heightCount - 2);
		const double f = index - static_cast<double>(k);
		if (grid.latitudesDeg.size() == 1 && grid.longitudesDeg.size() == 1)
			return mix(grid.samples[k], grid.samples[k + 1], f);
		const auto lat = locate(grid.latitudesDeg, position.latitudeRad / detail::deg);
		const auto lon = locate(grid.longitudesDeg, position.longitudeRad / detail::deg);
		const auto column = [&](std::size_t i, std::size_t j) {
			const auto offset = (i * grid.longitudesDeg.size() + j) * grid.heightCount + k;
			return mix(grid.samples[offset], grid.samples[offset + 1], f);
		};
		return mix(mix(column(lat.low, lon.low), column(lat.low, lon.high), lon.fraction),
		           mix(column(lat.high, lon.low), column(lat.high, lon.high), lon.fraction), lat.fraction);
	}

	double Environment::ground(const GeodeticPosition& position) const
	{
		const auto lat = locate(m_terrain.latitudesDeg, position.latitudeRad / detail::deg);
		const auto lon = locate(m_terrain.longitudesDeg, position.longitudeRad / detail::deg);
		const auto height = [&](std::size_t i, std::size_t j) {
			return m_terrain.heightsWgs84M[i * m_terrain.longitudesDeg.size() + j];
		};
		const double a =
		    height(lat.low, lon.low) + lon.fraction * (height(lat.low, lon.high) - height(lat.low, lon.low));
		const double b =
		    height(lat.high, lon.low) + lon.fraction * (height(lat.high, lon.high) - height(lat.high, lon.low));
		return a + lat.fraction * (b - a);
	}

	double Environment::clearance(const GeodeticPosition& position) const
	{
		if (position.heightWgs84M > m_maxGroundHeightM + 1000.0)
			return position.heightWgs84M - m_maxGroundHeightM;
		return position.heightWgs84M - ground(position);
	}

}  // namespace strewn
