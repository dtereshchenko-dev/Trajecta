#include "EntrySummary.h"
#include "VectorMath.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace strewn::output {
	namespace {
		using Json = nlohmann::json;

		struct Value {
			double value, weight;
		};

		Json describe(const std::vector<Value>& values, bool circular = false)
		{
			Json result = {{"mean", nullptr}, {"sigma", nullptr}, {"samples", values.size()}};
			if (values.empty())
				return result;
			double weight = 0, mean = 0;
			const double origin = values.front().value;
			if (circular) {
				double x = 0, y = 0;
				for (const auto& v : values) {
					weight += v.weight;
					x += v.weight * std::cos(v.value * detail::deg);
					y += v.weight * std::sin(v.value * detail::deg);
				}
				if (std::hypot(x, y) < 1e-10 * weight)
					return result;
				mean = std::atan2(y, x) / detail::deg;
			} else {
				for (const auto& v : values) {
					weight += v.weight;
					mean += v.weight * (v.value - origin);
				}
				mean = origin + mean / weight;
			}
			result["mean"] = mean;
			if (values.size() > 1) {
				double variance = 0;
				for (const auto& v : values) {
					const double delta = circular ? std::remainder(v.value - mean, 360.) : v.value - mean;
					variance += v.weight * delta * delta;
				}
				result["sigma"] = std::sqrt(variance / weight);
			}
			return result;
		}

		double dot(const Vector3& a, const Vector3& b) { return std::inner_product(a.begin(), a.end(), b.begin(), 0.); }
	}  // namespace

	Json summarizeEntries(std::span<const EntryState> entries)
	{
		double maxWeight = 0;
		for (const auto& e : entries) {
			if (!std::isfinite(e.weight) || e.weight < 0)
				throw std::invalid_argument("Invalid entry summary weight");
			maxWeight = std::max(maxWeight, e.weight);
		}
		if (!(maxWeight > 0))
			throw std::invalid_argument("Entry summary has no positive weights");
		std::vector<Value> times, latitudes, longitudes, heights, speeds, slopes, headings;
		double weight = 0, weightSquared = 0, epoch = 0;
		for (const auto& e : entries) {
			if (e.weight == 0)
				continue;
			if (!detail::finite(e.positionEcefM) || !detail::finite(e.velocityEcefMps) ||
			    !std::isfinite(e.epochMjdUtc) || detail::norm(e.positionEcefM) < 1)
				throw std::invalid_argument("Invalid entry summary state");
			const double w = e.weight / maxWeight;
			if (times.empty())
				epoch = e.epochMjdUtc;
			const auto p = geodetic(e.positionEcefM);
			weight += w;
			weightSquared += w * w;
			times.push_back({(e.epochMjdUtc - epoch) * 86400, w});
			latitudes.push_back({p.latitudeRad / detail::deg, w});
			longitudes.push_back({p.longitudeRad / detail::deg, w});
			heights.push_back({p.heightWgs84M / 1000, w});
			const double speed = detail::norm(e.velocityEcefMps);
			speeds.push_back({speed / 1000, w});
			const double east = dot(e.velocityEcefMps, enuToEcef({1, 0, 0}, p.latitudeRad, p.longitudeRad));
			const double north = dot(e.velocityEcefMps, enuToEcef({0, 1, 0}, p.latitudeRad, p.longitudeRad));
			const double up = dot(e.velocityEcefMps, enuToEcef({0, 0, 1}, p.latitudeRad, p.longitudeRad));
			if (speed > 1e-9)
				slopes.push_back({std::atan2(up, std::hypot(east, north)) / detail::deg, w});
			if (std::hypot(east, north) > 1e-10 * std::max(1., speed) && std::abs(std::cos(p.latitudeRad)) > 1e-10)
				headings.push_back({std::atan2(east, north) / detail::deg, w});
		}
		const auto time = describe(times);
		auto heading = describe(headings, true);
		if (!heading["mean"].is_null())
			heading["mean"] = std::fmod(heading["mean"].get<double>() + 360., 360.);
		Json result = {{"source", "entries.json"},
		               {"total_entries", entries.size()},
		               {"positive_weight_entries", times.size()},
		               {"effective_entries", weight * weight / weightSquared},
		               {"epoch_mjd_utc", epoch + time["mean"].get<double>() / 86400},
		               {"time_sigma_s", time["sigma"]},
		               {"latitude_deg", describe(latitudes)},
		               {"longitude_deg", describe(longitudes, true)},
		               {"height_wgs84_km", describe(heights)},
		               {"ecef_speed_km_s", describe(speeds)},
		               {"flight_path_deg", describe(slopes)},
		               {"heading_deg", heading},
		               {"ellipse", {{"semi_major_m", nullptr}, {"semi_minor_m", nullptr}, {"azimuth_deg", nullptr}}},
		               {"warnings", Json::array()}};
		if (times.size() == 1)
			result["warnings"].push_back(
			    "Only one positive-weight entry: uncertainty cannot be estimated from this archive.");
		if (headings.size() != times.size())
			result["warnings"].push_back(
			    "Heading excludes stationary, vertical and polar states where azimuth is undefined.");
		if (slopes.size() != times.size())
			result["warnings"].push_back("Flight path angle excludes stationary states.");
		if (result["longitude_deg"]["mean"].is_null()) {
			result["warnings"].push_back(
			    "Longitude mean and local ellipse are undefined for this dispersed entry ensemble.");
			return result;
		}
		const double lat = result["latitude_deg"]["mean"].get<double>() * detail::deg;
		const double lon = result["longitude_deg"]["mean"].get<double>() * detail::deg;
		const auto east = enuToEcef({1, 0, 0}, lat, lon), north = enuToEcef({0, 1, 0}, lat, lon);
		Vector3 mean =
		    std::find_if(entries.begin(), entries.end(), [](const auto& e) { return e.weight > 0; })->positionEcefM;
		for (std::size_t axis = 0; axis < 3; ++axis) {
			double delta = 0;
			for (const auto& e : entries)
				if (e.weight > 0)
					delta += (e.weight / maxWeight) * (e.positionEcefM[axis] - mean[axis]);
			mean[axis] += delta / weight;
		}
		Eigen::Matrix2d covariance = Eigen::Matrix2d::Zero();
		double maxAngle = 0;
		for (const auto& e : entries)
			if (e.weight > 0) {
				Vector3 delta;
				for (std::size_t axis = 0; axis < 3; ++axis)
					delta[axis] = e.positionEcefM[axis] - mean[axis];
				const Eigen::Vector2d offset(dot(delta, east), dot(delta, north));
				covariance += (e.weight / maxWeight / weight) * (offset * offset.transpose());
				const auto p = geodetic(e.positionEcefM);
				const double angle =
				    std::acos(std::clamp(std::sin(lat) * std::sin(p.latitudeRad) +
				                             std::cos(lat) * std::cos(p.latitudeRad) * std::cos(p.longitudeRad - lon),
				                         -1., 1.));
				maxAngle = std::max(maxAngle, angle / detail::deg);
			}
		result["max_angular_radius_deg"] = maxAngle;
		if (maxAngle > 5)
			result["warnings"].push_back(
			    "Wide entry cloud: a local tangent-plane ellipse is only an approximation, not a global corridor.");
		if (times.size() > 1) {
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(covariance);
			if (solver.info() != Eigen::Success)
				throw std::runtime_error("Cannot calculate entry ellipse");
			const double small = std::max(0., solver.eigenvalues()[0]), large = std::max(0., solver.eigenvalues()[1]);
			result["ellipse"]["semi_major_m"] = std::sqrt(large);
			result["ellipse"]["semi_minor_m"] = std::sqrt(small);
			result["ellipse"]["covariance_en_m2"] = {{covariance(0, 0), covariance(0, 1)},
			                                         {covariance(1, 0), covariance(1, 1)}};
			if (large > 1e-16 && large - small > 1e-8 * large && std::abs(std::cos(lat)) > 1e-10)
				result["ellipse"]["azimuth_deg"] = std::fmod(
				    std::atan2(solver.eigenvectors()(0, 1), solver.eigenvectors()(1, 1)) / detail::deg + 180., 180.);
		}
		return result;
	}
}  // namespace strewn::output
