#pragma once

#include "orbital_core/hazard/SemilinearContours.h"
#include <nlohmann/json.hpp>

namespace orbital_core::hazard {
	inline nlohmann::json semilinearJson(const SemilinearContours& value)
	{
		nlohmann::json result = {{"status", value.status},
		                         {"reason", value.reason},
		                         {"method", "nominal covariance regression; nonlinear boundary propagation"},
		                         {"contours", nlohmann::json::array()}};
		for (const auto& contour : value.contours) {
			nlohmann::json points = nlohmann::json::array();
			for (const auto& point : contour.points) {
				if (point)
					points.push_back({{"epoch", point->mjdTdb},
					                  {"latitude_deg", point->latitudeDeg},
					                  {"longitude_deg", point->longitudeDeg}});
				else
					points.push_back(nullptr);
			}
			result["contours"].push_back({{"sigma", contour.sigma},
			                              {"missed", contour.missed},
			                              {"failed", contour.failed},
			                              {"points", std::move(points)}});
		}
		return result;
	}
}  // namespace orbital_core::hazard
