#pragma once

#include <cmath>

namespace orbital_core::obs {

	// Terrestrial WGS84 position supplied with this observation, not a catalog
	// station or an inertial satellite offset. Longitude/latitude are degrees,
	// height is metres above the ellipsoid. Preserve the submitted representation.
	struct RovingObserver {
		double longitudeDeg{};
		double latitudeDeg{};
		double heightMeters{};

		[[nodiscard]] bool valid() const noexcept
		{
			return std::isfinite(longitudeDeg) && longitudeDeg >= -180.0 && longitudeDeg <= 360.0 &&
			       std::isfinite(latitudeDeg) && std::abs(latitudeDeg) <= 90.0 && std::isfinite(heightMeters);
		}

		bool operator==(const RovingObserver&) const = default;
	};

}  // namespace orbital_core::obs
