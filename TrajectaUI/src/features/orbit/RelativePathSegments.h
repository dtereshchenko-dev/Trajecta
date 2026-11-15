#pragma once

#include <algorithm>
#include <array>
#include <cstddef>

namespace relative_path {
	enum class Kind { Unbound, Capture, DistantNegativeEnergy };

	struct Part {
		double start{};
		double end{};
		Kind kind{};
	};

	struct Parts {
		std::array<Part, 3> values{};
		std::size_t count{};
	};

	// Interpolate the same energy and distance-minus-3-Hill values used by capture detection.
	inline Parts split(double energyA, double energyB, double distanceA, double distanceB)
	{
		std::array<double, 4> boundaries{0, 1};
		std::size_t count = 2;
		const auto addCrossing = [&](double a, double b) {
			if ((a < 0 && b > 0) || (a > 0 && b < 0))
				boundaries[count++] = a / (a - b);
		};
		addCrossing(energyA, energyB);
		if (energyA < 0 || energyB < 0)
			addCrossing(distanceA, distanceB);
		std::sort(boundaries.begin(), boundaries.begin() + count);
		Parts result;
		for (std::size_t i = 1; i < count; ++i) {
			const double start = boundaries[i - 1], end = boundaries[i];
			if (end <= start)
				continue;
			const double f = (start + end) / 2;
			const double energy = energyA * (1 - f) + energyB * f;
			const double distance = distanceA * (1 - f) + distanceB * f;
			const auto kind = energy >= 0 ? Kind::Unbound : distance < 0 ? Kind::Capture : Kind::DistantNegativeEnergy;
			if (result.count && result.values[result.count - 1].kind == kind)
				result.values[result.count - 1].end = end;
			else
				result.values[result.count++] = {start, end, kind};
		}
		return result;
	}
}  // namespace relative_path
