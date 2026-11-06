#pragma once

#include <cmath>

namespace orbital_core::detail {

	// TwoSum retains the bits lost when an integration offset is added to its epoch.
	// Only collapse the pair after subtracting the nearby ephemeris record origin.
	struct SplitTime {
		double high;
		double low;

		static SplitTime sum(double a, double b) noexcept
		{
			const double high = a + b;
			const double bVirtual = high - a;
			return {high, (a - (high - bVirtual)) + (b - bVirtual)};
		}

		SplitTime shifted(double origin) const noexcept
		{
			const auto delta = sum(high, -origin);
			return sum(delta.high, delta.low + low);
		}

		double since(double origin) const noexcept
		{
			const auto delta = sum(high, -origin);
			return delta.high + (delta.low + low);
		}

		bool atOrAfter(double origin) const noexcept { return high > origin || (high == origin && low >= 0.0); }

		bool atOrBefore(double origin) const noexcept { return high < origin || (high == origin && low <= 0.0); }

		SplitTime scaled(double factor) const noexcept
		{
			const double product = high * factor;
			return sum(product, std::fma(high, factor, -product) + low * factor);
		}
	};

}  // namespace orbital_core::detail
