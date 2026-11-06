#pragma once

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace orbital_core::hazard::sampling {

	[[nodiscard]] inline uint64_t splitMix64(uint64_t x) noexcept
	{
		x += 0x9e3779b97f4a7c15ULL;
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
		x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
		return x ^ (x >> 31);
	}

	[[nodiscard]] inline double uniform01(uint64_t key) noexcept
	{
		constexpr double kScale = 1.0 / 9007199254740992.0;
		const uint64_t bits = splitMix64(key) >> 11;
		return (static_cast<double>(bits) + 0.5) * kScale;
	}

	[[nodiscard]] inline double standardNormal(uint64_t seed, int sampleIndex, int dimensionIndex) noexcept
	{
		constexpr double kTwoPi = 6.283185307179586476925286766559;
		const uint64_t base = seed ^ (static_cast<uint64_t>(sampleIndex) * 0xd1b54a32d192ed03ULL) ^
		                      (static_cast<uint64_t>(dimensionIndex) * 0xabc98388fb8fac03ULL);
		const double u1 = std::max(uniform01(base), 1.0e-300);
		const double u2 = uniform01(base ^ 0x8cb92ba72f3d8dd7ULL);
		return std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
	}

	template <typename Derived>
	void fillStandardNormal(Eigen::MatrixBase<Derived>& out, uint64_t seed, int sampleIndex)
	{
		for (Eigen::Index i = 0; i < out.size(); ++i)
			out(i) = standardNormal(seed, sampleIndex, static_cast<int>(i));
	}

}  // namespace orbital_core::hazard::sampling
