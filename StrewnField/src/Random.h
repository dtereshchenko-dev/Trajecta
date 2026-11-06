#pragma once

#include "strewn_core/Types.h"
#include <random>

namespace strewn::detail {

	inline std::uint64_t mixSeed(std::uint64_t value)
	{
		value += 0x9e3779b97f4a7c15ULL;
		value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
		value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
		return value ^ (value >> 31);
	}

	class Random {
	   public:
		Random(std::uint64_t seed, std::uint64_t run, std::uint64_t stream)
		    : m_generator(mixSeed(mixSeed(seed ^ mixSeed(run)) ^ mixSeed(stream + 0xd1b54a32d192ed03ULL)))
		{
		}

		double uniform() { return static_cast<double>(m_generator() >> 11) * 0x1.0p-53; }

		double uniform(Range range) { return range.min + uniform() * (range.max - range.min); }

	   private:
		std::mt19937_64 m_generator;
	};

}  // namespace strewn::detail
