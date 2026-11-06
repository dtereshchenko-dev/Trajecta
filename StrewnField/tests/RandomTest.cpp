#include "Random.h"
#include <gtest/gtest.h>
#include <set>

namespace {

	TEST(Random, RunAndStreamRolesCannotCancelOrSwap)
	{
		std::set<double> firstDraws;
		for (std::uint64_t run = 0; run < 32; ++run)
			for (std::uint64_t stream = 0; stream < 3; ++stream) {
				strewn::detail::Random random(20230213, run, stream);
				EXPECT_TRUE(firstDraws.insert(random.uniform()).second) << run << ':' << stream;
			}
	}

	TEST(Random, SequenceIsRepeatableAndHalfOpen)
	{
		strewn::detail::Random a(42, 7, 2), b(42, 7, 2);
		for (unsigned i = 0; i < 1000; ++i) {
			const double value = a.uniform();
			EXPECT_DOUBLE_EQ(value, b.uniform());
			EXPECT_GE(value, 0.0);
			EXPECT_LT(value, 1.0);
		}
	}

}  // namespace
