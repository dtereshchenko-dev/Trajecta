#include <gtest/gtest.h>
#include "orbital_core/Constants.h"
#include "Gauss.h"
#include <cmath>

using Calculations::detail::findInitialGuesses;

TEST(FindInitialGuesses, SimplePolynomialWithKnownRoot)
{
	auto guesses = findInitialGuesses(0.0, 0.0, -1.0, -0.5, 2.0, 0.001);
	ASSERT_FALSE(guesses.empty());

	bool foundNearOne = false;
	for (double g : guesses) {
		if (std::abs(g - 1.0) < 0.01) {
			foundNearOne = true;
			break;
		}
	}

	EXPECT_TRUE(foundNearOne) << "Expected a root near 1.0";
}

TEST(FindInitialGuesses, NoRootsInRange)
{
	auto guesses = findInitialGuesses(0.0, 0.0, 1.0, 0.0, 10.0, 0.01);
	EXPECT_TRUE(guesses.empty());
}

TEST(FindInitialGuesses, InvalidStepReturnsEmpty)
{
	auto guesses = findInitialGuesses(0.0, 0.0, -1.0, 0.0, 10.0, 0.0);
	EXPECT_TRUE(guesses.empty());

	guesses = findInitialGuesses(0.0, 0.0, -1.0, 0.0, 10.0, -1.0);
	EXPECT_TRUE(guesses.empty());
}

TEST(FindInitialGuesses, ReversedRangeReturnsEmpty)
{
	auto guesses = findInitialGuesses(0.0, 0.0, -1.0, 10.0, 0.0, 0.01);
	EXPECT_TRUE(guesses.empty());
}

TEST(FindInitialGuesses, MultipleRoots)
{
	auto guesses = findInitialGuesses(-10.0, 5.0, -1.0, 0.0, 5.0, 0.001);
	EXPECT_GE(guesses.size(), 1u);
	for (double g : guesses) {
		EXPECT_GT(g, 0.0);
	}
}

TEST(FindInitialGuesses, PositiveRootsOnly)
{
	auto guesses = findInitialGuesses(-10.0, 20.0, -12.0, 0.0, 10.0, 0.001);
	ASSERT_FALSE(guesses.empty());
	for (double g : guesses) {
		EXPECT_GE(g, 0.0);
		EXPECT_LE(g, 10.0 + 0.001);
	}
}

TEST(Gauss, InvalidObservationsReturnEmpty)
{
	ObsNode obs1, obs2, obs3;
	obs1.JD = 2458634.0;
	obs2.JD = 2458635.0;
	obs3.JD = 2458636.0;

	obs1.vector = obs2.vector = obs3.vector = Vec3D(0.0, 0.0, 0.0);
	obs1.obsVector = obs2.obsVector = obs3.obsVector = Vec3D(0.0, 0.0, 0.0);

	auto result = Calculations::gauss(obs1, obs2, obs3);
	EXPECT_FALSE(result.isValid());
}

TEST(Gauss, CollinearObservationsReturnEmpty)
{
	ObsNode obs1, obs2, obs3;
	obs1.JD = 2458634.0;
	obs2.JD = 2458635.0;
	obs3.JD = 2458636.0;

	Vec3D dir(1.0, 0.0, 0.0);
	obs1.vector = obs2.vector = obs3.vector = dir;
	obs1.obsVector = obs2.obsVector = obs3.obsVector = Vec3D(1.0, 0.0, 0.0);

	auto result = Calculations::gauss(obs1, obs2, obs3);
	EXPECT_FALSE(result.isValid());
}
