#include "orbital_core/RuntimeData.h"
#include <gtest/gtest.h>

int main(int argc, char** argv)
{
	orbital_core::setObservatoryDataDirectory(ORBITAL_TEST_OBSERVATORY_DIR);
	orbital_core::setEphemerisDataDirectory(ORBITAL_TEST_OBSERVATORY_DIR);
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
