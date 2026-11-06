#include "orbital_core/obs/SbdbPhysical.h"
#include <gtest/gtest.h>

TEST(SbdbPhysical, UsesKilometreDiameterForRequestedObject)
{
	const auto radius = orbital_core::obs::parseSbdbRadiusKm(
	    R"({"object":{"des":"11468"},"phys_par":[{"name":"H","value":"14.4"},{"name":"diameter","value":"6.472","units":"km"}]})",
	    11468);
	ASSERT_TRUE(radius);
	EXPECT_DOUBLE_EQ(*radius, 3.236);
}

TEST(SbdbPhysical, MissingOrInvalidDiameterIsUnavailable)
{
	using orbital_core::obs::parseSbdbRadiusKm;
	EXPECT_FALSE(parseSbdbRadiusKm(R"({"object":{"des":"42"},"phys_par":[]})", 42));
	EXPECT_FALSE(
	    parseSbdbRadiusKm(R"({"object":{"des":"42"},"phys_par":[{"name":"diameter","value":"0","units":"km"}]})", 42));
	EXPECT_FALSE(
	    parseSbdbRadiusKm(R"({"object":{"des":"42"},"phys_par":[{"name":"diameter","value":"3","units":"m"}]})", 42));
	EXPECT_FALSE(
	    parseSbdbRadiusKm(R"({"object":{"des":"42"},"phys_par":[{"name":"diameter","value":"3x","units":"km"}]})", 42));
}

TEST(SbdbPhysical, RejectsResponseForAnotherObject)
{
	EXPECT_THROW((void)orbital_core::obs::parseSbdbRadiusKm(R"({"object":{"des":"43"},"phys_par":[]})", 42),
	             std::runtime_error);
}
