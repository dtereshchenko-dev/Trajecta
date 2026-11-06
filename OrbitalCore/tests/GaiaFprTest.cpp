#include <gtest/gtest.h>
#include "orbital_core/obs/GaiaFpr.h"
#include <cmath>
#include <string>

namespace {
	const std::string metadata =
	    R"([{ "name":"transit_id" },{ "name":"number_mp" },{ "name":"epoch_utc" },{ "name":"ra" },{ "name":"dec" },{ "name":"ra_error_systematic" },{ "name":"dec_error_systematic" },{ "name":"ra_dec_correlation_systematic" },{ "name":"ra_error_random" },{ "name":"dec_error_random" },{ "name":"ra_dec_correlation_random" },{ "name":"x_gaia_geocentric" },{ "name":"y_gaia_geocentric" },{ "name":"z_gaia_geocentric" }])";
}

TEST(GaiaFpr, PreservesTransitIdentityAndCorrelatedUncertainty)
{
	const std::string rows = R"([[90071992547409921,28751,1000,120,20,2,3,0.5,4,5,-0.2,0.01,0.02,0.03],
	                           [90071992547409921,28751,1001,121,21,2,3,0.5,4,5,-0.2,0.01,0.02,0.03],
	                           [90071992547409922,28751,1002,122,22,2,3,0.5,4,5,-0.2,0.01,0.02,0.03]])";
	const auto result =
	    orbital_core::obs::parseGaiaFprResponse("{\"metadata\":" + metadata + ",\"data\":" + rows + "}", 28751);
	ASSERT_EQ(result.size(), 3u);
	EXPECT_DOUBLE_EQ(result[0].epochMjdUtc, 56197.0);
	EXPECT_NEAR(result[0].sigmaRaArcsec, std::sqrt(24.0) / 1000.0, 1e-15);
	EXPECT_NEAR(result[2].sigmaRaArcsec, std::sqrt(20.0) / 1000.0, 1e-15);
	const double cov = (2 * 0.5 * 2 * 3 - 0.2 * 4 * 5) / 1e6;
	EXPECT_NEAR(result[0].correlation, cov / (result[0].sigmaRaArcsec * result[0].sigmaDecArcsec), 1e-15);
	EXPECT_NEAR(result[0].observerOffsetEquatorialAu.x, 0.01 * (1.0 - 1.550519768e-8), 1e-16);
}

TEST(GaiaFpr, RejectsWrongObjectOrTruncatedData)
{
	const std::string row = R"([[42,1,1000,120,20,2,3,0.5,4,5,-0.2,0.01,0.02,0.03]])";
	EXPECT_THROW(
	    (void)orbital_core::obs::parseGaiaFprResponse("{\"metadata\":" + metadata + ",\"data\":" + row + "}", 28751),
	    std::runtime_error);
	EXPECT_THROW((void)orbital_core::obs::parseGaiaFprResponse("{\"metadata\":[],\"data\":[]}", 28751),
	             std::runtime_error);
	const std::string badCoordinates = R"([[42,28751,1000,360,20,2,3,0.5,4,5,-0.2,0.01,0.02,0.03]])";
	EXPECT_THROW((void)orbital_core::obs::parseGaiaFprResponse(
	                 "{\"metadata\":" + metadata + ",\"data\":" + badCoordinates + "}", 28751),
	             std::runtime_error);
}
