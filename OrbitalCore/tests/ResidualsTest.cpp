#include <gtest/gtest.h>
#include "orbital_core/Residuals.h"
#include <cmath>
#include <sstream>
#include <filesystem>

TEST(ResidualsTest, EmptyRMS)
{
	Residuals r;
	EXPECT_EQ(r.count(), 0u);
	EXPECT_DOUBLE_EQ(r.rms(), 0.0);
}

TEST(ResidualsTest, AddAndCount)
{
	Residuals r;
	r.add({2460000.5, 1.0, 0.5, 0.1, 0.2, 0.224, "585"});
	r.add({2460001.5, 1.1, 0.6, -0.3, 0.1, 0.316, "585"});

	EXPECT_EQ(r.count(), 2u);
}

TEST(ResidualsTest, RMSCalculation)
{
	Residuals r;
	r.add({2460000.5, 0, 0, 0.6, 0.8, 1.0, "585"});
	r.add({2460001.5, 0, 0, 0.8, 0.6, 1.0, "585"});

	EXPECT_NEAR(r.rms(), std::sqrt(0.5), 0.001);
}

TEST(ResidualsTest, RMSNonTrivial)
{
	Residuals r;
	r.add({0, 0, 0, 0.3, 0.4, 0.5, "585"});
	r.add({0, 0, 0, 1.2, 1.6, 2.0, "585"});

	EXPECT_NEAR(r.rms(), std::sqrt(1.0625), 0.001);
}

TEST(ResidualsTest, Clear)
{
	Residuals r;
	r.add({0, 0, 0, 1.0, 1.0, 1.414, "585"});
	EXPECT_EQ(r.count(), 1u);

	r.clear();
	EXPECT_EQ(r.count(), 0u);
	EXPECT_DOUBLE_EQ(r.rms(), 0.0);
}

TEST(ResidualsTest, DataAccess)
{
	Residuals r;
	r.add({2460000.5, 1.5, 0.3, 0.5, -0.2, 0.539, "T08"});

	const auto& data = r.data();
	ASSERT_EQ(data.size(), 1u);
	EXPECT_DOUBLE_EQ(data[0].jd, 2460000.5);
	EXPECT_EQ(data[0].obsCode, "T08");
	EXPECT_NEAR(data[0].dRA_arcsec, 0.5, 1e-10);
	EXPECT_NEAR(data[0].dDEC_arcsec, -0.2, 1e-10);
}

TEST(ResidualsTest, PrintDoesNotCrash)
{
	Residuals r;
	r.add({2460000.5, 1.0, 0.5, 0.1, 0.2, 0.224, "585"});
	r.add({2460001.5, 1.1, 0.6, -0.3, 0.1, 0.316, "703"});

	std::ostringstream oss;
	r.print(oss);
	EXPECT_FALSE(oss.str().empty());
	EXPECT_NE(oss.str().find("585"), std::string::npos);
	EXPECT_NE(oss.str().find("703"), std::string::npos);
}

TEST(ResidualsTest, SaveToFile)
{
	Residuals r;
	r.add({2460000.5, 1.0, 0.5, 0.1, 0.2, 0.224, "585"});

	auto tmpPath = std::filesystem::temp_directory_path() / "test_residuals.txt";
	EXPECT_TRUE(r.save(tmpPath));
	EXPECT_TRUE(std::filesystem::exists(tmpPath));

	std::filesystem::remove(tmpPath);
}
