#include <gtest/gtest.h>
#include "orbital_core/SpkReader.h"
#include "orbital_core/RuntimeData.h"
#include "orbital_core/Types.h"
#include "orbital_core/Vec3D.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <thread>
#include <vector>

class SpkReaderTest : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (getExecutableDir().empty()) {
			auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
			std::filesystem::path candidates[] = {
			    srcDir.parent_path() / "data",
			    srcDir / "data",
			    srcDir / "release",
			    srcDir / "debug",
			};
			for (const auto& dir : candidates) {
				if (std::filesystem::exists(dir / "sb441-n16.bsp")) {
					getExecutableDir() = dir;
					break;
				}
			}
		}
		reader = std::make_unique<SpkReader>();
	}

	std::unique_ptr<SpkReader> reader;
};

TEST_F(SpkReaderTest, TargetIdsMatchExpected)
{
	if (!reader->isAvailable())
		GTEST_SKIP();

	const std::vector<int> expected = {
	    2000001, 2000002, 2000003, 2000004, 2000007, 2000010, 2000015, 2000016,
	    2000031, 2000052, 2000065, 2000087, 2000088, 2000107, 2000511, 2000704,
	};

	auto ids = reader->targetIds();
	std::sort(ids.begin(), ids.end());
	EXPECT_EQ(ids, expected) << "unexpected NAIF id set in sb441-n16.bsp";
}

TEST_F(SpkReaderTest, CeresPositionAtJ2000)
{
	if (!reader->isAvailable())
		GTEST_SKIP();

	Vec3D R, V;
	ASSERT_TRUE(reader->queryState(2000001, Constants::J2000_JD, R, V));

	EXPECT_NEAR(R.x, -2.379327705915647e+00, 1e-7);
	EXPECT_NEAR(R.y, 7.954860388931395e-01, 1e-7);
	EXPECT_NEAR(R.z, 4.630055715902157e-01, 1e-7);

	EXPECT_NEAR(V.x, -3.584228273182221e-03, 1e-9);
	EXPECT_NEAR(V.y, -1.058594861536097e-02, 1e-9);
	EXPECT_NEAR(V.z, 3.338536192580540e-04, 1e-9);
}

TEST_F(SpkReaderTest, VelocityMatchesFiniteDifference)
{
	if (!reader->isAvailable())
		GTEST_SKIP();

	const double jd = Constants::J2000_JD;
	const double h = 60.0 / 86400.0;

	Vec3D R0, V0;
	ASSERT_TRUE(reader->queryState(2000001, jd, R0, V0));
	Vec3D Rm, dummy;
	ASSERT_TRUE(reader->queryState(2000001, jd - h, Rm, dummy));
	Vec3D Rp;
	ASSERT_TRUE(reader->queryState(2000001, jd + h, Rp, dummy));

	Vec3D V_fd;
	V_fd.x = (Rp.x - Rm.x) / (2.0 * h);
	V_fd.y = (Rp.y - Rm.y) / (2.0 * h);
	V_fd.z = (Rp.z - Rm.z) / (2.0 * h);

	EXPECT_NEAR(V0.x, V_fd.x, 1e-9);
	EXPECT_NEAR(V0.y, V_fd.y, 1e-9);
	EXPECT_NEAR(V0.z, V_fd.z, 1e-9);
}

TEST_F(SpkReaderTest, UnknownTargetReturnsFalse)
{
	if (!reader->isAvailable())
		GTEST_SKIP();
	Vec3D R, V;
	EXPECT_FALSE(reader->queryState(2000005, Constants::J2000_JD, R, V));
}

TEST_F(SpkReaderTest, OutOfRangeJdReturnsFalse)
{
	if (!reader->isAvailable())
		GTEST_SKIP();

	Vec3D R, V;

	EXPECT_FALSE(reader->queryState(2000001, 1.0e8, R, V)) << "query past end";

	EXPECT_FALSE(reader->queryState(2000001, -1.0e7, R, V)) << "query before start";
}

TEST_F(SpkReaderTest, FallbackWhenFileAbsent)
{
	SpkReader noFile("C:/__nonexistent__/sb441-n16.bsp");
	EXPECT_FALSE(noFile.isAvailable());
	EXPECT_TRUE(noFile.targetIds().empty());

	Vec3D R, V;
	EXPECT_FALSE(noFile.queryState(2000001, Constants::J2000_JD, R, V));
}

TEST_F(SpkReaderTest, LoadsFromEphemerisDataDirectory)
{
	if (!reader->isAvailable())
		GTEST_SKIP() << "SB441 test data is unavailable";

	const auto sourcePath = reader->sourcePath();
	const auto kernel = std::filesystem::path(std::u8string(sourcePath.begin(), sourcePath.end()));
	const auto previousEphemeris = orbital_core::ephemerisDataDirectory();
	const auto previousExecutable = getExecutableDir();

	struct RestoreDirectories {
		std::filesystem::path ephemeris;
		std::filesystem::path executable;

		~RestoreDirectories()
		{
			orbital_core::setEphemerisDataDirectory(ephemeris);
			getExecutableDir() = executable;
		}
	} restore{previousEphemeris, previousExecutable};

	orbital_core::setEphemerisDataDirectory(kernel.parent_path());
	getExecutableDir() = kernel.parent_path() / "missing-directory";
	SpkReader fromDataDirectory;
	ASSERT_TRUE(fromDataDirectory.isAvailable());
	const auto actualPath = fromDataDirectory.sourcePath();
	EXPECT_EQ(std::filesystem::path(std::u8string(actualPath.begin(), actualPath.end())), kernel);
}

TEST_F(SpkReaderTest, ConcurrentQueriesProduceIdenticalResults)
{
	if (!reader->isAvailable())
		GTEST_SKIP();

	constexpr int nThreads = 8;
	constexpr int nQueries = 100;
	std::vector<std::thread> workers;
	std::vector<std::vector<Vec3D>> results(nThreads);
	std::atomic<bool> stop{false};

	for (int t = 0; t < nThreads; ++t) {
		results[t].reserve(nQueries);
		workers.emplace_back([&, t]() {
			for (int q = 0; q < nQueries; ++q) {
				if (stop.load())
					break;
				Vec3D R, V;
				reader->queryState(2000001, Constants::J2000_JD + q * 0.001, R, V);
				results[t].push_back(R);
			}
		});
	}
	for (auto& w : workers)
		w.join();

	for (int q = 0; q < nQueries; ++q) {
		for (int t = 1; t < nThreads; ++t) {
			ASSERT_EQ(results[t][q].x, results[0][q].x) << "thread " << t << " q=" << q;
			ASSERT_EQ(results[t][q].y, results[0][q].y);
			ASSERT_EQ(results[t][q].z, results[0][q].z);
		}
	}
}
