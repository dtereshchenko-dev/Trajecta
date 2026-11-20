#include <gtest/gtest.h>
#include "orbital_core/SpkReader.h"
#include "orbital_core/RuntimeData.h"
#include "orbital_core/Types.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/perturbations/AsteroidPerturbations.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
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
		if (const char* path = std::getenv("TRAJECTA_SB441"))
			reader = std::make_unique<SpkReader>(path);
		else
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
	R = {1.0, 2.0, 3.0};
	EXPECT_FALSE(reader->queryPositionOffsetMjd(2000001, 1.0e8, 0.0, R));
	EXPECT_DOUBLE_EQ(R.x, 1.0);
	EXPECT_DOUBLE_EQ(R.y, 2.0);
	EXPECT_DOUBLE_EQ(R.z, 3.0);
}

TEST_F(SpkReaderTest, AsteroidForceFailsOutsideKernelCoverage)
{
	if (!reader->isAvailable())
		GTEST_SKIP() << "SB441 test data is unavailable";
	orbital_core::forces::ForceContext ctx{};
	ctx.R = {1.0, 0.0, 0.0};
	ctx.asteroidSource = reader.get();
	ctx.useSb441 = true;
	Vec3D acceleration{};
	Eigen::Matrix3d dr = Eigen::Matrix3d::Zero(), dv = Eigen::Matrix3d::Zero();
	orbital_core::forces::AsteroidPerturbations force;
	ctx.setEpochMjd(51544.5);
	EXPECT_NO_THROW(force.evaluate(ctx, acceleration, dr, dv));
	EXPECT_GT(acceleration.magnitude(), 0.0);
	ctx.setEpochMjd(1.0e8);
	EXPECT_THROW(force.evaluate(ctx, acceleration, dr, dv), EphemerisError);
}

TEST(SpkReaderForceTest, ExplicitSb441DoesNotSilentlyFallBackWhenUnavailable)
{
	orbital_core::forces::ForceContext ctx{};
	ctx.R = {1.0, 0.0, 0.0};
	ctx.useSb441 = true;
	ctx.setEpochMjd(51544.5);
	Vec3D acceleration{};
	Eigen::Matrix3d dr = Eigen::Matrix3d::Zero(), dv = Eigen::Matrix3d::Zero();
	orbital_core::forces::AsteroidPerturbations force;
	EXPECT_THROW(force.evaluate(ctx, acceleration, dr, dv), EphemerisError);
	ctx.useSb441 = false;
	EXPECT_NO_THROW(force.evaluate(ctx, acceleration, dr, dv));
}

TEST_F(SpkReaderTest, FallbackWhenFileAbsent)
{
	SpkReader noFile("C:/__nonexistent__/sb441-n16.bsp");
	EXPECT_FALSE(noFile.isAvailable());
	EXPECT_TRUE(noFile.targetIds().empty());

	Vec3D R, V;
	EXPECT_FALSE(noFile.queryState(2000001, Constants::J2000_JD, R, V));
}

TEST(SpkReaderFileTest, InvalidSegmentDoesNotExposePartialTargetIndex)
{
	std::array<unsigned char, 3072> bytes{};
	const auto writeInt = [&](std::size_t offset, std::int32_t value) {
		std::memcpy(bytes.data() + offset, &value, sizeof(value));
	};
	const auto writeDouble = [&](std::size_t offset, double value) {
		std::memcpy(bytes.data() + offset, &value, sizeof(value));
	};
	std::memcpy(bytes.data(), "DAF/SPK ", 8);
	std::memcpy(bytes.data() + 88, "LTL-IEEE", 8);
	writeInt(8, 2);
	writeInt(12, 6);
	writeInt(76, 2);
	writeDouble(1024 + 16, 1.0);
	constexpr std::size_t summary = 1024 + 24;
	writeDouble(summary, 0.0);
	writeDouble(summary + 8, 86400.0);
	writeInt(summary + 16, 2000001);
	writeInt(summary + 20, 10);
	writeInt(summary + 24, 1);
	writeInt(summary + 28, 2);
	writeInt(summary + 32, 257);
	writeInt(summary + 36, 265);
	// The segment tail has a zero interval length and must be rejected.
	writeDouble(2088 + 16, 5.0);
	writeDouble(2088 + 24, 1.0);
	const auto path = std::filesystem::temp_directory_path() /
	                  ("trajecta-invalid-spk-" +
	                   std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".bsp");

	struct RemoveFile {
		std::filesystem::path path;

		~RemoveFile()
		{
			std::error_code error;
			std::filesystem::remove(path, error);
		}
	} remove{path};

	{
		std::ofstream output(path, std::ios::binary);
		ASSERT_TRUE(output);
		output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
		ASSERT_TRUE(output);
	}
	const auto utf8 = path.u8string();
	SpkReader reader(std::string(utf8.begin(), utf8.end()));
	EXPECT_FALSE(reader.isAvailable());
	EXPECT_TRUE(reader.targetIds().empty());
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
