#include <gtest/gtest.h>
#include "io/IDataSource.h"
#include "io/MemoryDataSource.h"
#include "io/LoadedFile.h"
#include "io/MmappedFile.h"
#include "io/DataSourceFactory.h"

#include <filesystem>
#include <string>
#include <vector>

namespace {
	std::string writeTempFile(std::string_view payload, std::string_view name)
	{
		std::filesystem::path full = std::filesystem::temp_directory_path() / name;
		std::FILE* fp = nullptr;
#ifdef _MSC_VER
		if (fopen_s(&fp, full.string().c_str(), "wb") != 0 || fp == nullptr)
			return {};
#else
		fp = std::fopen(full.string().c_str(), "wb");
		if (!fp)
			return {};
#endif
		std::fwrite(payload.data(), 1, payload.size(), fp);
		std::fclose(fp);
		return full.string();
	}

	void removeTempFile(const std::string& path) noexcept
	{
		std::error_code ec;
		std::filesystem::remove(path, ec);
	}
}  // namespace

TEST(MemoryDataSourceTest, FromStringRoundTrips)
{
	constexpr std::string_view kPayload = "hello world\n";
	auto src = MemoryDataSource::fromString(kPayload, "fixture-A");

	EXPECT_TRUE(src.isOpen());
	EXPECT_EQ(src.size(), kPayload.size());
	EXPECT_EQ(src.text(), kPayload);
	EXPECT_EQ(src.path(), "fixture-A");
}

TEST(MemoryDataSourceTest, FromBufferRoundTrips)
{
	std::vector<std::byte> buf(8);
	for (std::size_t i = 0; i < buf.size(); ++i)
		buf[i] = static_cast<std::byte>(i + 1);

	MemoryDataSource src(std::move(buf), "binary-fixture");
	EXPECT_TRUE(src.isOpen());
	ASSERT_EQ(src.size(), 8u);
	const auto bytes = src.bytes();
	for (std::size_t i = 0; i < bytes.size(); ++i)
		EXPECT_EQ(static_cast<int>(bytes[i]), static_cast<int>(i + 1));
}

TEST(LoadedFileTest, ReadsSmallFile)
{
	constexpr std::string_view kPayload = "loaded-file-content-12345";
	const auto path = writeTempFile(kPayload, "OrbitalCoreTest_LoadedFile.bin");
	if (path.empty())
		GTEST_SKIP() << "could not write temp file";

	LoadedFile src(path);
	ASSERT_TRUE(src.isOpen());
	EXPECT_EQ(src.size(), kPayload.size());
	EXPECT_EQ(src.text(), kPayload);
	EXPECT_EQ(src.path(), path);

	removeTempFile(path);
}

TEST(LoadedFileTest, MissingFileReportsClosed)
{
	LoadedFile src("/this/path/does/not/exist/anywhere.bin");
	EXPECT_FALSE(src.isOpen());
	EXPECT_EQ(src.size(), 0u);
}

TEST(MmappedFileTest, MapsAndReads)
{
	std::string payload(8192, 'X');
	for (std::size_t i = 0; i < payload.size(); ++i)
		payload[i] = static_cast<char>('a' + (i % 26));

	const auto path = writeTempFile(payload, "OrbitalCoreTest_MmappedFile.bin");
	if (path.empty())
		GTEST_SKIP() << "could not write temp file";

	MmappedFile src(path);
	if (!src.isOpen()) {
		removeTempFile(path);
		GTEST_SKIP() << "mmap unavailable on this filesystem";
	}

	EXPECT_EQ(src.size(), payload.size());
	EXPECT_EQ(src.text(), payload);
	EXPECT_EQ(src.path(), path);

	removeTempFile(path);
}

TEST(MmappedFileTest, MissingFileReportsClosed)
{
	MmappedFile src("/this/path/does/not/exist/anywhere.bin");
	EXPECT_FALSE(src.isOpen());
	EXPECT_EQ(src.size(), 0u);
}

TEST(MmappedFileTest, MovedDonorBecomesInert)
{
	constexpr std::string_view kPayload = "movable";
	const auto path = writeTempFile(kPayload, "OrbitalCoreTest_MmappedFile_Move.bin");
	if (path.empty())
		GTEST_SKIP() << "could not write temp file";

	MmappedFile a(path);
	if (!a.isOpen()) {
		removeTempFile(path);
		GTEST_SKIP() << "mmap unavailable";
	}

	MmappedFile b = std::move(a);
	EXPECT_TRUE(b.isOpen());
	EXPECT_FALSE(a.isOpen());
	EXPECT_EQ(b.text(), kPayload);

	removeTempFile(path);
}

TEST(DataSourceFactoryTest, SmallFileChoosesLoaded)
{
	constexpr std::string_view kPayload = "small-payload-under-threshold";
	const auto path = writeTempFile(kPayload, "OrbitalCoreTest_FactorySmall.bin");
	if (path.empty())
		GTEST_SKIP() << "could not write temp file";

	auto src = DataSourceFactory::open(path);
	ASSERT_NE(src, nullptr);
	EXPECT_TRUE(src->isOpen());
	EXPECT_EQ(src->text(), kPayload);

	removeTempFile(path);
}

TEST(DataSourceFactoryTest, MissingFileNonNullButClosed)
{
	auto src = DataSourceFactory::open("/this/path/does/not/exist/anywhere.bin");
	ASSERT_NE(src, nullptr);
	EXPECT_FALSE(src->isOpen());
}

TEST(DataSourceFactoryTest, ForceLoadedAndMmappedSameContent)
{
	constexpr std::string_view kPayload = "forced-impl-comparison-fixture";
	const auto path = writeTempFile(kPayload, "OrbitalCoreTest_FactoryForce.bin");
	if (path.empty())
		GTEST_SKIP() << "could not write temp file";

	auto loaded = DataSourceFactory::openLoaded(path);
	auto mapped = DataSourceFactory::openMmapped(path);

	ASSERT_TRUE(loaded->isOpen());
	if (!mapped->isOpen()) {
		removeTempFile(path);
		GTEST_SKIP() << "mmap unavailable on this filesystem";
	}

	EXPECT_EQ(loaded->text(), mapped->text());
	removeTempFile(path);
}
