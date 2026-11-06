#include "DataSourceFactory.h"
#include <memory>
#include <string>
#include "MmappedFile.h"
#include "LoadedFile.h"

#include <filesystem>

std::unique_ptr<IDataSource> DataSourceFactory::open(const std::string& filePath)
{
	std::error_code ec;
	const auto sz = std::filesystem::file_size(filePath, ec);
	if (ec)
		return std::make_unique<LoadedFile>(filePath);

	static constexpr std::size_t kMmapThresholdBytes = 16ull * 1024 * 1024;
	if (sz >= kMmapThresholdBytes) {
		auto mapped = std::make_unique<MmappedFile>(filePath);
		if (mapped->isOpen())
			return mapped;
	}

	return std::make_unique<LoadedFile>(filePath);
}

std::unique_ptr<IDataSource> DataSourceFactory::openLoaded(const std::string& filePath)
{
	return std::make_unique<LoadedFile>(filePath);
}

std::unique_ptr<IDataSource> DataSourceFactory::openMmapped(const std::string& filePath)
{
	return std::make_unique<MmappedFile>(filePath);
}
