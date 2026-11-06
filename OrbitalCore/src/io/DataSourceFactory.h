#pragma once

#include "IDataSource.h"
#include <memory>
#include <string>

class ORBITAL_CORE_API DataSourceFactory {
   public:
	[[nodiscard]] static std::unique_ptr<IDataSource> open(const std::string& filePath);

	// Force `LoadedFile` regardless of size.  Useful for tests where mmap
	// would be brittle (Windows file-locking) and for very-small fixtures
	// where mmap setup time dominates a microbenchmark.
	[[nodiscard]] static std::unique_ptr<IDataSource> openLoaded(const std::string& filePath);

	// Force `MmappedFile`.  Returns a source with `isOpen() == false` if
	// mmap fails — caller can fall back to `openLoaded`.  Mostly useful in
	// tests that explicitly want to verify the mmap codepath.
	[[nodiscard]] static std::unique_ptr<IDataSource> openMmapped(const std::string& filePath);

   private:
	DataSourceFactory() = delete;
	~DataSourceFactory() = delete;
};
