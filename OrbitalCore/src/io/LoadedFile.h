#pragma once

#include "IDataSource.h"
#include <string>
#include <string_view>
#include <vector>

// ---------------------------------------------------------------------------
// LoadedFile
//
// Reads a file in full into memory at construction time.  Used for small
// resources (under `DataSourceFactory::kMmapThresholdBytes`) where mmap
// setup overhead exceeds savings, and for any platform/format where mmap
// is unavailable or undesired.
//
// Failure mode: if the file cannot be opened or read, `isOpen()` returns
// false and `bytes()` is empty — no exception thrown so the caller can
// fall back to a default ("file-not-found means no entries").  Callers
// that REQUIRE the file should check `isOpen()` and emit their own error.
// ---------------------------------------------------------------------------
class ORBITAL_CORE_API LoadedFile final : public IDataSource {
   public:
	explicit LoadedFile(std::string filePath);

	LoadedFile(LoadedFile&&) noexcept = default;
	LoadedFile& operator=(LoadedFile&&) noexcept = default;

	[[nodiscard]] raw_file_data bytes() const noexcept override;
	[[nodiscard]] std::string_view path() const noexcept override;

	[[nodiscard]] bool isOpen() const noexcept override { return m_open; }

   private:
	std::vector<std::byte> m_buffer;
	std::string m_path;
	bool m_open{false};
};
