#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <string_view>
#include <span>

// ---------------------------------------------------------------------------
// IDataSource
//
// Uniform abstraction for read-only file access across OrbitalCore.  Every
// concrete reader (`SpkReader`, `FileReader`, …) consumes
// data through this interface so we can:
//
//   * memory-map large binaries (DE440 / SB441) for zero-copy hot
//     paths instead of using `fread`+`fseek` syscalls per Chebyshev block;
//   * route small text resources (OBSCODE.dat, weights/vfcc17.rules) through
//     the same factory without forcing them onto mmap (cheap full-load);
// ---------------------------------------------------------------------------

using raw_file_data = std::span<const std::byte>;

class ORBITAL_CORE_API IDataSource {
   public:
	virtual ~IDataSource() = default;

	IDataSource(const IDataSource&) = delete;
	IDataSource& operator=(const IDataSource&) = delete;

	[[nodiscard]] virtual raw_file_data bytes() const noexcept = 0;
	[[nodiscard]] virtual std::string_view text() const noexcept;
	[[nodiscard]] virtual std::string_view path() const noexcept = 0;
	[[nodiscard]] virtual std::size_t size() const noexcept;
	[[nodiscard]] virtual bool isOpen() const noexcept = 0;

   protected:
	IDataSource() = default;
	IDataSource(IDataSource&&) noexcept = default;
	IDataSource& operator=(IDataSource&&) noexcept = default;
};
