#pragma once

#include "IDataSource.h"
#include <string_view>
#include <string>

// ---------------------------------------------------------------------------
// MmappedFile
//
// Memory-mapped read-only file.  RAII: opens the file + maps the region
// in the constructor, releases both in the destructor.  Move-only.
//
// Use for large binary kernels (DE440, SB441) where the LM-loop
// hot path performs millions of random-access reads — mmap eliminates the
// per-read syscall and lets the OS page-cache amortise across worker
// threads.  For small files (<16 MB by default) prefer `LoadedFile` —
// mmap setup cost outweighs the saving.
// ---------------------------------------------------------------------------
class ORBITAL_CORE_API MmappedFile final : public IDataSource {
   public:
	explicit MmappedFile(std::string filePath);
	~MmappedFile() override;

	MmappedFile(MmappedFile&&) noexcept;
	MmappedFile& operator=(MmappedFile&&) noexcept;

	[[nodiscard]] raw_file_data bytes() const noexcept override;

	[[nodiscard]] std::string_view path() const noexcept override { return m_path; }

	[[nodiscard]] bool isOpen() const noexcept override;

   private:
	struct Impl;
	Impl* p_impl{nullptr};
	std::string m_path;

	void release() noexcept;
};
