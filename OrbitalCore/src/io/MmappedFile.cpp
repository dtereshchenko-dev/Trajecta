#include "MmappedFile.h"
#include <string>
#include <utility>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

struct MmappedFile::Impl {
#ifdef _WIN32
	HANDLE fileHandle = INVALID_HANDLE_VALUE;
	HANDLE mappingHandle = nullptr;
	void* view = nullptr;
#else
	int fd = -1;
	void* view = nullptr;
#endif
	std::size_t size = 0;
};

namespace {
#ifdef _WIN32
	[[nodiscard]] std::wstring widenUtf8(const std::string& s) noexcept
	{
		if (s.empty())
			return {};
		const int needed = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
		if (needed <= 0)
			return {};

		std::wstring out(static_cast<std::size_t>(needed), L'\0');
		MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), out.data(), needed);

		return out;
	}
#endif
}  // namespace

MmappedFile::MmappedFile(std::string filePath) : p_impl(new Impl), m_path(std::move(filePath))
{
#ifdef _WIN32
	const std::wstring wide = widenUtf8(m_path);
	p_impl->fileHandle = CreateFileW(wide.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
	                                 FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS, nullptr);

	if (p_impl->fileHandle == INVALID_HANDLE_VALUE)
		return;

	LARGE_INTEGER sizeLi;
	if (!GetFileSizeEx(p_impl->fileHandle, &sizeLi) || sizeLi.QuadPart <= 0) {
		release();
		return;
	}

	p_impl->size = static_cast<std::size_t>(sizeLi.QuadPart);

	p_impl->mappingHandle = CreateFileMappingW(p_impl->fileHandle, nullptr, PAGE_READONLY, 0, 0, nullptr);

	if (p_impl->mappingHandle == nullptr) {
		release();
		return;
	}

	p_impl->view = MapViewOfFile(p_impl->mappingHandle, FILE_MAP_READ, 0, 0, p_impl->size);

	if (p_impl->view == nullptr) {
		release();
		return;
	}
#else
	p_impl->fd = ::open(m_path.c_str(), O_RDONLY);
	if (p_impl->fd < 0)
		return;

	struct stat st{};
	if (::fstat(p_impl->fd, &st) != 0 || st.st_size <= 0) {
		release();
		return;
	}

	p_impl->size = static_cast<std::size_t>(st.st_size);

	p_impl->view = ::mmap(nullptr, p_impl->size, PROT_READ, MAP_PRIVATE, p_impl->fd, 0);
	if (p_impl->view == MAP_FAILED) {
		p_impl->view = nullptr;
		release();
		return;
	}
#endif
}

MmappedFile::~MmappedFile()
{
	release();
	delete p_impl;
}

void MmappedFile::release() noexcept
{
	if (p_impl == nullptr)
		return;

#ifdef _WIN32
	if (p_impl->view) {
		UnmapViewOfFile(p_impl->view);
		p_impl->view = nullptr;
	}

	if (p_impl->mappingHandle) {
		CloseHandle(p_impl->mappingHandle);
		p_impl->mappingHandle = nullptr;
	}

	if (p_impl->fileHandle != INVALID_HANDLE_VALUE) {
		CloseHandle(p_impl->fileHandle);
		p_impl->fileHandle = INVALID_HANDLE_VALUE;
	}

#else
	if (p_impl->view) {
		::munmap(p_impl->view, p_impl->size);
		p_impl->view = nullptr;
	}
	if (p_impl->fd >= 0) {
		::close(p_impl->fd);
		p_impl->fd = -1;
	}
#endif
	p_impl->size = 0;
}

MmappedFile::MmappedFile(MmappedFile&& other) noexcept : p_impl(other.p_impl), m_path(std::move(other.m_path))
{
	other.p_impl = nullptr;
}

MmappedFile& MmappedFile::operator=(MmappedFile&& other) noexcept
{
	if (this == &other)
		return *this;

	release();
	delete p_impl;

	p_impl = other.p_impl;
	m_path = std::move(other.m_path);
	other.p_impl = nullptr;
	return *this;
}

raw_file_data MmappedFile::bytes() const noexcept
{
	if (!isOpen())
		return {};

	return raw_file_data(static_cast<const std::byte*>(p_impl->view), p_impl->size);
}

bool MmappedFile::isOpen() const noexcept { return p_impl != nullptr && p_impl->view != nullptr; }
