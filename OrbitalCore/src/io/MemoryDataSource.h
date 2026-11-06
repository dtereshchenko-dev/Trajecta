#pragma once

#include "IDataSource.h"
#include <utility>
#include <string>
#include <string_view>
#include <vector>

// ---------------------------------------------------------------------------
// MemoryDataSource
//
// In-memory `IDataSource` impl backed by `std::vector<std::byte>`. Designed
// for unit tests — parsers can be exercised with synthetic byte fixtures
// without touching the filesystem.
// ---------------------------------------------------------------------------
class MemoryDataSource final : public IDataSource {
   public:
	MemoryDataSource(std::vector<std::byte> buffer, std::string name = {}) noexcept
	    : m_buffer(std::move(buffer)), m_name(std::move(name))
	{
		///<
	}

	static MemoryDataSource fromString(std::string_view text, std::string name = {})
	{
		std::vector<std::byte> buf(text.size());
		for (std::size_t i = 0; i < text.size(); ++i)
			buf[i] = static_cast<std::byte>(text[i]);

		return MemoryDataSource(std::move(buf), std::move(name));
	}

	[[nodiscard]] raw_file_data bytes() const noexcept override
	{
		return raw_file_data(m_buffer.data(), m_buffer.size());
	}

	[[nodiscard]] std::string_view path() const noexcept override { return m_name; }

	[[nodiscard]] bool isOpen() const noexcept override
	{
		return true;  // memory buffer is always "open"
	}

	MemoryDataSource(MemoryDataSource&&) noexcept = default;
	MemoryDataSource& operator=(MemoryDataSource&&) noexcept = default;

   private:
	std::vector<std::byte> m_buffer;
	std::string m_name;
};
