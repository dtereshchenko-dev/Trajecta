#include "IDataSource.h"
#include <string_view>

std::string_view IDataSource::text() const noexcept
{
	const auto raw = bytes();
	return std::string_view(reinterpret_cast<const char*>(raw.data()), raw.size());
}

std::size_t IDataSource::size() const noexcept { return bytes().size(); }
