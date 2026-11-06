#include "LoadedFile.h"
#include <string>
#include <string_view>
#include <utility>
#include <cstdio>

LoadedFile::LoadedFile(std::string filePath) : m_path(std::move(filePath))
{
#ifdef _MSC_VER
	std::FILE* fp = nullptr;
	if (fopen_s(&fp, m_path.c_str(), "rb") != 0 || fp == nullptr)
		return;
#else
	std::FILE* fp = std::fopen(m_path.c_str(), "rb");
	if (!fp)
		return;
#endif

	if (std::fseek(fp, 0, SEEK_END) != 0) {
		std::fclose(fp);
		return;
	}

	const long size = std::ftell(fp);
	if (size < 0) {
		std::fclose(fp);
		return;
	}

	std::rewind(fp);

	m_buffer.resize(static_cast<std::size_t>(size));
	const std::size_t read = std::fread(m_buffer.data(), 1, m_buffer.size(), fp);
	std::fclose(fp);

	if (read != m_buffer.size()) {
		m_buffer.clear();
		return;
	}

	m_open = true;
}

raw_file_data LoadedFile::bytes() const noexcept { return raw_file_data(m_buffer.data(), m_buffer.size()); }

std::string_view LoadedFile::path() const noexcept { return m_path; }
