#include "FileReader.h"
#include <string>

#include <algorithm>
#include "io/DataSourceFactory.h"
#include "io/IDataSource.h"
#include <filesystem>

FileReader::FileReader() = default;
FileReader::~FileReader() = default;

void FileReader::read(const std::filesystem::path& filename)
{
	p_source.reset();
	m_data.clear();

	const std::string narrowPath = filename.string();

	if (narrowPath.empty()) {
		m_error = Error::FileNotFound;
		return;
	}

	p_source = DataSourceFactory::open(narrowPath);
	if (!p_source || !p_source->isOpen()) {
		m_error = std::filesystem::exists(narrowPath) ? Error::ReadError : Error::FileNotFound;
		p_source.reset();
		return;
	}

	const auto view = p_source->text();
	m_data.assign(view.data(), view.size());
	m_data.erase(std::remove(m_data.begin(), m_data.end(), '\r'), m_data.end());

	m_error = Error::NoError;
}

std::string FileReader::data() const { return m_data; }

FileReader::Error FileReader::error() const noexcept { return m_error; }

std::string FileReader::errorString() const noexcept
{
	switch (m_error) {
		case Error::FileNotFound:
			return "File not found";
		case Error::ReadError:
			return "Error reading file";
		default:
			return "";
	}
}
