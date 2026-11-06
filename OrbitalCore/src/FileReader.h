#ifndef FILEREADER_H
#define FILEREADER_H
#include <memory>
#include <filesystem>
#include <string>
#include "orbital_core/OrbitalCoreExport.h"

class IDataSource;

// ---------------------------------------------------------------------------
// FileReader — small text/binary file loader.  Thin wrapper around the
// `io/IDataSource` layer (Sprint 1) so every existing consumer (MPECParser,
// AstrometryBias, the OBSCODE.dat loader, weights/vfcc17.rules,
// …) routes through the same factory-selected backend (`MmappedFile` for
// large binaries, `LoadedFile` for small text).  Public API stays
// unchanged — `read(filename)` + `data()` + `error()` — so the migration
// is invisible to callers.
// ---------------------------------------------------------------------------
class ORBITAL_CORE_API FileReader {
   public:
	enum class Error { NoError, FileNotFound, ReadError };

	FileReader();
	~FileReader();

	FileReader(const FileReader&) = delete;
	FileReader& operator=(const FileReader&) = delete;

	// Takes a path, not a raw character pointer: path::value_type is wchar_t
	// on Windows and char elsewhere, so a pointer parameter fixes the build to
	// one platform. Every existing call site still compiles -- path converts
	// implicitly from both spellings.
	void read(const std::filesystem::path& filename);
	std::string data() const;

	Error error() const noexcept;
	std::string errorString() const noexcept;

   private:
	Error m_error{Error::NoError};
	std::string m_data;
	std::unique_ptr<IDataSource> p_source;
};

#endif  // FILEREADER_H
