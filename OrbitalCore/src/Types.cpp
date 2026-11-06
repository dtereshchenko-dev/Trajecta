#include "orbital_core/Types.h"
#include <filesystem>
#include "orbital_core/IEphemerisSource.h"

std::filesystem::path& getExecutableDir()
{
	static std::filesystem::path executableDir;
	return executableDir;
}

static IEphemerisSource* s_defaultEphemerisSource = nullptr;

IEphemerisSource* getDefaultEphemerisSource() { return s_defaultEphemerisSource; }

void setDefaultEphemerisSource(IEphemerisSource* source) { s_defaultEphemerisSource = source; }
