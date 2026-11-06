#pragma once

#include "JobFiles.h"

namespace strewn::jobfiles {

	Json downloadEnvironment(const EnvironmentUrls& sources, const std::filesystem::path& directory,
	                         std::stop_token cancel);

}
