#pragma once

#include "strewn_core/StrewnJob.h"
#include <functional>

namespace strewn {
	JobResult runConfiguredJob(const std::filesystem::path& config, const RunControl& control = {},
	                           std::function<void(int, int)> orbitalProgress = {});

}  // namespace strewn
