#pragma once

#include <functional>
#include <stop_token>
#include <string>

namespace orbital_core {
	struct FitRunControl {
		std::stop_token cancel;
		std::function<void(const std::string&)> progress;
	};
}  // namespace orbital_core
