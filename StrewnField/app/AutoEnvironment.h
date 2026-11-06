#pragma once

#include "GfsSource.h"
#include <span>

namespace strewn::sources {

	void validateAutomatic(const AutomaticEnvironment& options);
	Bounds environmentBounds(const AutomaticEnvironment& options, std::span<const EntryState> entries);
	Json automaticEnvironment(const AutomaticEnvironment& options, std::span<const EntryState> entries,
	                          const std::filesystem::path& directory, std::stop_token cancel);

}  // namespace strewn::sources
