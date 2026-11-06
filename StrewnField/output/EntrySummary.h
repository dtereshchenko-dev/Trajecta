#pragma once

#include "strewn_core/Environment.h"
#include <nlohmann/json.hpp>
#include <span>

namespace strewn::output {

	// Original entry weights, normalized conditional on entry, never atmospheric run weights.
	// Population spreads use wrapped residuals for angles. A single entry has no estimated sigma.
	// Horizontal ellipse axes are sqrt(eigenvalues) of local EN covariance at entry height;
	// they do not describe a 68% joint confidence contour or uncertainty of the mean.
	nlohmann::json summarizeEntries(std::span<const EntryState> entries);

}  // namespace strewn::output
