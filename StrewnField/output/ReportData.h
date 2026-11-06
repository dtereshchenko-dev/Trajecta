#pragma once

#include "FieldOutputInternal.h"
#include <span>
#include <string_view>

namespace strewn::output {

	std::string base64(std::string_view bytes);
	std::string encodeText(std::string_view bytes);
	// Visualization payload: float32 LE, byte-plane shuffle, zlib Deflate, base64.
	std::string encodeFloats(std::span<const double> values);
	Json buildReportData(const std::filesystem::path& input, const ArchivedField& archive, ArtifactReport& report,
	                     std::stop_token cancel);

}  // namespace strewn::output
