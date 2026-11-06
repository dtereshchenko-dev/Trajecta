#pragma once
#include "obs/IObservationReader.h"

namespace orbital_core::obs {
	// Reads RWO v2 and ESA NEOCC v3. JPL radar IDs resolve through RADCODE
	// station identities; unknown IDs fail instead of losing geometry.
	// Station vocabulary: https://ssd-api.jpl.nasa.gov/doc/sb_radar.html
	class ORBITAL_CORE_API RwoReader final : public IObservationReader {
	   public:
		const char* formatName() const noexcept override { return "RWO"; }

		bool canRead(std::string_view text) const noexcept override;
		std::vector<ParsedObservation> read(const std::string& text, std::string& error) const override;
	};
}  // namespace orbital_core::obs
