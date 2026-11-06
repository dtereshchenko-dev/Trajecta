#pragma once

#include "orbital_core/ObsNode.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <istream>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>

namespace orbital_core::obs {

	using ObservatoryCatalog = std::unordered_map<std::string, ObservatoryRecord>;

	// Internal parsed RADCODE v2: one geometry per station, indices derived from
	// its JPL/MPC columns. '-' denotes an identifier not established by a source.
	struct RadarCatalog {
		ObservatoryCatalog stations;
		std::unordered_map<std::string, std::string> aliases;
	};

	// Requires '! RADCODE format=2', followed by rows:
	// code jpl_code mpc_code longitude_deg rho_xy_1e-10AU rho_z_1e-10AU name.
	// A JPL-only station can use its negative JPL identifier as the canonical code.
	// Rejects malformed rows, duplicate identifiers and ambiguous namespaces.
	// No network access; coordinates are converted to Earth-radius units once.
	[[nodiscard]] ORBITAL_CORE_API RadarCatalog readRadarCatalog(std::istream& stream);

	// MPC observatory list (`OBSCODE.dat` in the configured data directory), parsed on first
	// use and immutable for the rest of the process.  The observation loader
	// (geodetic data copied into ObsNode) and ObservationPipeline (observer
	// geometry for residuals) both resolve stations here, so the file is read
	// once rather than once per fit pass.  Throws std::runtime_error when the
	// file is missing or empty; a later call retries the load.
	[[nodiscard]] ORBITAL_CORE_API const ObservatoryCatalog& observatoryCatalog();
	[[nodiscard]] std::filesystem::path activeOpticalCatalogPath();

	// Dedicated radar coordinates from RADCODE.dat. Returned parallax constants
	// use the same Earth-radius units as ObservatoryRecord. Unknown radar stations
	// must be added to this data file rather than silently using optical rounding.
	[[nodiscard]] ORBITAL_CORE_API const ObservatoryCatalog& radarObservatoryCatalog();
	[[nodiscard]] std::filesystem::path activeRadarCatalogPath();

	// Resolves a canonical RADCODE key or its JPL/MPC identifier from the same
	// station row. All identifiers return the same canonical record.
	// The returned record carries the canonical RADCODE key. Legacy radar keys
	// may differ from current MPC codes; this lookup must not serve optical data.
	[[nodiscard]] ORBITAL_CORE_API const ObservatoryRecord& radarObservatory(std::string_view code);

}  // namespace orbital_core::obs
