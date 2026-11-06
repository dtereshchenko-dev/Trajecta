#pragma once

#include "orbital_core/FitResult.h"
#include <filesystem>
#include <optional>
#include <string>

namespace orbital_core {

	struct OrbitSolution {
		std::string objectName;
		Orbit nominal;
		// Already scaled, physical layout [R AU, V AU/day, fitted non-grav].
		FitParameterCovariance covariance;
		NonGravModel nonGrav;
		std::optional<double> absoluteMagnitude;
		double slopeG{0.15};
	};

	// Reads standalone OEF2.0 ML (CAR/EQU/KEP/COT), including LSP/NGR
	// and full cross-covariance, plus Core's OC-OEF1.0 for nonstandard force laws.
	// OC-OEF1.0 requires a reader that supports this extension.
	// A directory selects orbit_cartesian_oef.txt.
	// Epochs are converted to model TDB. Unsupported models/frames and incomplete
	// matrices fail explicitly; sidecars are not read.
	// No fitting, covariance rescaling or repair is performed.
	[[nodiscard]] ORBITAL_CORE_API OrbitSolution loadOrbitSolution(const std::filesystem::path& directory);

}  // namespace orbital_core
