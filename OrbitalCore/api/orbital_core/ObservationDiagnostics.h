#pragma once

#include "orbital_core/ObservationIdentity.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <span>
#include <vector>

namespace orbital_core {

	enum class ObservationIssue {
		CoincidentPosition,
		NearbyObservationConflict,
		BelowHorizon,
		Daylight,
	};

	struct ObservationVisibilityDiagnostic {
		ObservationId observation;
		std::vector<ObservationIssue> issues;
		std::optional<double> objectAltitudeDeg;
		std::optional<double> sunAltitudeDeg;
	};

	struct ObservationTimeGroup {
		std::vector<ObservationId> observations;
		std::vector<ObservationIssue> issues;
		// The largest separation from the first record, not a fit residual.
		double maxSeparationArcsec{};
	};

	// Reports exact-position duplicate candidates and groups conflicting under
	// the three-second / two-arcsecond coincidence check. No fit mask is changed.
	[[nodiscard]] ORBITAL_CORE_API std::vector<ObservationTimeGroup> findCoincidentObservations(
	    std::span<const ObsNode> observations);

	// Geometric altitudes use the submitted epoch and ignore refraction.
	// Geocentric and spacecraft records have no terrestrial horizon.
	[[nodiscard]] ORBITAL_CORE_API std::vector<ObservationVisibilityDiagnostic> diagnoseObservationVisibility(
	    std::span<const ObsNode> observations);

}  // namespace orbital_core
