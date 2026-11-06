#pragma once

#include "fit/OrbitFitConfig.h"
#include "obs/PreparedObservation.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <cstdint>
#include <span>
#include <vector>

namespace orbital_core::fit {

	// Fit-ready observations: the per-observation preparation of
	// ObservationPipeline followed by the vfcc17 same-station relaxation.  The
	// relaxation counts batches over the whole set it is given, so the set that
	// defines the batches is prepared once and every fit of a subset selects rows
	// from it; preparing a subset again would recount its batches.  Input must be
	// JD-sorted; the output keeps the input order.
	[[nodiscard]] ORBITAL_CORE_API std::vector<obs::PreparedObservation> prepareFitObservations(
	    std::span<const ObsNode> observations, const OrbitOptimizerConfig& config);

	// Relaxation counted over exactly `observations`; active for vfcc17 with
	// observation weighting on.  Logs a summary when config.logDiagnostics is set.
	ORBITAL_CORE_API void relaxSameStationBatches(std::vector<obs::PreparedObservation>& observations,
	                                              const OrbitOptimizerConfig& config);

	// Rows whose mask entry is non-zero, in input order.
	[[nodiscard]] ORBITAL_CORE_API std::vector<obs::PreparedObservation> selectObservations(
	    std::span<const obs::PreparedObservation> observations, std::span<const std::uint8_t> included);

}  // namespace orbital_core::fit
