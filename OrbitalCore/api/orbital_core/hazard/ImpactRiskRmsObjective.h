#pragma once

#include "orbital_core/FitConfig.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/hazard/ImpactRisk.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <cstddef>
#include <span>

namespace orbital_core::hazard {

	// `fitObservations` is the set the fit ran on (FitPipeline gives FitSession the prepared
	// ObservationSet `all`) and `keptIndices` the fit's selection (FitResult::keptObsIndices;
	// empty keeps every row).  The objective prepares that set once and evaluates the kept
	// rows, so its weights, including the vfcc17 same-station relaxation that counts batches
	// over the whole set, are the weights of the fit.
	[[nodiscard]] ORBITAL_CORE_API ImpactRiskRmsScaling
	makeImpactRiskRmsScaling(Orbit_p fitEpochOrbit, std::span<const ObsNode> fitObservations,
	                         std::span<const std::size_t> keptIndices, const FitConfig& fitConfig,
	                         const NonGravModel& nominalNonGrav, double objectiveScale = 1.0, int numSamples = 1000);

}  // namespace orbital_core::hazard
