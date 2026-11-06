#pragma once

#include "PerturberSelector.h"
#include "orbital_core/FitConfig.h"

// Internal single-pass configuration. Selection/clipping belong to FitSession.
struct OrbitOptimizerConfig {
	orbital_core::NBodySettings nbody;
	orbital_core::NonGravSettings nonGrav;
	// Compatible NG warm start or numerical gravity ready for the next stage.
	bool skipNonGravBootstrap{false};
	bool useWeights{false};
	orbital_core::AstrometricWeightProfile weightProfile{orbital_core::AstrometricWeightProfile::Vfcc17};
	orbital_core::ObservationWeightingOptions obsWeighting;
	std::optional<double> gaiaPhotocenterRadiusKm;
	bool evaluateOnly{false};

	// JD TDB, resolved by the caller; zero retains the fit epoch.
	double outputEpochTdb{0.0};
	bool computeUncertainty{true};
	bool populateResiduals{true};
	bool logDiagnostics{true};
	bool logCacheStats{true};

	// ATT-EXT policy, consumed before constructing the numerical model.
	Calculations::PerturberSelectionMode perturberSelectionMode{Calculations::PerturberSelectionMode::AutoByArcLength};
	std::uint16_t alwaysIncludedPerturberMask{0};
};
