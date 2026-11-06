#pragma once

#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/OrbitalCoreExport.h"

class ObsNode;

namespace orbital_core::obs {
	// Loads native weights/vfcc17.rules next to the executable once.
	// All lookups are thread-safe; vfcc17 uses station, technique, catalog,
	// program and UTC date with reference precedence and (start, end] intervals.
	ORBITAL_CORE_API void ensureSigmaTable(AstrometricWeightProfile profile = AstrometricWeightProfile::Vfcc17);

	[[nodiscard]] ORBITAL_CORE_API ObservationUncertainty lookupUncertainty(
	    const ObsNode& obs, double jd, AstrometricWeightProfile profile = AstrometricWeightProfile::Vfcc17);

	// Whether the observer's ADES rmsRA/rmsDec replace the table sigma.  They do
	// when enabled, and always for occultation astrometry ('E'): no station
	// table describes chord timings at the milliarcsecond level. Deprecated
	// and replaced records ('X', 'x') keep their exclusion weight either way.
	[[nodiscard]] ORBITAL_CORE_API bool usesObserverSigma(const ObsNode& obs, bool enabled) noexcept;

}  // namespace orbital_core::obs
