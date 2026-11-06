#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include <vector>

struct ORBITAL_CORE_API SizeEstimate {
	double albedoLo;
	double albedoHi;
	double diameterLoKm;  // at albedoHi (smaller diameter)
	double diameterHiKm;  // at albedoLo (larger diameter)
};

using SizeEstimate_arr = std::vector<SizeEstimate>;

class ORBITAL_CORE_API PhysicalParams {
   public:
	/// Estimate the diameter bounds at fixed H for 0 < albedoLo <= albedoHi.
	static SizeEstimate sizeEstimate(double H, double albedoLo, double albedoHi);

	/// Compute diameter estimates for standard albedo ranges
	static SizeEstimate_arr sizeEstimates(double H);
};
