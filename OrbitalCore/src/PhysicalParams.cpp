#include "PhysicalParams.h"
#include "orbital_core/hazard/ImpactEnergy.h"
#include <vector>

SizeEstimate PhysicalParams::sizeEstimate(double H, double albedoLo, double albedoHi)
{
	return {albedoLo, albedoHi, orbital_core::hazard::diameterFromMagnitudeKm(H, albedoHi),
	        orbital_core::hazard::diameterFromMagnitudeKm(H, albedoLo)};
}

SizeEstimate_arr PhysicalParams::sizeEstimates(double H)
{
	static const struct {
		double lo, hi;
	} ranges[] = {
	    {0.05, 0.10}, {0.10, 0.15}, {0.15, 0.25}, {0.25, 0.50}, {0.50, 1.00},
	};

	std::vector<SizeEstimate> result;
	for (const auto& r : ranges) {
		result.push_back(sizeEstimate(H, r.lo, r.hi));
	}
	return result;
}
