#pragma once

#include "orbital_core/FitRunControl.h"
#include "orbital_core/hazard/CloseApproach.h"

namespace orbital_core::hazard {
	// Checked, cancellable variant of the shared nominal minimum finder. Epochs
	// have the same JD TDB convention as CloseApproachEvent::jd.
	std::vector<CloseApproachEvent> findCloseApproachesChecked(const Vec3D& r, const Vec3D& v, double epochJdTdb,
	                                                           const NBodyPropagator& propagator,
	                                                           const NonGravModel& ng,
	                                                           const CloseApproachConfig& config, IEphemerisSource* eph,
	                                                           const FitRunControl& control);
}  // namespace orbital_core::hazard
