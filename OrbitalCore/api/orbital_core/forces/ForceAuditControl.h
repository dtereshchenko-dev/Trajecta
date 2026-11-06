#pragma once

#include "orbital_core/OrbitalCoreExport.h"

// Control surface for the per-force acceleration audit. The audit itself lives
// inside the library (src/propagator/ForceAudit.h); only phase arming needs to
// be reachable from a driver, and it must cross the DLL boundary as a real
// exported symbol so it acts on the library's own budgets.
namespace orbital_core::forces::audit {

	// Re-grants the per-phase record budgets and labels subsequent records
	// with `phase`. No effect unless ORBITAL_FORCE_AUDIT=1.
	ORBITAL_CORE_API void arm(const char* phase) noexcept;
}  // namespace orbital_core::forces::audit
