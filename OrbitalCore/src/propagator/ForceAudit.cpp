#include "orbital_core/forces/ForceAuditControl.h"

#include "propagator/ForceAudit.h"

namespace orbital_core::forces::audit {

	void arm(const char* phase) noexcept { armInternal(phase ? phase : "unnamed"); }

}  // namespace orbital_core::forces::audit
