#pragma once

#include "IForceModel.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <string_view>
#include <Eigen/Core>

namespace orbital_core::forces {
	class EarthJ2 {
	   public:
		ORBITAL_CORE_API void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr,
		                               Eigen::Matrix3d& /*dadv*/) const noexcept;

		[[nodiscard]] std::string_view name() const noexcept { return "EarthJ2"; }

	   private:
		[[nodiscard]] ORBITAL_CORE_API static Eigen::Matrix3d equatorialToEclipticAtMjd(double mjdTdb) noexcept;
	};

	static_assert(ForceModel<EarthJ2>);
}  // namespace orbital_core::forces
