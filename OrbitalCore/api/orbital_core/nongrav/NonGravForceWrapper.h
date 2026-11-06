#pragma once

#include "orbital_core/forces/IForceModel.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <Eigen/Core>
#include <string_view>

namespace orbital_core::forces {

	// ForceModel-concept adapter for NonGravModel.  Reads the per-call non-grav
	// pointer carried by `ForceContext`, adds the resulting acceleration to the
	// running accumulators. State partials drive both STM and sensitivity transport;
	// direct A/DT partials are supplied separately by the sensitivity forcing.
	class NonGravForceWrapper {
	   public:
		void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& dadv) const noexcept
		{
			if (!ctx.nonGrav)
				return;
			a = a + ctx.nonGrav->acceleration(ctx.R, ctx.V, ctx.jd);
			if (ctx.computePartials) {
				Eigen::Matrix3d dr, dv;
				ctx.nonGrav->accelerationPartials(ctx.R, ctx.V, dr, dv);
				dadr += dr;
				dadv += dv;
			}
		}

		[[nodiscard]] std::string_view name() const noexcept { return "NonGravForceWrapper"; }
	};

	static_assert(ForceModel<NonGravForceWrapper>);

}  // namespace orbital_core::forces
