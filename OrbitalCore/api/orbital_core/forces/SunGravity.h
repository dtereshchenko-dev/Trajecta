#pragma once

#include "IForceModel.h"
#include <string_view>
#include <cmath>
#include "PointMassHelpers.h"
#include "orbital_core/Constants.h"

namespace orbital_core::forces {
	class SunGravity {
	   public:
		void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& /*dadv*/) const noexcept
		{
			const double r2 = ctx.R.dot(ctx.R);
			if (r2 < 1e-24)
				return;

			const double rMag = std::sqrt(r2);
			const double r3 = r2 * rMag;
			const double gm = (ctx.centralNewtonianGm > 0.0) ? ctx.centralNewtonianGm : Constants::GM_SUN;

			a = a + ctx.R * (-gm / r3);

			if (ctx.computePartials) {
				accumulatePointMass_dR(Vec3D(-ctx.R.x, -ctx.R.y, -ctx.R.z), gm, dadr);
			}
		}

		[[nodiscard]] std::string_view name() const noexcept { return "SunGravity"; }
	};

	static_assert(ForceModel<SunGravity>);
}  // namespace orbital_core::forces
