#pragma once

#include "IForceModel.h"
#include "orbital_core/Constants.h"

#include <cmath>
#include <string_view>

namespace orbital_core::forces {

	class SunJ2 {
	   public:
		// DE440 J2/reference-radius pair. Pole RA/Dec are in the J2000 equator;
		// the force and its position derivative use the same ecliptic J2000 axis.
		static constexpr double kJ2 = 2.1961391516529825e-7;
		static constexpr double kRadiusAu = 696000.0 / 149597870.7;

		[[nodiscard]] static const Vec3D& poleEcliptic() noexcept
		{
			static const Vec3D pole = [] {
				const double ra = 286.13 * Constants::DEG_TO_RAD;
				const double dec = 63.87 * Constants::DEG_TO_RAD;
				const double y = std::cos(dec) * std::sin(ra);
				const double z = std::sin(dec);
				const double c = std::cos(Constants::EARTH_INCL);
				const double s = std::sin(Constants::EARTH_INCL);
				return Vec3D(std::cos(dec) * std::cos(ra), c * y + s * z, -s * y + c * z);
			}();
			return pole;
		}

		void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& /*dadv*/) const noexcept
		{
			const double r2 = ctx.R.dot(ctx.R);
			if (r2 < 1.0e-24)
				return;

			const double r = std::sqrt(r2);
			const double r5 = r2 * r2 * r;
			const double sunGm = ctx.sunGm > 0.0 ? ctx.sunGm : Constants::GM_SUN;

			const Vec3D& pole = poleEcliptic();
			const double z = ctx.R.dot(pole);
			const double s2 = z * z / r2;
			const double ratio = 1.5 * sunGm * kJ2 * kRadiusAu * kRadiusAu / r5;
			const double bracket = 5.0 * s2 - 1.0;
			a += (ctx.R * bracket - pole * (2.0 * z)) * ratio;
			if (ctx.computePartials) {
				for (int i = 0; i < 3; ++i) {
					for (int j = 0; j < 3; ++j) {
						dadr(i, j) +=
						    ratio *
						    ((i == j ? bracket : 0.0) + (10.0 * z / r2) * (ctx.R[i] * pole[j] + pole[i] * ctx.R[j]) +
						     ((5.0 - 35.0 * s2) / r2) * ctx.R[i] * ctx.R[j] - 2.0 * pole[i] * pole[j]);
					}
				}
			}
		}

		[[nodiscard]] std::string_view name() const noexcept { return "SunJ2"; }
	};

	static_assert(ForceModel<SunJ2>);

}  // namespace orbital_core::forces
