#pragma once

#include "orbital_core/Orbit.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <Eigen/Core>

namespace orbital_core::uncertainty {

	using ElementCov6x6 = Eigen::Matrix<double, 6, 6, Eigen::DontAlign>;

	struct DerivedSigmas {
		// ---- Raw element sigmas (sqrt of covOrb diagonal) -----
		double sigma_a{0.0};   // semi-major axis, AU
		double sigma_e{0.0};   // eccentricity (dimensionless)
		double sigma_i{0.0};   // inclination, rad
		double sigma_Om{0.0};  // longitude of ascending node, rad
		double sigma_w{0.0};   // argument of pericentre, rad
		double sigma_nu{0.0};  // true anomaly, rad

		// ---- Derived (linear error propagation through analytic Jacobians) ----
		double sigma_inv_a{0.0};  // 1/a, per-AU
		double sigma_M{0.0};      // mean anomaly, rad — propagated from σ_e, σ_ν, cov(e,ν)
		double sigma_Tp{0.0};     // time of perihelion, days
		double sigma_q{0.0};      // perihelion distance q = a(1-e), AU
		double sigma_Q{0.0};      // apogee Q = a(1+e), AU — NaN for hyperbolic
		double sigma_n{0.0};      // mean motion, rad/day
		double sigma_P{0.0};      // orbital period, days — NaN for hyperbolic
		double sigma_Py{0.0};     // orbital period, years — NaN for hyperbolic
	};

	class ORBITAL_CORE_API DerivedSigmasComputer {
	   public:
		[[nodiscard]] static DerivedSigmas compute(const Orbit& orbit, const ElementCov6x6& covOrb);
	};
}  // namespace orbital_core::uncertainty
