#include "orbital_core/uncertainty/DerivedSigmas.h"
#include "orbital_core/Constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbital_core::uncertainty {
	DerivedSigmas DerivedSigmasComputer::compute(const Orbit& orbit, const ElementCov6x6& covOrb)
	{
		DerivedSigmas s;

		constexpr double kFdStep = 1.0e-7;
		constexpr double kDaysYr = 365.25;

		const double a = orbit.getSemiMajorAxis();
		const double e = orbit.getExcentricity();
		const double nu = orbit.getTrueAnomaly();
		const double absA = std::fabs(a);
		const bool isElliptical = (orbit.getType() == OrbitType::Elliptical);

		s.sigma_a = std::sqrt(covOrb(0, 0));
		s.sigma_e = std::sqrt(covOrb(1, 1));
		s.sigma_i = std::sqrt(covOrb(2, 2));
		s.sigma_Om = std::sqrt(covOrb(3, 3));
		s.sigma_w = std::sqrt(covOrb(4, 4));
		s.sigma_nu = std::sqrt(covOrb(5, 5));
		s.sigma_inv_a = s.sigma_a / (a * a);

		auto computeM = [](double ee, double nuu) noexcept {
			if (ee < 1.0) {
				const double tanHalfE = std::sqrt((1.0 - ee) / (1.0 + ee)) * std::tan(0.5 * nuu);
				const double E = 2.0 * std::atan(tanHalfE);
				return E - ee * std::sin(E);
			}

			else {
				constexpr double kClamp = 1.0 - 1e-12;
				const double arg = std::sqrt((ee - 1.0) / (ee + 1.0)) * std::tan(0.5 * nuu);
				const double argClamped = std::max(-kClamp, std::min(kClamp, arg));
				const double H = 2.0 * std::atanh(argClamped);
				return ee * std::sinh(H) - H;
			}
		};

		const double dM_de = (computeM(e + kFdStep, nu) - computeM(e - kFdStep, nu)) / (2.0 * kFdStep);
		const double dM_dnu = (computeM(e, nu + kFdStep) - computeM(e, nu - kFdStep)) / (2.0 * kFdStep);
		const double sigmaM_sq =
		    dM_de * dM_de * covOrb(1, 1) + dM_dnu * dM_dnu * covOrb(5, 5) + 2.0 * dM_de * dM_dnu * covOrb(1, 5);
		s.sigma_M = (sigmaM_sq > 0.0) ? std::sqrt(sigmaM_sq) : 0.0;

		const double per_days = 2.0 * Constants::PI * std::sqrt(absA * absA * absA / Constants::GM_SUN);
		if (isElliptical && a > 0.0) {
			const double invN = per_days / (2.0 * Constants::PI);
			Eigen::Matrix<double, 1, 6> dTp = Eigen::Matrix<double, 1, 6>::Zero();
			dTp(0, 0) = 1.5 * (orbit.getTimeOfPerihelion() - orbit.getEpoch()) / a;
			dTp(0, 1) = -dM_de * invN;
			dTp(0, 5) = -dM_dnu * invN;

			const double sigmaTpSq = (dTp * covOrb * dTp.transpose())(0, 0);
			s.sigma_Tp = (sigmaTpSq > 0.0) ? std::sqrt(sigmaTpSq) : 0.0;
		}

		else {
			s.sigma_Tp = s.sigma_M * per_days / (2.0 * Constants::PI);
		}

		{
			const double dq_da = 1.0 - e;
			const double dq_de = -a;
			const double sq =
			    dq_da * dq_da * covOrb(0, 0) + dq_de * dq_de * covOrb(1, 1) + 2.0 * dq_da * dq_de * covOrb(0, 1);
			s.sigma_q = (sq > 0.0) ? std::sqrt(sq) : 0.0;
		}

		if (isElliptical) {
			const double dQ_da = 1.0 + e;
			const double dQ_de = a;
			const double sQ =
			    dQ_da * dQ_da * covOrb(0, 0) + dQ_de * dQ_de * covOrb(1, 1) + 2.0 * dQ_da * dQ_de * covOrb(0, 1);
			s.sigma_Q = (sQ > 0.0) ? std::sqrt(sQ) : 0.0;
		}

		else {
			s.sigma_Q = std::numeric_limits<double>::quiet_NaN();
		}

		const double n_radPerDay = std::sqrt(Constants::GM_SUN / (absA * absA * absA));
		s.sigma_n = 1.5 * n_radPerDay * s.sigma_a / absA;

		if (isElliptical) {
			s.sigma_P = 1.5 * per_days * s.sigma_a / absA;
			s.sigma_Py = s.sigma_P / kDaysYr;
		}

		else {
			s.sigma_P = std::numeric_limits<double>::quiet_NaN();
			s.sigma_Py = std::numeric_limits<double>::quiet_NaN();
		}

		return s;
	}
}  // namespace orbital_core::uncertainty
