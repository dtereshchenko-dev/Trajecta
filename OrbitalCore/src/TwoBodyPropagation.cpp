#include "TwoBodyPropagation.h"
#include "TwoBodyStumpff.h"
#include "orbital_core/Constants.h"

#include <cmath>
#include <stdexcept>

namespace Calculations {
	namespace {
		constexpr double kCostSentinel = 1.0e30;
		constexpr int kKeplerMaxIter = 80;
		constexpr double kKeplerConvTol = 1.0e-11;
		constexpr double kEllipticAlphaThreshold = 1.0e-9;
		constexpr double kHyperbolicAlphaThreshold = -1.0e-9;
		constexpr double kMinStateMagAU = 1.0e-12;
		constexpr double kMinKeplerR = 1.0e-14;
	}  // namespace

	namespace detail {
		[[nodiscard]] bool keplerPropUV(const Vec3D& R0, const Vec3D& V0, double tau, Vec3D& Rout, Vec3D& Vout,
		                                const StateVectors* initialRate, double tauRate, Vec3D* positionRate,
		                                Vec3D* velocityRate) noexcept
		{
			const double r0 = R0.magnitude();
			if (r0 < kMinStateMagAU) {
				return false;
			}

			const double sigma0 = R0.dot(V0);
			const double v02 = V0.dot(V0);
			const double alpha = 2.0 / r0 - v02;

			double chi;
			// The hyperbolic logarithmic seed is ill-conditioned for short steps far from the Sun.
			if (std::abs(tau) * std::sqrt(v02 + 1.0 / r0) < 0.5 * r0) {
				chi = tau / r0;
			} else if (alpha > kEllipticAlphaThreshold) {
				chi = tau * alpha;
			} else if (alpha < kHyperbolicAlphaThreshold) {
				const double a = 1.0 / alpha;
				const double sign = (tau >= 0.0) ? 1.0 : -1.0;
				const double num = -2.0 * alpha * tau;
				const double den = sigma0 + sign * std::sqrt(-a) * (1.0 - r0 * alpha);

				chi = (den != 0.0 && num / den > 0.0) ? sign * std::sqrt(-a) * std::log(num / den) : tau / r0;
			} else {
				chi = tau / r0;
			}

			double c2 = 0.5;
			double c3 = 1.0 / 6.0;
			double r = r0;
			bool converged = false;

			for (int iter = 0; iter < kKeplerMaxIter; ++iter) {

				const double psi = chi * chi * alpha;
				c2 = stumpffC(psi);
				c3 = stumpffS(psi);

				r = chi * chi * c2 + sigma0 * chi * (1.0 - psi * c3) + r0 * (1.0 - psi * c2);

				const double F = sigma0 * chi * chi * c2 + (1.0 - r0 * alpha) * chi * chi * chi * c3 + r0 * chi - tau;

				if (std::fabs(F) < kKeplerConvTol) {
					converged = true;
					break;
				}

				if (std::fabs(r) < kMinKeplerR) {
					return false;
				}

				chi -= F / r;
			}

			if (!converged) {
				return false;
			}

			const double psi = chi * chi * alpha;
			c2 = stumpffC(psi);
			c3 = stumpffS(psi);

			const double f = 1.0 - (chi * chi / r0) * c2;
			const double g = tau - chi * chi * chi * c3;
			const double gd = 1.0 - (chi * chi / r) * c2;
			const double fd = (chi / (r * r0)) * (psi * c3 - 1.0);

			Rout = R0 * f + V0 * g;
			Vout = R0 * fd + V0 * gd;

			if (initialRate && (positionRate || velocityRate)) {
				const double dr0 = R0.dot(initialRate->first) / r0;
				const double dsigma = initialRate->first.dot(V0) + R0.dot(initialRate->second);
				const double dalpha = -2.0 * dr0 / (r0 * r0) - 2.0 * V0.dot(initialRate->second);
				const double dc2 = std::abs(psi) <= kStumpffTaylorPsiThreshold
				                       ? -1.0 / 24.0 + psi / 360.0
				                       : (1.0 - psi * c3 - 2.0 * c2) / (2.0 * psi);
				const double dc3 = std::abs(psi) <= kStumpffTaylorPsiThreshold ? -1.0 / 120.0 + psi / 2520.0
				                                                               : (c2 - 3.0 * c3) / (2.0 * psi);
				const double partial = dsigma * chi * chi * c2 + sigma0 * chi * chi * dc2 * chi * chi * dalpha -
				                       (dr0 * alpha + r0 * dalpha) * chi * chi * chi * c3 +
				                       (1.0 - r0 * alpha) * chi * chi * chi * dc3 * chi * chi * dalpha + dr0 * chi -
				                       tauRate;
				const double dchi = -partial / r;
				const double dpsi = 2.0 * chi * dchi * alpha + chi * chi * dalpha;
				const double df =
				    -(2.0 * chi * dchi * c2 + chi * chi * dc2 * dpsi) / r0 + chi * chi * c2 * dr0 / (r0 * r0);
				const double dg = tauRate - 3.0 * chi * chi * dchi * c3 - chi * chi * chi * dc3 * dpsi;
				const Vec3D dposition = initialRate->first * f + R0 * df + initialRate->second * g + V0 * dg;
				if (positionRate)
					*positionRate = dposition;
				if (velocityRate) {
					const double dr = Rout.dot(dposition) / r;
					const double dgd =
					    -(2.0 * chi * dchi * c2 + chi * chi * dc2 * dpsi) / r + chi * chi * c2 * dr / (r * r);
					const double dfd =
					    (dchi * (psi * c3 - 1.0) + chi * dpsi * (c3 + psi * dc3)) / (r * r0) - fd * (dr / r + dr0 / r0);
					*velocityRate = initialRate->first * fd + R0 * dfd + initialRate->second * gd + V0 * dgd;
				}
			}

			return std::fabs(Rout.x) < kCostSentinel && std::fabs(Vout.x) < kCostSentinel;
		}
	}  // namespace detail

	StateVectors propagateKeplerState(const StateVectors& initial, double dtDays)
	{
		if (!initial.isValid())
			return {};

		Vec3D R;
		Vec3D V;
		if (!detail::keplerPropUV(initial.first, initial.second, Constants::GAUSS_K * dtDays, R, V)) {
			return {};
		}

		return StateVectors(R, V);
	}

	StateVectors propagateKeplerStateRate(const StateVectors& initial, double dtDays, const StateVectors& initialRate,
	                                      double dtRate)
	{
		Vec3D r, v, dr, dv;
		if (!detail::keplerPropUV(initial.first, initial.second, Constants::GAUSS_K * dtDays, r, v, &initialRate,
		                          Constants::GAUSS_K * dtRate, &dr, &dv))
			throw std::runtime_error("Kepler state derivative failed");
		return {dr, dv};
	}

}  // namespace Calculations
