#include "Herget.h"
#include "IodScoring.h"
#include "TwoBodyStumpff.h"
#include <vector>
#include "orbital_core/Constants.h"
#include "logger/Logger.h"

#include <algorithm>
#include <cmath>

namespace Calculations {

	namespace {
		constexpr double kCostSentinel = 1.0e30;

		constexpr int kLambertMaxIter = 150;
		constexpr double kLambertConvTol = 1.0e-9;
		constexpr double kPsiRampStep = 0.1;
		constexpr int kPsiRampMaxIter = 80;

		constexpr double kMinStateMagAU = 1.0e-12;
		constexpr double kMinDenomG = 1.0e-14;
		constexpr double kMinArcDays = 0.01;
		constexpr double kMinRhoAU = 0.01;
		constexpr double kMinFinalRhoAU = 0.05;
		constexpr double kMinFinalHelioAU = 1.0e-6;

		constexpr double kGridRhoMin = 0.5;
		constexpr double kGridRhoMax = 5.0;
		constexpr double kGridRhoStep = 0.25;

		constexpr double kFdStepFloor = 1.0e-4;
		constexpr double kFdStepRelScale = 1.0e-3;
		constexpr double kStepClampFrac = 0.3;
		constexpr int kLineSearchMaxIt = 6;
		constexpr double kHessianPosFloor = 1.0e-12;
		constexpr double kGradientFallback = 0.1;
		constexpr double kConvergenceTol = 1.0e-11;

		using detail::stumpffC;
		using detail::stumpffS;

		[[nodiscard]] bool lambertUV(const Vec3D& R1, const Vec3D& R2, double tau, Vec3D& V1out) noexcept
		{
			const double r1 = R1.magnitude();
			const double r2 = R2.magnitude();

			if (r1 < kMinStateMagAU || r2 < kMinStateMagAU || tau <= 0.0) {
				return false;
			}

			const double cosDnu = std::clamp(R1.dot(R2) / (r1 * r2), -1.0, 1.0);
			const Vec3D cr = R1.cross(R2);

			const double tm = (cr.z >= 0.0) ? 1.0 : -1.0;
			const double A = tm * std::sqrt(r1 * r2 * (1.0 + cosDnu));

			if (std::fabs(A) < kMinStateMagAU) {
				return false;
			}

			double psiLo = -4.0 * Constants::PI * Constants::PI;
			double psiHi = 4.0 * Constants::PI * Constants::PI;
			double psi = 0.0;
			double c2 = 0.5;
			double c3 = 1.0 / 6.0;
			double y = 0.0;

			for (int iter = 0; iter < kLambertMaxIter; ++iter) {

				c2 = stumpffC(psi);
				c3 = stumpffS(psi);
				y = r1 + r2 + A * (psi * c3 - 1.0) / std::sqrt(c2);

				if (A > 0.0 && y < 0.0) {

					double psiAdj = psi;
					for (int k = 0; k < kPsiRampMaxIter && y < 0.0; ++k) {

						psiAdj += kPsiRampStep;
						if (psiAdj > psiHi) {
							return false;
						}

						c2 = stumpffC(psiAdj);
						c3 = stumpffS(psiAdj);
						y = r1 + r2 + A * (psiAdj * c3 - 1.0) / std::sqrt(c2);
					}

					if (y < 0.0) {
						return false;
					}

					psiLo = psiAdj;
					psi = psiAdj;
				}

				const double x = std::sqrt(y / c2);
				const double tCurr = x * x * x * c3 + A * std::sqrt(y);

				if (std::fabs(tCurr - tau) < kLambertConvTol) {
					break;
				}

				if (tCurr < tau) {
					psiLo = psi;
				} else {
					psiHi = psi;
				}

				psi = 0.5 * (psiLo + psiHi);
			}

			const double f = 1.0 - y / r1;
			const double g = A * std::sqrt(y);

			if (std::fabs(g) < kMinDenomG) {
				return false;
			}

			V1out = (R2 - R1 * f) * (1.0 / g);
			return true;
		}

	}  // namespace

	StateVectors transferOrbitFromEndpointDistances(const ObsNode& obs1, const ObsNode& obs2, double rho1, double rho2)
	{
		if (rho1 <= kMinRhoAU || rho2 <= kMinRhoAU)
			return {};

		const double dtDays = obs2.epochJdTdb() - obs1.epochJdTdb();
		if (std::fabs(dtDays) < kMinArcDays)
			return {};

		const Vec3D R1 = obs1.vector * rho1 - obs1.obsVector;
		const Vec3D R2 = obs2.vector * rho2 - obs2.obsVector;
		Vec3D V1;
		if (!lambertUV(R1, R2, Constants::GAUSS_K * dtDays, V1))
			return {};

		return StateVectors(R1, V1);
	}

	namespace {
		inline double residualCost(const Vec3D& R1, const Vec3D& V1, double t1JD,
		                           const std::vector<ObsNode>& obs) noexcept
		{
			return keplerResidualCost(R1, V1, t1JD, obs);
		}
	}  // namespace

	StateVectors herget(const std::vector<ObsNode>& obs, double rho1Init, double rhoNInit, int maxIter)
	{
		if (obs.size() < 3) {
			LOG_DEBUG << "Herget: need ≥ 3 observations, got " << obs.size();
			return {};
		}

		const ObsNode& o1 = obs.front();
		const ObsNode& oN = obs.back();
		const double dtDays = oN.epochJdTdb() - o1.epochJdTdb();

		if (std::fabs(dtDays) < kMinArcDays) {
			LOG_DEBUG << "Herget: arc too short (Δt=" << dtDays << " d)";
			return {};
		}

		auto eval = [&](double rho1, double rhoN, Vec3D* R1out = nullptr, Vec3D* V1out = nullptr) -> double {
			if (rho1 <= kMinRhoAU || rhoN <= kMinRhoAU) {
				return kCostSentinel;
			}

			const StateVectors transfer = transferOrbitFromEndpointDistances(o1, oN, rho1, rhoN);
			if (!transfer.isValid()) {
				return kCostSentinel;
			}
			const Vec3D& R1 = transfer.first;
			const Vec3D& V1 = transfer.second;

			const double c = residualCost(R1, V1, o1.epochJdTdb(), obs);

			if (R1out != nullptr) {
				*R1out = R1;
			}
			if (V1out != nullptr) {
				*V1out = V1;
			}

			return c;
		};

		double rho1 = rho1Init;
		double rhoN = rhoNInit;

		if (rho1 <= 0.0 || rhoN <= 0.0) {

			double bestCost = kCostSentinel;
			const long long nCells = std::llround((kGridRhoMax - kGridRhoMin) / kGridRhoStep);

			for (long long i = 0; i <= nCells; ++i) {
				const double r1 = kGridRhoMin + static_cast<double>(i) * kGridRhoStep;

				for (long long j = 0; j <= nCells; ++j) {
					const double rN = kGridRhoMin + static_cast<double>(j) * kGridRhoStep;
					const double c = eval(r1, rN);

					if (c < bestCost) {
						bestCost = c;
						rho1 = r1;
						rhoN = rN;
					}
				}
			}

			if (bestCost >= kCostSentinel) {
				LOG_DEBUG << "Herget: grid search found no valid (ρ₁, ρₙ)";
				return {};
			}

			LOG_DEBUG << "Herget grid init: ρ₁=" << rho1 << ", ρₙ=" << rhoN << ", cost=" << bestCost;
		}

		double prevCost = kCostSentinel;
		Vec3D R1;
		Vec3D V1;

		for (int iter = 0; iter < maxIter; ++iter) {

			const double cost = eval(rho1, rhoN, &R1, &V1);
			if (cost >= kCostSentinel) {
				LOG_DEBUG << "Herget: non-finite cost at iter " << iter;
				return {};
			}

			const double h = std::max(kFdStepFloor, kFdStepRelScale * std::min(rho1, rhoN));
			const double c1p = eval(rho1 + h, rhoN);
			const double c1m = eval(rho1 - h, rhoN);
			const double cNp = eval(rho1, rhoN + h);
			const double cNm = eval(rho1, rhoN - h);

			if (c1p >= kCostSentinel || c1m >= kCostSentinel || cNp >= kCostSentinel || cNm >= kCostSentinel) {
				LOG_DEBUG << "Herget: FD neighbour non-finite at iter " << iter;
				break;
			}

			const double g1 = (c1p - c1m) / (2.0 * h);
			const double gN = (cNp - cNm) / (2.0 * h);
			const double h1 = (c1p - 2.0 * cost + c1m) / (h * h);
			const double hN = (cNp - 2.0 * cost + cNm) / (h * h);

			double d1 = (h1 > kHessianPosFloor) ? -g1 / h1 : -kGradientFallback * g1;
			double dN = (hN > kHessianPosFloor) ? -gN / hN : -kGradientFallback * gN;

			auto clampStep = [](double d, double rho) {
				const double lim = kStepClampFrac * rho;
				return std::max(-lim, std::min(lim, d));
			};
			d1 = clampStep(d1, rho1);
			dN = clampStep(dN, rhoN);

			double r1New = rho1 + d1;
			double rNNew = rhoN + dN;
			double newCost = eval(r1New, rNNew);

			for (int ls = 0; ls < kLineSearchMaxIt && !(newCost < cost); ++ls) {
				d1 *= 0.5;
				dN *= 0.5;
				r1New = rho1 + d1;
				rNNew = rhoN + dN;
				newCost = eval(r1New, rNNew);
			}

			if (!(newCost < cost)) {
				LOG_DEBUG << "Herget iter " << iter << ": line search stalled";
				break;
			}

			rho1 = std::max(kMinFinalRhoAU, r1New);
			rhoN = std::max(kMinFinalRhoAU, rNNew);

			LOG_DEBUG << "Herget iter " << iter << ": ρ₁=" << rho1 << ", ρₙ=" << rhoN << ", cost=" << newCost;

			if (std::fabs(prevCost - newCost) < kConvergenceTol) {
				break;
			}

			prevCost = newCost;
		}

		const double finalCost = eval(rho1, rhoN, &R1, &V1);
		if (finalCost >= kCostSentinel || R1.magnitude() < kMinFinalHelioAU) {
			LOG_DEBUG << "Herget: failed to produce valid state";
			return {};
		}

		LOG_DEBUG << "Herget done: ρ₁=" << rho1 << ", ρₙ=" << rhoN << ", cost=" << finalCost
		          << ", |R|=" << R1.magnitude() << ", |V|=" << V1.magnitude();

		return StateVectors(R1, V1);
	}

}  // namespace Calculations
