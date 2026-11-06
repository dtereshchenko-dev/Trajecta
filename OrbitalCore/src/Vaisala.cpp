#include "Vaisala.h"
#include "orbital_core/Constants.h"
#include "logger/Logger.h"

#include <algorithm>
#include <cmath>

namespace Calculations {

	namespace {


		constexpr double kMinArcDays = 0.01;

		constexpr double kMinRhoAU = 0.01;

		constexpr double kCoarseGridStepAU = 0.1;

		constexpr double kFineGridHalfWindowAU = 2.0;
		constexpr double kFineGridStepAU = 0.001;

		constexpr int kDescendingRunForMin = 3;

		constexpr double kMinExitErrorRatio = 1.05;

		constexpr double kExpectedMotionFloor = 1.0e-10;

		constexpr double kCostSentinel = 1.0e30;

		[[nodiscard]] inline double solveRho(const Vec3D& L, const Vec3D& Re, double r) noexcept
		{
			const double b = 2.0 * L.dot(Re);
			const double c = Re.dot(Re) - r * r;
			const double disc = b * b - 4.0 * c;

			if (disc < 0.0) {
				return -1.0;
			}

			const double sq = std::sqrt(disc);
			const double r1 = (-b + sq) * 0.5;
			const double r2 = (-b - sq) * 0.5;

			if (r1 > kMinRhoAU && r2 > kMinRhoAU) {
				return std::min(r1, r2);
			}

			if (r1 > kMinRhoAU) {
				return r1;
			}

			if (r2 > kMinRhoAU) {
				return r2;
			}

			return -1.0;
		}

	}  // namespace

	StateVectors vaisala(const ObsNode& obs1, const ObsNode& obs2, double rMin, double rMax)
	{
		const double dt = obs2.epochJdTdb() - obs1.epochJdTdb();
		if (std::abs(dt) < kMinArcDays) {
			LOG_DEBUG << "Väisälä: observations too close in time";
			return {};
		}

		const Vec3D& L1 = obs1.vector;
		const Vec3D& L2 = obs2.vector;
		const Vec3D Re1 = obs1.obsVector * (-1.0);
		const Vec3D Re2 = obs2.obsVector * (-1.0);

		const double absDt = std::abs(dt);

		auto relError = [&](double r) -> double {
			const double rho1 = solveRho(L1, Re1, r);
			const double rho2 = solveRho(L2, Re2, r);

			if (rho1 < 0.0 || rho2 < 0.0) {
				return kCostSentinel;
			}

			const Vec3D r1 = L1 * rho1 + Re1;
			const Vec3D r2 = L2 * rho2 + Re2;

			const double expected = Constants::GAUSS_K / (r * std::sqrt(r)) * absDt;
			if (expected < kExpectedMotionFloor) {
				return kCostSentinel;
			}

			const double cosT = std::clamp(r1.dot(r2) / (r1.magnitude() * r2.magnitude()), -1.0, 1.0);
			return std::abs(std::acos(cosT) - expected) / expected;
		};

		double bestErr = kCostSentinel;
		double bestR = -1.0;
		int descCount = 0;

		const long long nCoarse = std::llround((rMax - rMin) / kCoarseGridStepAU) + 1;
		for (long long i = 0; i < nCoarse; ++i) {
			const double r = rMin + static_cast<double>(i) * kCoarseGridStepAU;
			const double err = relError(r);

			if (err < bestErr) {
				bestErr = err;
				bestR = r;
				++descCount;
			}

			else if (bestR > 0.0 && err > bestErr * kMinExitErrorRatio && descCount >= kDescendingRunForMin) {
				break;
			}
		}

		if (bestR < 0.0) {
			LOG_DEBUG << "Väisälä: no solution in " << rMin << "–" << rMax << " AU";
			return {};
		}

		const double lo = std::max(rMin, bestR - kFineGridHalfWindowAU);
		const double hi = std::min(rMax, bestR + kFineGridHalfWindowAU);

		double fineErr = kCostSentinel;
		const long long nFine = std::llround((hi - lo) / kFineGridStepAU) + 1;

		for (long long i = 0; i < nFine; ++i) {

			const double r = lo + static_cast<double>(i) * kFineGridStepAU;
			const double err = relError(r);

			if (err < fineErr) {
				fineErr = err;
				bestR = r;
			}
		}

		bestErr = fineErr;

		const double rho1 = solveRho(L1, Re1, bestR);
		const double rho2 = solveRho(L2, Re2, bestR);
		const Vec3D R1 = L1 * rho1 + Re1;
		const Vec3D R2 = L2 * rho2 + Re2;

		const double r1mag = R1.magnitude();
		const double tau = Constants::GAUSS_K * dt;
		const double u = 1.0 / (r1mag * r1mag * r1mag);
		const double f = 1.0 - 0.5 * u * tau * tau;
		const double g = tau - (1.0 / 6.0) * u * tau * tau * tau;
		const Vec3D V1 = (R2 - R1 * f) / g;

		LOG_DEBUG << "Väisälä r=" << bestR << " AU, relErr=" << bestErr << ", |R|=" << r1mag
		          << ", |V|=" << V1.magnitude();
		return StateVectors(R1, V1);
	}

}  // namespace Calculations
