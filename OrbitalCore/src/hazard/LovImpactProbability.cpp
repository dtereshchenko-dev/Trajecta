#include "hazard/LovImpactProbability.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace orbital_core::hazard {
	namespace {
		constexpr double kInvSqrtTwoPi = 0.398942280401432677939946059934381868;
		constexpr double kAbsTolerance = 1.0e-18;
		constexpr double kRelTolerance = 1.0e-8;
		constexpr int kMaxDepth = 20;

		[[nodiscard]] double normalInterval(double lo, double hi) noexcept
		{
			const double rootTwo = std::sqrt(2.0);
			if (lo >= 0.0)
				return 0.5 * (std::erfc(lo / rootTwo) - std::erfc(hi / rootTwo));
			if (hi <= 0.0)
				return 0.5 * (std::erfc(-hi / rootTwo) - std::erfc(-lo / rootTwo));
			return 0.5 * (std::erf(hi / rootTwo) - std::erf(lo / rootTwo));
		}

		[[nodiscard]] double simpson(double fa, double fm, double fb, double a, double b) noexcept
		{
			return (b - a) * (fa + 4.0 * fm + fb) / 6.0;
		}

		double adaptiveSimpson(const std::function<double(double)>& f, double a, double b, double fa, double fm,
		                       double fb, double whole, double tolerance, int depth)
		{
			const double m = 0.5 * (a + b);
			const double lm = 0.5 * (a + m);
			const double rm = 0.5 * (m + b);
			const double flm = f(lm);
			const double frm = f(rm);
			const double left = simpson(fa, flm, fm, a, m);
			const double right = simpson(fm, frm, fb, m, b);
			const double delta = left + right - whole;

			if (depth <= 0 || std::abs(delta) <= 15.0 * tolerance)
				return left + right + delta / 15.0;

			return adaptiveSimpson(f, a, m, fa, flm, fm, left, tolerance * 0.5, depth - 1) +
			       adaptiveSimpson(f, m, b, fm, frm, fb, right, tolerance * 0.5, depth - 1);
		}

		[[nodiscard]] double integrateAdaptive(const std::function<double(double)>& f, double a, double b)
		{
			if (!(a < b))
				return 0.0;

			const double m = 0.5 * (a + b);
			const double fa = f(a);
			const double fm = f(m);
			const double fb = f(b);
			const double whole = simpson(fa, fm, fb, a, b);
			const double tolerance = std::max(kAbsTolerance, kRelTolerance * std::abs(whole));
			return adaptiveSimpson(f, a, b, fa, fm, fb, whole, tolerance, kMaxDepth);
		}
	}  // namespace

	LovImpactProbabilityResult lovImpactProbability(const LovImpactProbabilityInput& input)
	{
		LovImpactProbabilityResult out;
		const TargetPlaneTracePoint& p = input.tracePoint;
		const double r = input.impactCrossSectionBodyRadii;

		if (!p.valid || !(r > 0.0) || !(p.stretchBodyRadii > 0.0) || !(p.widthBodyRadii > 0.0))
			return out;

		const double width = std::max(p.widthBodyRadii, 1.0e-4);
		const double stretch = p.stretchBodyRadii;
		const double c = std::cos(p.alphaRad);
		const double s = std::sin(p.alphaRad);
		const double csi = p.coordsKm.x() / p.bodyRadiusKm * c + p.coordsKm.y() / p.bodyRadiusKm * s;
		const double zeta = -p.coordsKm.x() / p.bodyRadiusKm * s + p.coordsKm.y() / p.bodyRadiusKm * c;
		const double b = p.coordsKm.norm() / p.bodyRadiusKm;

		out.sigmaImpact = std::max((std::abs(csi) - r) / width, 0.0);
		out.candidate = (b - 10.0 * width) < r;
		if (!out.candidate) {
			out.valid = true;
			return out;
		}

		const double xMin = std::max(-8.0, (-r - csi) / width);
		const double xMax = std::min(8.0, (r - csi) / width);
		if (!(xMin <= xMax)) {
			out.valid = true;
			return out;
		}

		const auto integrand = [&](double x) -> double {
			const double u = csi + width * x;
			const double cap2 = std::max(r * r - u * u, 0.0);
			const double cap = std::sqrt(cap2);
			const double normalX = kInvSqrtTwoPi * std::exp(-0.5 * x * x);
			const double lo = p.sigma + (-cap - zeta) / stretch;
			const double hi = p.sigma + (cap - zeta) / stretch;
			return normalX * normalInterval(lo, hi);
		};

		out.probability = std::clamp(integrateAdaptive(integrand, xMin, xMax), 0.0, 1.0);

		if (!std::isfinite(out.probability))
			out.probability = 0.0;

		out.valid = true;
		return out;
	}

}  // namespace orbital_core::hazard
