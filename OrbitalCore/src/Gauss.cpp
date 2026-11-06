#include "Gauss.h"
#include <vector>
#include "orbital_core/Types.h"
#include "orbital_core/Orbit.h"
#include "logger/Logger.h"

#include <cmath>
#include <cstddef>

namespace Calculations {

	namespace {

		constexpr double kSlantRangeRejectAU = 10000.0;

		constexpr double kRhoConvTolAU = 1.0e-5;

		constexpr std::size_t kMaxRhoIter = 96;

		constexpr double kRootSearchStart = 0.0;
		constexpr double kRootSearchEnd = 10.0;
		constexpr double kRootSearchStep = 1.0e-3;

		constexpr double kDenomMinThreshold = 1.0e-14;

		constexpr double kMinHelioDistAU = 0.1;

		constexpr double kCostSentinel = 1.0e30;

		constexpr std::size_t kSmoothingItersCutoff = 4;

	}  // namespace

	namespace detail {

		std::vector<double> findInitialGuesses(double a, double b, double c, double start, double end, double step)
		{
			constexpr auto evalPoly = [](double aa, double bb, double cc, double x) noexcept -> double {
				const double x2 = x * x;
				const double x3 = x2 * x;
				const double x6 = x3 * x3;
				const double x8 = x6 * x2;
				return x8 + aa * x6 + bb * x3 + cc;
			};

			std::vector<double> guesses;

			if (step <= 0.0 || end <= start) {
				return guesses;
			}

			const long long nSteps = std::llround((end - start) / step);

			double x_prev = start;
			double f_prev = evalPoly(a, b, c, x_prev);

			for (long long i = 1; i <= nSteps; ++i) {

				const double x = start + static_cast<double>(i) * step;
				const double f = evalPoly(a, b, c, x);

				if (f_prev * f < 0.0 || f == 0.0) {
					guesses.push_back(0.5 * (x_prev + x));
				}

				x_prev = x;
				f_prev = f;
			}

			return guesses;
		}

	}  // namespace detail

	namespace {

		struct DMatrix {
			double d11, d12, d13;
			double d21, d22, d23;
			double d31, d32, d33;
			double d0;
		};

		struct CharacteristicPoly {
			double Ra;
			double Rb;
			double Rc;
		};

		struct IterationState {
			double r2 = 0.0;
			double rval = 0.0;
			Vec3D Rv;
			Vec3D Sv;
			double f1 = 0.0, g1 = 0.0;
			double f3 = 0.0, g3 = 0.0;
		};

		[[nodiscard]] DMatrix buildDMatrix(const ObsNode& obs1, const ObsNode& obs2, const ObsNode& obs3) noexcept
		{
			DMatrix m;

			m.d11 = -obs2.vector.dot(obs1.obsVector.cross(obs3.vector));
			m.d12 = -obs2.vector.dot(obs2.obsVector.cross(obs3.vector));
			m.d13 = -obs2.vector.dot(obs3.obsVector.cross(obs3.vector));

			m.d21 = obs1.obsVector.dot(obs3.vector.cross(obs1.vector));
			m.d22 = obs2.obsVector.dot(obs3.vector.cross(obs1.vector));
			m.d23 = obs3.obsVector.dot(obs3.vector.cross(obs1.vector));

			m.d31 = obs1.obsVector.dot(obs1.vector.cross(obs2.vector));
			m.d32 = obs2.obsVector.dot(obs1.vector.cross(obs2.vector));
			m.d33 = obs3.obsVector.dot(obs1.vector.cross(obs2.vector));

			m.d0 = obs1.vector.dot(obs2.vector.cross(obs3.vector));
			return m;
		}

		[[nodiscard]] CharacteristicPoly buildCharacteristicPoly(const DMatrix& m, const ObsNode& obs2, double tau1,
		                                                         double tau3) noexcept
		{
			const double tau = tau3 - tau1;

			const double A1 = tau3 / tau;
			const double A3 = -tau1 / tau;
			const double B1 = A1 * ((tau * tau) - (tau3 * tau3)) / 6.0;
			const double B3 = A3 * ((tau * tau) - (tau1 * tau1)) / 6.0;

			const double A = ((A1 * m.d21) - m.d22 + (A3 * m.d23)) / (-m.d0);
			const double B = (B1 * m.d21 + B3 * m.d23) / (-m.d0);

			const double E = -2.0 * obs2.obsVector.dot(obs2.vector);
			const double F = obs2.obsVector.dot(obs2.obsVector);

			CharacteristicPoly p;
			p.Ra = -(A * A + A * E + F);
			p.Rb = -Constants::MU_IN_UNIT * (2.0 * A * B + B * E);
			p.Rc = -Constants::MU_IN_UNIT * Constants::MU_IN_UNIT * B * B;
			return p;
		}

		[[nodiscard]] std::size_t filterPositiveRoots(std::vector<double>& roots) noexcept
		{
			std::size_t count = roots.size();

			for (std::size_t i = count; i-- > 0;) {
				if (roots[i] < 0.0) {
					roots[i] = roots[count - 1];
					--count;
				}
			}

			roots.resize(count);
			return count;
		}

		void applyHighOrderFGCorrections(double u2, double z2, double q, double tau1, double tau3, double& f1_new,
		                                 double& f3_new, double& g1_new, double& g3_new) noexcept
		{
			const double zs = z2 * z2;
			const double tau1_sq = tau1 * tau1;
			const double tau3_sq = tau3 * tau3;

			const double u2_z2_tau1_cubed = u2 * z2 * tau1_sq * tau1;
			const double u2_z2_tau3_cubed = u2 * z2 * tau3_sq * tau3;

			const double quartic = u2 * (3.0 * q - 15.0 * zs + u2) / 24.0;

			const double quint_f = u2 * (3.0 * q + 7.0 * zs - u2) * z2 / 8.0;
			const double quint_g = u2 * (9.0 * q - 45.0 * zs + u2) / 120.0;

			const double sixth_f =
			    u2 * (zs * (630.0 * q + 210.0 * u2 - 945.0 * zs) - u2 * (24.0 * q + u2) - 45.0 * q * q) / 720.0;
			const double sixth_g = u2 * z2 * (210.0 * zs - 90.0 * q - 15.0 * u2) / 360.0;

			const double tau1_fourth = tau1_sq * tau1_sq;
			const double tau3_fourth = tau3_sq * tau3_sq;

			f1_new += u2_z2_tau1_cubed / 2.0 + tau1_fourth * quartic;
			f3_new += u2_z2_tau3_cubed / 2.0 + tau3_fourth * quartic;
			g1_new += u2_z2_tau1_cubed * tau1 / 4.0;
			g3_new += u2_z2_tau3_cubed * tau3 / 4.0;

			f1_new += quint_f * tau1_fourth * tau1;
			f3_new += quint_f * tau3_fourth * tau3;
			g1_new += quint_g * tau1_fourth * tau1;
			g3_new += quint_g * tau3_fourth * tau3;

			f1_new += sixth_f * tau1_fourth * tau1_sq;
			f3_new += sixth_f * tau3_fourth * tau3_sq;
			g1_new += sixth_g * tau1_fourth * tau1_sq;
			g3_new += sixth_g * tau3_fourth * tau3_sq;
		}

		[[nodiscard]] IterationState iterateSingleRoot(double initialR2, const DMatrix& m, const ObsNode& obs1,
		                                               const ObsNode& obs2, const ObsNode& obs3, double tau1_init,
		                                               double tau3_init)
		{
			IterationState st;
			st.r2 = initialR2;
			st.rval = 0.0;

			LOG_DEBUG << "Initial rho: " << st.r2;

			double tau1 = tau1_init;
			double tau3 = tau3_init;
			bool keepIterating = true;

			for (std::size_t iter = 0; keepIterating && iter < kMaxRhoIter; ++iter) {

				const double r2_cubed = st.r2 * st.r2 * st.r2;

				if (r2_cubed <= 0.0) {
					LOG_DEBUG << "Gauss: non-positive r2 encountered";
					st.rval = 0.0;
					break;
				}

				const double u2 = Constants::MU_IN_UNIT / r2_cubed;
				const double tau1_sq = tau1 * tau1;
				const double tau3_sq = tau3 * tau3;

				double f1_new = 1.0 - u2 * tau1_sq / 2.0;
				double f3_new = 1.0 - u2 * tau3_sq / 2.0;
				double g1_new = tau1 * (1.0 - u2 * tau1_sq / 6.0);
				double g3_new = tau3 * (1.0 - u2 * tau3_sq / 6.0);

				if (iter > 0) {
					const double z2 = st.Rv.dot(st.Sv) / (st.r2 * st.r2);
					const double q = st.Sv.dot(st.Sv) / (st.r2 * st.r2) - u2;

					applyHighOrderFGCorrections(u2, z2, q, tau1, tau3, f1_new, f3_new, g1_new, g3_new);
				}

				if (iter > 0 && iter < kSmoothingItersCutoff) {
					st.f1 = 0.5 * (st.f1 + f1_new);
					st.f3 = 0.5 * (st.f3 + f3_new);
					st.g1 = 0.5 * (st.g1 + g1_new);
					st.g3 = 0.5 * (st.g3 + g3_new);
				} else {
					st.f1 = f1_new;
					st.f3 = f3_new;
					st.g1 = g1_new;
					st.g3 = g3_new;
				}

				const double denom = st.f1 * st.g3 - st.f3 * st.g1;
				if (std::abs(denom) < kDenomMinThreshold) {
					LOG_DEBUG << "Gauss: f·g denominator too small";
					st.rval = 0.0;
					break;
				}

				const double c1 = st.g3 / denom;
				const double c2 = -1.0;
				const double c3 = -st.g1 / denom;

				double p1 = (c1 * m.d11 + c2 * m.d12 + c3 * m.d13) / (c1 * m.d0);
				double p2 = (c1 * m.d21 + c2 * m.d22 + c3 * m.d23) / (c2 * m.d0);
				double p3 = (c1 * m.d31 + c2 * m.d32 + c3 * m.d33) / (c3 * m.d0);

				const double d1 = -st.f3 / denom;
				const double d3 = st.f1 / denom;

				LOG_DEBUG << "Iteration " << iter << ", r2: " << st.r2 << ", p1: " << p1 << ", p2: " << p2
				          << ", p3: " << p3;

				p1 = std::fabs(p1);
				p2 = std::fabs(p2);
				p3 = std::fabs(p3);

				if (p1 > kSlantRangeRejectAU || p2 > kSlantRangeRejectAU || p3 > kSlantRangeRejectAU) {
					st.rval = 0.0;
					keepIterating = false;
				} else {
					double r2_new_sq = 0.0;

					for (std::size_t i = 0; i < 3; ++i) {
						const double r1_comp = p1 * obs1.vector[i] - obs1.obsVector[i];
						const double r2_comp = p2 * obs2.vector[i] - obs2.obsVector[i];
						const double r3_comp = p3 * obs3.vector[i] - obs3.obsVector[i];

						st.Rv[i] = r2_comp;
						st.Sv[i] = d1 * r1_comp + d3 * r3_comp;

						r2_new_sq += r2_comp * r2_comp;
					}

					st.rval = p2;
					const double new_r2 = std::sqrt(r2_new_sq);

					const double rel = std::abs(new_r2 - st.r2) / std::max(1.0, std::abs(st.r2));
					if (std::abs(new_r2 - st.r2) < kRhoConvTolAU || rel < kRhoConvTolAU) {
						keepIterating = false;
					}

					st.r2 = new_r2;
				}

				tau1 = Constants::GAUSS_K *
				       (obs1.epochJdTdb() - obs2.epochJdTdb() - (p1 - st.rval) / Constants::AU_PER_DAY);
				tau3 = Constants::GAUSS_K *
				       (obs3.epochJdTdb() - obs2.epochJdTdb() - (p3 - st.rval) / Constants::AU_PER_DAY);
			}

			return st;
		}

	}  // namespace

	StateVectors gauss(const ObsNode& obs1, const ObsNode& obs2, const ObsNode& obs3)
	{
		const double tau1 = Constants::GAUSS_K * (obs1.epochJdTdb() - obs2.epochJdTdb());
		const double tau3 = Constants::GAUSS_K * (obs3.epochJdTdb() - obs2.epochJdTdb());

		const DMatrix m = buildDMatrix(obs1, obs2, obs3);
		const CharacteristicPoly poly = buildCharacteristicPoly(m, obs2, tau1, tau3);

		LOG_DEBUG << "Gauss poly coefficients: A=" << poly.Ra << " B=" << poly.Rb << " C=" << poly.Rc;

		std::vector<double> roots =
		    detail::findInitialGuesses(poly.Ra, poly.Rb, poly.Rc, kRootSearchStart, kRootSearchEnd, kRootSearchStep);

		LOG_DEBUG << "Gauss bracket scan: " << roots.size() << " sign-change root(s) in [" << kRootSearchStart << ", "
		          << kRootSearchEnd << "] AU";

		const std::size_t nRoots = filterPositiveRoots(roots);

		if (nRoots == 0) {
			LOG_DEBUG << "Gauss: no positive real roots survived (degenerate geometry)";
			return {};
		}

		double bestCost = kCostSentinel;
		double bestRval = 0.0;
		StateVectors bestSolution;

		for (std::size_t rootIdx = 0; rootIdx < nRoots; ++rootIdx) {
			const IterationState st = iterateSingleRoot(roots[rootIdx], m, obs1, obs2, obs3, tau1, tau3);

			if (st.rval <= 0.0 || st.Rv.magnitude() <= kMinHelioDistAU) {
				continue;
			}

			const auto testOrbit = make_orbit_from_states_vectors(StateVectors(st.Rv, st.Sv), obs2.epochJdTdb());

			if (!testOrbit || testOrbit->getType() == OrbitType::Undef) {
				continue;
			}

			const double cost = st.Rv.magnitude();
			LOG_DEBUG << "Root " << rootIdx << ": r2=" << roots[rootIdx] << " → r_helio=" << cost
			          << ", rval=" << st.rval;

			if (bestCost >= kCostSentinel || cost > bestCost) {
				bestCost = cost;
				bestSolution = StateVectors(st.Rv, st.Sv);
				bestRval = st.rval;
			}
		}

		if (bestRval <= 0.0) {
			return {};
		}

		return bestSolution;
	}

	std::vector<StateVectors> gaussAllRoots(const ObsNode& obs1, const ObsNode& obs2, const ObsNode& obs3)
	{
		const double tau1 = Constants::GAUSS_K * (obs1.epochJdTdb() - obs2.epochJdTdb());
		const double tau3 = Constants::GAUSS_K * (obs3.epochJdTdb() - obs2.epochJdTdb());

		const DMatrix m = buildDMatrix(obs1, obs2, obs3);
		const CharacteristicPoly poly = buildCharacteristicPoly(m, obs2, tau1, tau3);

		std::vector<double> roots =
		    detail::findInitialGuesses(poly.Ra, poly.Rb, poly.Rc, kRootSearchStart, kRootSearchEnd, kRootSearchStep);

		const std::size_t nRoots = filterPositiveRoots(roots);
		LOG_DEBUG << "gaussAllRoots: " << nRoots << " positive root(s) to iterate";

		std::vector<StateVectors> out;
		if (nRoots == 0)
			return out;
		out.reserve(nRoots);

		for (std::size_t rootIdx = 0; rootIdx < nRoots; ++rootIdx) {
			const IterationState st = iterateSingleRoot(roots[rootIdx], m, obs1, obs2, obs3, tau1, tau3);

			if (st.rval <= 0.0 || st.Rv.magnitude() <= kMinHelioDistAU)
				continue;

			const auto testOrbit = make_orbit_from_states_vectors(StateVectors(st.Rv, st.Sv), obs2.epochJdTdb());
			if (!testOrbit || testOrbit->getType() == OrbitType::Undef)
				continue;

			LOG_DEBUG << "gaussAllRoots root " << rootIdx << ": r2=" << roots[rootIdx]
			          << " → r_helio=" << st.Rv.magnitude();
			out.emplace_back(st.Rv, st.Sv);
		}
		return out;
	}

}  // namespace Calculations
