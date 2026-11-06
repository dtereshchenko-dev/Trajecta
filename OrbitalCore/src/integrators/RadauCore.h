#pragma once

// Gauss-Radau spacing and the coefficient machinery Everhart's RA15 and
// Rein & Spiegel's IAS15 have in common: both use the same eight-node
// quadrature, the same b <-> g relations and the same way of carrying the
// polynomial from one step into the next.  They differ only in how a step is
// judged - the error estimate, the convergence test, whether a step can be
// rejected, and whether the state is accumulated with compensation.  Keeping
// the common part here means the two integrators cannot drift apart in the
// parts that must stay identical.
//
// Indices are one-based throughout, which is how Everhart wrote them and how
// the reference implementation is laid out; index 0 of each array is unused.

#include "IIntegrator.h"

#include <Eigen/Dense>
#include <array>

namespace orbital_core::integrators::radau {

	struct RadauConstants {
		std::array<double, 9> h{};
		std::array<double, 8> w{};
		std::array<double, 8> u{};
		std::array<double, 22> c{};
		std::array<double, 22> d{};
		std::array<double, 22> r{};
		double w1 = 0.5;
	};

	inline const RadauConstants& radauConstants()
	{
		static const RadauConstants values = [] {
			RadauConstants out;
			out.h[1] = 0.0;
			out.h[2] = 0.05626256053692215;
			out.h[3] = 0.18024069173689236;
			out.h[4] = 0.35262471711316964;
			out.h[5] = 0.54715362633055538;
			out.h[6] = 0.73421017721541053;
			out.h[7] = 0.88532094683909577;
			out.h[8] = 0.97752061356128750;

			for (int n = 2; n <= 8; ++n) {
				out.w[n - 1] = 1.0 / static_cast<double>(n + n * n);
				out.u[n - 1] = 1.0 / static_cast<double>(n);
			}

			constexpr std::array<int, 9> nw{0, 0, 0, 1, 3, 6, 10, 15, 21};
			out.c[1] = -out.h[2];
			out.d[1] = out.h[2];
			out.r[1] = 1.0 / (out.h[3] - out.h[2]);

			int la = 1;
			int lc = 1;
			for (int k = 3; k <= 7; ++k) {
				const int lb = la;
				la = lc + 1;
				lc = nw[k + 1];
				out.c[la] = -out.h[k] * out.c[lb];
				out.c[lc] = out.c[la - 1] - out.h[k];
				out.d[la] = out.h[2] * out.d[lb];
				out.d[lc] = -out.c[lc];
				out.r[la] = 1.0 / (out.h[k + 1] - out.h[2]);
				out.r[lc] = 1.0 / (out.h[k + 1] - out.h[k]);
				if (k == 3)
					continue;

				for (int l = 4; l <= k; ++l) {
					const int ld = la + l - 3;
					const int le = lb + l - 4;
					out.c[ld] = out.c[le] - out.h[k] * out.c[le + 1];
					out.d[ld] = out.d[le] + out.h[l - 1] * out.d[le + 1];
					out.r[ld] = 1.0 / (out.h[k + 1] - out.h[l - 1]);
				}
			}
			return out;
		}();
		return values;
	}

	template <int N>
	using RadauVector = Eigen::Matrix<double, N, 1>;

	template <int N>
	void clearRadauStorage(std::array<RadauVector<N>, 8>& b, std::array<RadauVector<N>, 8>& e,
	                       std::array<RadauVector<N>, 8>& bd)
	{
		for (int l = 1; l <= 7; ++l) {
			b[l].setZero();
			e[l].setZero();
			bd[l].setZero();
		}
	}

	template <int N>
	void computeRadauBeta(const std::array<RadauVector<N>, 8>& b, std::array<RadauVector<N>, 8>& g,
	                      const RadauConstants& rc)
	{
		for (int k = 0; k < N; ++k) {
			g[1](k) = b[1](k) + rc.d[1] * b[2](k) + rc.d[2] * b[3](k) + rc.d[4] * b[4](k) + rc.d[7] * b[5](k) +
			          rc.d[11] * b[6](k) + rc.d[16] * b[7](k);

			g[2](k) = b[2](k) + rc.d[3] * b[3](k) + rc.d[5] * b[4](k) + rc.d[8] * b[5](k) + rc.d[12] * b[6](k) +
			          rc.d[17] * b[7](k);

			g[3](k) = b[3](k) + rc.d[6] * b[4](k) + rc.d[9] * b[5](k) + rc.d[13] * b[6](k) + rc.d[18] * b[7](k);

			g[4](k) = b[4](k) + rc.d[10] * b[5](k) + rc.d[14] * b[6](k) + rc.d[19] * b[7](k);

			g[5](k) = b[5](k) + rc.d[15] * b[6](k) + rc.d[20] * b[7](k);
			g[6](k) = b[6](k) + rc.d[21] * b[7](k);
			g[7](k) = b[7](k);
		}
	}

	template <int N>
	void extrapolateRadauCoefficients(double q, std::array<RadauVector<N>, 8>& b, std::array<RadauVector<N>, 8>& e,
	                                  std::array<RadauVector<N>, 8>& bd, int completedSteps)
	{
		const double q2 = q * q;
		const double q3 = q2 * q;
		const double q4 = q3 * q;
		const double q5 = q4 * q;
		const double q6 = q5 * q;
		const double q7 = q6 * q;

		for (int k = 0; k < N; ++k) {
			if (completedSteps != 1) {
				for (int l = 1; l <= 7; ++l)
					bd[l](k) = b[l](k) - e[l](k);
			}

			e[1](k) = q * (b[1](k) + 2.0 * b[2](k) + 3.0 * b[3](k) + 4.0 * b[4](k) + 5.0 * b[5](k) + 6.0 * b[6](k) +
			               7.0 * b[7](k));

			e[2](k) = q2 * (b[2](k) + 3.0 * b[3](k) + 6.0 * b[4](k) + 10.0 * b[5](k) + 15.0 * b[6](k) + 21.0 * b[7](k));

			e[3](k) = q3 * (b[3](k) + 4.0 * b[4](k) + 10.0 * b[5](k) + 20.0 * b[6](k) + 35.0 * b[7](k));

			e[4](k) = q4 * (b[4](k) + 5.0 * b[5](k) + 15.0 * b[6](k) + 35.0 * b[7](k));

			e[5](k) = q5 * (b[5](k) + 6.0 * b[6](k) + 21.0 * b[7](k));
			e[6](k) = q6 * (b[6](k) + 7.0 * b[7](k));
			e[7](k) = q7 * b[7](k);

			for (int l = 1; l <= 7; ++l)
				b[l](k) = e[l](k) + bd[l](k);
		}
	}

	template <int N>
	IntegratedStateFull radauState(const RadauVector<N>& x, const RadauVector<N>& v)
	{
		IntegratedStateFull state;
		state.R = Vec3D(x(0), x(1), x(2));
		state.V = Vec3D(v(0), v(1), v(2));

		if constexpr (N >= 21) {
			for (int c = 0; c < 6; ++c)
				for (int r = 0; r < 3; ++r) {
					state.Phi(r, c) = x(3 + 3 * c + r);
					state.Phi(r + 3, c) = v(3 + 3 * c + r);
				}
		}

		if constexpr (N >= 33) {
			for (int c = 0; c < 4; ++c)
				for (int r = 0; r < 3; ++r) {
					state.Sens(r, c) = x(21 + 3 * c + r);
					state.Sens(r + 3, c) = v(21 + 3 * c + r);
				}
		}

		return state;
	}

	template <int N>
	void setRadauVariations(RadauVector<N>& x, RadauVector<N>& v, const IntegratedStateFull& state)
	{
		if constexpr (N >= 21) {
			for (int c = 0; c < 6; ++c)
				for (int r = 0; r < 3; ++r) {
					x(3 + 3 * c + r) = state.Phi(r, c);
					v(3 + 3 * c + r) = state.Phi(r + 3, c);
				}
		}

		if constexpr (N >= 33) {
			for (int c = 0; c < 4; ++c)
				for (int r = 0; r < 3; ++r) {
					x(21 + 3 * c + r) = state.Sens(r, c);
					v(21 + 3 * c + r) = state.Sens(r + 3, c);
				}
		}
	}

	template <int N>
	Vec3D positionFromRadau(const RadauVector<N>& x)
	{
		return Vec3D(x(0), x(1), x(2));
	}

	template <int N>
	Vec3D velocityFromRadau(const RadauVector<N>& v)
	{
		return Vec3D(v(0), v(1), v(2));
	}

	template <int N>
	void packState(const Vec3D& r, const Vec3D& v, RadauVector<N>& x, RadauVector<N>& velocity)
	{
		x.setZero();
		velocity.setZero();
		x(0) = r.x;
		x(1) = r.y;
		x(2) = r.z;
		velocity(0) = v.x;
		velocity(1) = v.y;
		velocity(2) = v.z;
	}


}  // namespace orbital_core::integrators::radau
