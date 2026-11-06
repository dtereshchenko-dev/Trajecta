#include "orbital_core/hazard/Moid.h"

#include "orbital_core/Constants.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Types.h"
#include "orbital_core/Vec3D.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace orbital_core::hazard {
	namespace {
		struct OrbitState {
			double a{};
			double e{};
			Vec3D ux{};
			Vec3D uy{};
		};

		OrbitState makeState(const Orbit& o) noexcept
		{
			const double Om = o.getLongitudeOfAscendingNode();
			const double w = o.getArgumentOfPericentre();
			const double i = o.getInclination();

			const double cO = std::cos(Om);
			const double sO = std::sin(Om);
			const double cw = std::cos(w);
			const double sw = std::sin(w);
			const double ci = std::cos(i);
			const double si = std::sin(i);

			OrbitState s;
			s.a = o.getSemiMajorAxis();
			s.e = o.getExcentricity();
			s.ux = Vec3D(cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si);
			s.uy = Vec3D(-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si);
			return s;
		}

		Vec3D position(const OrbitState& s, double nu) noexcept
		{
			const double p = s.a * (1.0 - s.e * s.e);
			const double cf = std::cos(nu);
			const double sf = std::sin(nu);
			const double r = p / (1.0 + s.e * cf);
			return s.ux * (r * cf) + s.uy * (r * sf);
		}

		Vec3D positionDeriv(const OrbitState& s, double nu) noexcept
		{
			const double p = s.a * (1.0 - s.e * s.e);
			const double cf = std::cos(nu);
			const double sf = std::sin(nu);
			const double denom = 1.0 + s.e * cf;
			const double r = p / denom;
			const double dr = p * s.e * sf / (denom * denom);
			return s.ux * (dr * cf - r * sf) + s.uy * (dr * sf + r * cf);
		}

		Vec3D positionSecond(const OrbitState& s, double nu) noexcept
		{
			const double p = s.a * (1.0 - s.e * s.e);
			const double cf = std::cos(nu);
			const double sf = std::sin(nu);
			const double denom = 1.0 + s.e * cf;
			const double r = p / denom;
			const double dr = p * s.e * sf / (denom * denom);
			const double d2r = p * s.e * (cf * denom + 2.0 * s.e * sf * sf) / (denom * denom * denom);
			return s.ux * (d2r * cf - 2.0 * dr * sf - r * cf) + s.uy * (d2r * sf + 2.0 * dr * cf - r * sf);
		}

		constexpr int kCoarseN = 256;
		constexpr int kNewtonMaxIter = 50;
		constexpr double kNewtonTol = 1.0e-12;
		constexpr double kNewtonStepLimit = 0.3;

	}  // namespace

	MoidResult computeMoid(const Orbit& A, const Orbit& B) noexcept
	{
		MoidResult r;

		if (A.getType() != OrbitType::Elliptical || B.getType() != OrbitType::Elliptical) {
			return r;
		}

		const auto sa = makeState(A);
		const auto sb = makeState(B);

		double bestNu1 = 0.0;
		double bestNu2 = 0.0;
		double bestSq = std::numeric_limits<double>::max();

		const double step = Constants::TWOPI / static_cast<double>(kCoarseN);

		for (int i = 0; i < kCoarseN; ++i) {
			const double nu1 = i * step;
			const Vec3D R1 = position(sa, nu1);

			for (int j = 0; j < kCoarseN; ++j) {
				const double nu2 = j * step;
				const Vec3D D = R1 - position(sb, nu2);
				const double d2 = D.dot(D);

				if (d2 < bestSq) {
					bestSq = d2;
					bestNu1 = nu1;
					bestNu2 = nu2;
				}
			}
		}

		double nu1 = bestNu1;
		double nu2 = bestNu2;
		int iter = 0;
		bool converged = false;

		for (; iter < kNewtonMaxIter; ++iter) {
			const Vec3D D = position(sa, nu1) - position(sb, nu2);
			const Vec3D dR1 = positionDeriv(sa, nu1);
			const Vec3D dR2 = positionDeriv(sb, nu2);
			const double g1 = D.dot(dR1);
			const double g2 = -D.dot(dR2);

			if (std::fabs(g1) < kNewtonTol && std::fabs(g2) < kNewtonTol) {
				converged = true;
				break;
			}

			const Vec3D ddR1 = positionSecond(sa, nu1);
			const Vec3D ddR2 = positionSecond(sb, nu2);
			const double H11 = dR1.dot(dR1) + D.dot(ddR1);
			const double H22 = dR2.dot(dR2) - D.dot(ddR2);
			const double H12 = -dR1.dot(dR2);
			const double det = H11 * H22 - H12 * H12;

			if (std::fabs(det) < 1.0e-18)
				break;

			double dNu1 = -(H22 * g1 - H12 * g2) / det;
			double dNu2 = -(H11 * g2 - H12 * g1) / det;

			const double s = std::max(std::fabs(dNu1), std::fabs(dNu2));
			if (s > kNewtonStepLimit) {
				const double k = kNewtonStepLimit / s;
				dNu1 *= k;
				dNu2 *= k;
			}

			nu1 += dNu1;
			nu2 += dNu2;
		}

		r.trueAnomaly1Rad = Helpers::wrapTwoPi(nu1);
		r.trueAnomaly2Rad = Helpers::wrapTwoPi(nu2);
		r.moidAu = (position(sa, r.trueAnomaly1Rad) - position(sb, r.trueAnomaly2Rad)).magnitude();
		r.iterations = iter;
		r.converged = converged;
		return r;
	}

	namespace {
		double pointToOrbitDistance(const OrbitState& s, const Vec3D& P) noexcept
		{
			double bestNu = 0.0;
			double bestSq = std::numeric_limits<double>::max();
			const double step = Constants::TWOPI / static_cast<double>(kCoarseN);

			for (int i = 0; i < kCoarseN; ++i) {
				const double nu = i * step;
				const Vec3D D = position(s, nu) - P;
				const double d2 = D.dot(D);
				if (d2 < bestSq) {
					bestSq = d2;
					bestNu = nu;
				}
			}

			double nu = bestNu;
			for (int iter = 0; iter < kNewtonMaxIter; ++iter) {
				const Vec3D D = position(s, nu) - P;
				const Vec3D dR = positionDeriv(s, nu);
				const double g = D.dot(dR);
				if (std::fabs(g) < kNewtonTol)
					break;

				const Vec3D ddR = positionSecond(s, nu);
				const double H = dR.dot(dR) + D.dot(ddR);
				if (std::fabs(H) < 1.0e-18)
					break;

				double dNu = -g / H;
				if (std::fabs(dNu) > kNewtonStepLimit)
					dNu = (dNu > 0 ? 1 : -1) * kNewtonStepLimit;
				nu += dNu;
			}

			return (position(s, Helpers::wrapTwoPi(nu)) - P).magnitude();
		}

		double computeMoonMoid(const Orbit& target, double jdEpochTdb, IEphemerisSource& eph) noexcept
		{
			if (target.getType() != OrbitType::Elliptical)
				return std::numeric_limits<double>::quiet_NaN();

			const OrbitState sa = makeState(target);

			constexpr int kSamples = 730;
			constexpr double kYearDay = 365.25;
			const double dt = kYearDay / static_cast<double>(kSamples - 1);

			double best = std::numeric_limits<double>::max();
			for (int i = 0; i < kSamples; ++i) {
				EphemerisRecord r;
				try {
					r = eph.queryJD(EphemerisPlanet::Moon, jdEpochTdb + i * dt);
				} catch (...) {
					continue;
				}
				const double d = pointToOrbitDistance(sa, r.R);
				if (d < best)
					best = d;
			}
			return std::isfinite(best) ? best : std::numeric_limits<double>::quiet_NaN();
		}
	}  // namespace

	PlanetMoidTable computeMoidsToAllPlanets(const Orbit& target, double jdEpochTdb, IEphemerisSource* eph) noexcept
	{
		if (!eph)
			eph = getDefaultEphemerisSource();

		PlanetMoidTable table;
		if (!eph)
			return table;

		for (std::size_t i = 0; i < kPlanetCount; ++i) {
			const Planet planet = static_cast<Planet>(i);

			if (planet == Planet::Moon) {
				const double m = computeMoonMoid(target, jdEpochTdb, *eph);
				if (std::isfinite(m))
					table[Planet::Moon] = m;
				continue;
			}

			EphemerisRecord state;

			try {
				state = eph->queryJD(toEphemerisPlanet(planet), jdEpochTdb);
			}

			catch (...) {
				continue;
			}

			const double kInvGauss = 1.0 / Constants::GAUSS_K;
			auto planetOrbit = make_orbit_from_states_vectors({state.R, state.V * kInvGauss}, jdEpochTdb);

			if (!planetOrbit)
				continue;

			const auto result = computeMoid(target, *planetOrbit);
			if (result.isValid())
				table[planet] = result.moidAu;
		}

		return table;
	}
}  // namespace orbital_core::hazard
