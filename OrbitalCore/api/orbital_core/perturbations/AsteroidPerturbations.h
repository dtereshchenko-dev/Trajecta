#pragma once

#include "orbital_core/forces/IForceModel.h"
#include <string_view>
#include "orbital_core/forces/PointMassHelpers.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "orbital_core/Vec3D.h"

#include <Eigen/Core>
#include <array>

namespace orbital_core::forces {
	namespace detail {
		// Rounded Aegis Table 2 ratios, only a fallback when the active DE
		// header lacks MAxxxx (also used by the approximate Kepler-3 model).
		inline constexpr double GM_CERES = Constants::GM_SUN * 4.7191422767e-10;
		inline constexpr double GM_PALLAS = Constants::GM_SUN * 1.0297360324e-10;
		inline constexpr double GM_JUNO = Constants::GM_SUN * 1.4471670475e-11;
		inline constexpr double GM_VESTA = Constants::GM_SUN * 1.3026836726e-10;
		inline constexpr double GM_IRIS = Constants::GM_SUN * 8.5890389994e-12;
		inline constexpr double GM_HYGIEA = Constants::GM_SUN * 4.2385986149e-11;
		inline constexpr double GM_EUNOMIA = Constants::GM_SUN * 1.5243642468e-11;
		inline constexpr double GM_PSYCHE = Constants::GM_SUN * 1.1978215785e-11;
		inline constexpr double GM_EUPHROSYNE = Constants::GM_SUN * 8.1331596146e-12;
		inline constexpr double GM_EUROPA = Constants::GM_SUN * 2.0216913527e-11;
		inline constexpr double GM_CYBELE = Constants::GM_SUN * 7.0687100328e-12;
		inline constexpr double GM_SYLVIA = Constants::GM_SUN * 1.6337820878e-11;
		inline constexpr double GM_THISBE = Constants::GM_SUN * 8.9653065564e-12;
		inline constexpr double GM_CAMILLA = Constants::GM_SUN * 1.0878696848e-11;
		inline constexpr double GM_DAVIDA = Constants::GM_SUN * 2.9345275748e-11;
		inline constexpr double GM_INTERAMNIA = Constants::GM_SUN * 2.1327387534e-11;

		// SB441-N16 asteroid table for SPK-driven perturbations.  NAIF ID =
		// 2000000 + MPC number.  Mirrors DE441's SB441-N16 set (16 most-perturbing
		// main-belt asteroids).
		struct Sb441Asteroid {
			int naifId;
			double gm;
			const char* name;
		};

		inline constexpr std::array<Sb441Asteroid, 16> kSb441Asteroids = {{
		    {2000001, GM_CERES, "Ceres"},
		    {2000002, GM_PALLAS, "Pallas"},
		    {2000003, GM_JUNO, "Juno"},
		    {2000004, GM_VESTA, "Vesta"},
		    {2000007, GM_IRIS, "Iris"},
		    {2000010, GM_HYGIEA, "Hygiea"},
		    {2000015, GM_EUNOMIA, "Eunomia"},
		    {2000016, GM_PSYCHE, "Psyche"},
		    {2000031, GM_EUPHROSYNE, "Euphrosyne"},
		    {2000052, GM_EUROPA, "Europa"},
		    {2000065, GM_CYBELE, "Cybele"},
		    {2000087, GM_SYLVIA, "Sylvia"},
		    {2000088, GM_THISBE, "Thisbe"},
		    {2000107, GM_CAMILLA, "Camilla"},
		    {2000511, GM_DAVIDA, "Davida"},
		    {2000704, GM_INTERAMNIA, "Interamnia"},
		}};

		// Fixed Keplerian elements for major asteroids (J2000 ecliptic, epoch JD 2460000.5).
		// Used as Kepler-3 fallback when SB441 SPK is not available.
		struct AsteroidElements {
			double epoch;
			double a;  // AU
			double e;
			double i;   // rad
			double Om;  // rad
			double w;   // rad
			double M0;  // rad
			double gm;  // AU³/day²
			int naifId;
		};

		inline constexpr int kNumAsteroids = 3;

		inline constexpr AsteroidElements kAsteroids[kNumAsteroids] = {
		    // 1 Ceres  (JPL#48, epoch 2460000.5)
		    {2460000.5, 2.7671817, 0.0788175, 0.18476649, 1.40080388, 1.28230143, 0.30046978, GM_CERES, 2000001},
		    // 2 Pallas (JPL#73, epoch 2460000.5)
		    {2460000.5, 2.7696317, 0.2300843, 0.60959167, 3.01799088, 5.42561411, 6.24565095, GM_PALLAS, 2000002},
		    // 4 Vesta  (JPL#48, epoch 2460000.5)
		    {2460000.5, 2.3630382, 0.0887575, 0.12460359, 1.81090651, 2.64590449, 2.00944975, GM_VESTA, 2000004},
		};

		[[nodiscard]] inline Vec3D asteroidPosition(const AsteroidElements& el, double jd) noexcept
		{
			return Calculations::positionAt(el.a, el.e, el.i, el.Om, el.w, el.M0, el.epoch, jd, Constants::GM_SUN,
			                                OrbitType::Elliptical);
		}
	}  // namespace detail

	class AsteroidPerturbations {
	   public:
		void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& /*dadv*/) const noexcept
		{
			const bool useSpk = ctx.useSb441 && ctx.asteroidSource && ctx.asteroidSource->isAvailable();
			const double sunGm = ctx.sunGm > 0.0 ? ctx.sunGm : Constants::GM_SUN;

			if (useSpk) {
				for (const auto& ast : detail::kSb441Asteroids) {
					if (ast.naifId == ctx.targetNaifId)
						continue;
					Vec3D rAst = ctx.asteroidSource->getPositionOffsetMjd(ast.naifId, ctx.epochMjd, ctx.offsetDays);
					const double kernelGm =
					    ctx.ephemerisSource ? ctx.ephemerisSource->asteroidGravitationalParameter(ast.naifId) : 0.0;
					const double gm = kernelGm > 0.0 ? kernelGm : ast.gm * (sunGm / Constants::GM_SUN);
					addContribution(ctx.R, rAst, gm, ctx.computePartials, a, dadr);
				}
			}

			else {
				for (int idx = 0; idx < detail::kNumAsteroids; ++idx) {
					if (!(ctx.bodyMask & (1u << (8 + idx))))
						continue;

					const auto& ast = detail::kAsteroids[idx];
					if (ast.naifId == ctx.targetNaifId)
						continue;
					Vec3D rAst = detail::asteroidPosition(ast, ctx.jd);
					const double kernelGm =
					    ctx.ephemerisSource ? ctx.ephemerisSource->asteroidGravitationalParameter(ast.naifId) : 0.0;
					addContribution(ctx.R, rAst, kernelGm > 0.0 ? kernelGm : ast.gm, ctx.computePartials, a, dadr);
				}
			}
		}

		[[nodiscard]] std::string_view name() const noexcept { return "AsteroidPerturbations"; }

	   private:
		static inline void addContribution(const Vec3D& r, const Vec3D& rAst, double gm, bool computePartials, Vec3D& a,
		                                   Eigen::Matrix3d& dadr) noexcept
		{
			const double raMag = rAst.magnitude();
			if (raMag < 1e-12)
				return;

			Vec3D q = rAst - r;
			const double qMag = q.magnitude();
			if (qMag < 1e-12)
				return;

			const double q2 = q.dot(q);
			const double q3 = q2 * qMag;
			a = a + (q * (gm / q3));
			if (computePartials)
				accumulatePointMass_dR(q, gm, dadr);
		}
	};

	static_assert(ForceModel<AsteroidPerturbations>);
}  // namespace orbital_core::forces
