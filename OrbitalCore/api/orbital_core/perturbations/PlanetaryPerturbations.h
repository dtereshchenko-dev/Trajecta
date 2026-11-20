#pragma once

#include "orbital_core/forces/IForceModel.h"
#include "orbital_core/forces/PointMassHelpers.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/PlanetaryBodyMask.h"
#include "orbital_core/DE440Source.h"

#include <Eigen/Core>
#include <cstdint>
#include <iterator>
#include <string_view>

namespace orbital_core::forces {

	namespace detail {
		// Exact DE440 kernel-header values in AU^3/day^2. Earth and Moon are
		// separate point masses, split from GMB with the DE440 EMRAT.
		inline constexpr double GM_MERCURY = 4.9125001948893182e-11;
		inline constexpr double GM_VENUS = 7.2434523326441187e-10;
		inline constexpr double GM_EARTH_PLUS_MOON = Constants::GM_EARTH_PLUS_MOON;
		inline constexpr double GM_EARTH = Constants::GM_EARTH_ALONE;
		inline constexpr double GM_MOON = 1.0931894624024351e-11;
		inline constexpr double GM_MARS = 9.5495488297258119e-11;
		inline constexpr double GM_JUPITER = 2.8253458252257918e-7;
		inline constexpr double GM_SATURN = 8.4597059933762903e-8;
		inline constexpr double GM_URANUS = 1.2920265649682399e-8;
		inline constexpr double GM_NEPTUNE = 1.5243573478851939e-8;
		inline constexpr double GM_PLUTO = 2.1750964648933581e-12;

		struct PerturbingBody {
			double gm;
			EphemerisPlanet planet;
			int naifId;
		};

		inline constexpr PerturbingBody kPerturbers[] = {
		    {GM_MERCURY, EphemerisPlanet::Mercury, 199}, {GM_VENUS, EphemerisPlanet::Venus, 299},
		    {GM_EARTH, EphemerisPlanet::Earth, 399},     {GM_MARS, EphemerisPlanet::Mars, 499},
		    {GM_JUPITER, EphemerisPlanet::Jupiter, 599}, {GM_SATURN, EphemerisPlanet::Saturn, 699},
		    {GM_URANUS, EphemerisPlanet::Uranus, 799},   {GM_NEPTUNE, EphemerisPlanet::Neptune, 899},
		    {GM_MOON, EphemerisPlanet::Moon, 301},       {GM_PLUTO, EphemerisPlanet::Pluto, 999},
		};

		[[nodiscard]] inline const Vec3D& bodyPosition(const ForceContext& ctx, std::size_t index) noexcept
		{
			if (index < 8)
				return ctx.planetPositions[index];
			return kPerturbers[index].planet == EphemerisPlanet::Moon ? ctx.moonPosition
			                                                          : ctx.ephemerisTerms.plutoPosition;
		}

		[[nodiscard]] inline const Vec3D& bodyVelocity(const ForceContext& ctx, std::size_t index) noexcept
		{
			if (index < 8)
				return ctx.planetVelocities[index];
			return kPerturbers[index].planet == EphemerisPlanet::Moon ? ctx.moonVelocity
			                                                          : ctx.ephemerisTerms.plutoVelocity;
		}

		[[nodiscard]] inline double bodyGm(const ForceContext& ctx, std::size_t index) noexcept
		{
			const auto& body = kPerturbers[index];
			const double gm = index < 8                              ? ctx.planetGms[index]
			                  : body.planet == EphemerisPlanet::Moon ? ctx.moonGm
			                                                         : ctx.plutoGm;
			return gm > 0.0 ? gm : body.gm;
		}
	}  // namespace detail

	class PlanetaryPerturbations {
	   public:
		void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& /*dadv*/) const
		{
			if (!ctx.ephemerisSource || !ctx.ephemerisSource->isAvailable())
				return;

			ctx.ephemerisSource->getAllPlanetStatesOffsetMjd(ctx.epochMjd, ctx.offsetDays, ctx.planetPositions,
			                                                 ctx.planetVelocities, &ctx.moonPosition, &ctx.moonVelocity,
			                                                 &ctx.ephemerisTerms);

			ctx.havePlanets = true;
			ctx.havePlanetVelocities = true;
			// Heliocentric origin follows DE independently of the direct-force mask.
			// Subtract it once; planets, Moon and asteroids contribute direct terms.
			a -= ctx.ephemerisTerms.sunAcceleration;

			for (std::size_t idx = 0; idx < std::size(detail::kPerturbers); ++idx) {
				const auto& body = detail::kPerturbers[idx];
				if (!(ctx.bodyMask & planetaryBodyMask(body.planet)))
					continue;
				const double bodyGm = detail::bodyGm(ctx, idx);

				const Vec3D& rPlanet = detail::bodyPosition(ctx, idx);
				const double rpMag = rPlanet.magnitude();

				if (rpMag < 1e-12)
					continue;

				Vec3D q = rPlanet - ctx.R;
				const double qMag = q.magnitude();
				if (qMag < 1e-12)
					continue;

				const double q2 = q.dot(q);
				const double q3 = q2 * qMag;
				a = a + (q * (bodyGm / q3));
				if (ctx.computePartials)
					accumulatePointMass_dR(q, bodyGm, dadr);
			}
		}

		[[nodiscard]] std::string_view name() const noexcept { return "PlanetaryPerturbations"; }
	};

	static_assert(ForceModel<PlanetaryPerturbations>);
}  // namespace orbital_core::forces
