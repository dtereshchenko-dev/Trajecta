#pragma once

#include "IForceModel.h"
#include "orbital_core/Constants.h"
#include "orbital_core/perturbations/PlanetaryPerturbations.h"

#include <string_view>
#include <array>

namespace orbital_core::forces {

	class EihRelativity {
	   public:
		void evaluate(ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& dadv) const noexcept
		{
			if (ctx.forceModel.relativity == RelativityMode::Disabled)
				return;
			constexpr std::size_t kMaxBodies = 10;
			constexpr double kMinDistance = 1.0e-12;

			std::array<Vec3D, kMaxBodies> positions{};
			std::array<Vec3D, kMaxBodies> velocities{};
			std::array<double, kMaxBodies> gms{};
			std::size_t bodyCount = 0;

			if (ctx.havePlanets && ctx.havePlanetVelocities) {
				for (std::size_t index = 0; index < std::size(detail::kPerturbers); ++index) {
					const auto& body = detail::kPerturbers[index];
					if ((ctx.bodyMask & planetaryBodyMask(body.planet)) == 0)
						continue;

					positions[bodyCount] = detail::bodyPosition(ctx, index);
					velocities[bodyCount] = detail::bodyVelocity(ctx, index);
					gms[bodyCount] = detail::bodyGm(ctx, index);
					++bodyCount;
				}
			}

			const double asteroidRadius = ctx.R.magnitude();
			if (asteroidRadius < kMinDistance)
				return;

			const double sunGm = ctx.sunGm > 0.0 ? ctx.sunGm : Constants::GM_SUN;
			const Vec3D sunVelocity = ctx.ephemerisTerms.sunVelocity;

			const Vec3D asteroidVelocity = sunVelocity + ctx.V;
			const double sunVelocity2 = sunVelocity.dot(sunVelocity);
			const double asteroidVelocity2 = asteroidVelocity.dot(asteroidVelocity);

			std::array<double, kMaxBodies> heliocentricDistances{};
			std::array<double, kMaxBodies> inverseHeliocentricDistance3{};
			std::array<double, kMaxBodies> asteroidDistances{};
			std::array<double, kMaxBodies> inverseAsteroidDistance3{};
			std::array<Vec3D, kMaxBodies> barycentricVelocities{};
			std::array<Vec3D, kMaxBodies> newtonianAccelerations{};

			for (std::size_t body = 0; body < bodyCount; ++body) {
				heliocentricDistances[body] = positions[body].magnitude();
				asteroidDistances[body] = (ctx.R - positions[body]).magnitude();
				if (heliocentricDistances[body] < kMinDistance || asteroidDistances[body] < kMinDistance)
					return;
				inverseHeliocentricDistance3[body] =
				    1.0 / (heliocentricDistances[body] * heliocentricDistances[body] * heliocentricDistances[body]);
				inverseAsteroidDistance3[body] =
				    1.0 / (asteroidDistances[body] * asteroidDistances[body] * asteroidDistances[body]);
				barycentricVelocities[body] = sunVelocity + velocities[body];
			}

			Vec3D correction;
			Eigen::Matrix3d positionPartials = Eigen::Matrix3d::Zero();
			Eigen::Matrix3d velocityPartials = Eigen::Matrix3d::Zero();
			Vec3D potentialGradient;
			if (ctx.computePartials) {
				potentialGradient = ctx.R * (4.0 * sunGm / (asteroidRadius * asteroidRadius * asteroidRadius));
				for (std::size_t body = 0; body < bodyCount; ++body)
					potentialGradient += (ctx.R - positions[body]) * (4.0 * gms[body] * inverseAsteroidDistance3[body]);
			}

			// Planet terms: rows 1-3 of DESCANSO equation 4-26.
			for (std::size_t body = 0; body < bodyCount && ctx.forceModel.relativity == RelativityMode::Full; ++body) {
				double term = -4.0 * sunGm / asteroidRadius - sunGm / heliocentricDistances[body];
				for (std::size_t other = 0; other < bodyCount; ++other)
					term -= 4.0 * gms[other] / asteroidDistances[other];

				Vec3D bodyAcceleration = positions[body] * (-sunGm * inverseHeliocentricDistance3[body]);
				for (std::size_t other = 0; other < bodyCount; ++other) {
					if (other == body)
						continue;
					const Vec3D separation = positions[other] - positions[body];
					const double separationDistance = separation.magnitude();
					if (separationDistance < kMinDistance)
						return;
					term -= gms[other] / separationDistance;
					const double inverseSeparationDistance3 =
					    1.0 / (separationDistance * separationDistance * separationDistance);
					bodyAcceleration += separation * (gms[other] * inverseSeparationDistance3);
				}
				newtonianAccelerations[body] = bodyAcceleration;

				const Vec3D asteroidToBody = ctx.R - positions[body];
				const Vec3D& bodyVelocity = barycentricVelocities[body];
				const double projectedVelocity = asteroidToBody.dot(bodyVelocity) / asteroidDistances[body];
				term += asteroidVelocity2 + 2.0 * bodyVelocity.dot(bodyVelocity) -
				        4.0 * asteroidVelocity.dot(bodyVelocity) - 1.5 * projectedVelocity * projectedVelocity +
				        0.5 * (positions[body] - ctx.R).dot(bodyAcceleration);

				correction += (positions[body] - ctx.R) * (term * gms[body] * inverseAsteroidDistance3[body]);
				if (ctx.computePartials)
					accumulatePartials(asteroidToBody, asteroidVelocity, bodyVelocity, bodyAcceleration,
					                   potentialGradient, asteroidDistances[body],
					                   gms[body] * inverseAsteroidDistance3[body], term, positionPartials,
					                   velocityPartials);
			}

			// Planet terms: rows 4 and 5.
			for (std::size_t body = 0; body < bodyCount && ctx.forceModel.relativity == RelativityMode::Full; ++body) {
				const Vec3D asteroidToBody = ctx.R - positions[body];
				const Vec3D& bodyVelocity = barycentricVelocities[body];
				const double velocityProjection = asteroidToBody.dot(asteroidVelocity * 4.0 - bodyVelocity * 3.0);

				correction += (asteroidVelocity - bodyVelocity) *
				              (gms[body] * inverseAsteroidDistance3[body] * velocityProjection);
				correction += newtonianAccelerations[body] * (3.5 * gms[body] / asteroidDistances[body]);
			}

			// Sun terms.
			double sunTerm = -4.0 * sunGm / asteroidRadius;
			Vec3D sunAcceleration;
			for (std::size_t body = 0; body < bodyCount; ++body) {
				sunTerm -= 4.0 * gms[body] / asteroidDistances[body] + gms[body] / heliocentricDistances[body];
				sunAcceleration += positions[body] * (gms[body] * inverseHeliocentricDistance3[body]);
			}
			const double projectedSunVelocity = ctx.R.dot(sunVelocity) / asteroidRadius;
			sunTerm += asteroidVelocity2 + 2.0 * sunVelocity2 - 4.0 * asteroidVelocity.dot(sunVelocity) -
			           1.5 * projectedSunVelocity * projectedSunVelocity + 0.5 * (ctx.R * -1.0).dot(sunAcceleration);

			const double inverseAsteroidRadius3 = 1.0 / (asteroidRadius * asteroidRadius * asteroidRadius);
			correction -= ctx.R * (sunTerm * sunGm * inverseAsteroidRadius3);
			correction += (asteroidVelocity - sunVelocity) *
			              (sunGm * inverseAsteroidRadius3 * ctx.R.dot(asteroidVelocity * 4.0 - sunVelocity * 3.0));
			correction += sunAcceleration * (3.5 * sunGm / asteroidRadius);

			const double c2 = Constants::AU_PER_DAY * Constants::AU_PER_DAY;
			a += correction / c2;

			if (ctx.computePartials) {
				accumulatePartials(ctx.R, asteroidVelocity, sunVelocity, sunAcceleration, potentialGradient,
				                   asteroidRadius, sunGm * inverseAsteroidRadius3, sunTerm, positionPartials,
				                   velocityPartials);
				dadr += positionPartials / c2;
				dadv += velocityPartials / c2;
			}
		}

		[[nodiscard]] std::string_view name() const noexcept { return "EihRelativity"; }

	   private:
		// For each gravitating body the acceleration is
		// -f*T*r + f*(v-vj)*[r.(4v-3vj)] + 3.5*mu*aj/d, f=mu/d^3.
		// T includes the shared -4*sum(mu/d) potential and the projected velocity.
		static void accumulatePartials(const Vec3D& relativePosition, const Vec3D& velocity, const Vec3D& bodyVelocity,
		                               const Vec3D& bodyAcceleration, const Vec3D& potentialGradient, double distance,
		                               double factor, double term, Eigen::Matrix3d& dadr,
		                               Eigen::Matrix3d& dadv) noexcept
		{
			const Eigen::Vector3d r(relativePosition.x, relativePosition.y, relativePosition.z);
			const Eigen::Vector3d v(velocity.x, velocity.y, velocity.z);
			const Eigen::Vector3d vj(bodyVelocity.x, bodyVelocity.y, bodyVelocity.z);
			const Eigen::Vector3d aj(bodyAcceleration.x, bodyAcceleration.y, bodyAcceleration.z);
			const Eigen::Vector3d potential(potentialGradient.x, potentialGradient.y, potentialGradient.z);
			const Eigen::Vector3d w = v - vj;
			const Eigen::Vector3d velocityTerm = 4.0 * v - 3.0 * vj;
			const double inverseDistance2 = 1.0 / (distance * distance);
			const double projectedVelocity = r.dot(vj) / distance;
			const double projection = r.dot(velocityTerm);
			const Eigen::Vector3d termGradient =
			    potential - 0.5 * aj -
			    3.0 * projectedVelocity * (vj / distance - projectedVelocity * inverseDistance2 * r);
			dadr += factor * (-term * Eigen::Matrix3d::Identity() +
			                  r * (3.0 * term * inverseDistance2 * r - termGradient).transpose() +
			                  w * (velocityTerm - 3.0 * projection * inverseDistance2 * r).transpose() -
			                  3.5 * aj * r.transpose());

			dadv += factor * (-r * (2.0 * v - 4.0 * vj).transpose() + projection * Eigen::Matrix3d::Identity() +
			                  4.0 * w * r.transpose());
		}
	};

	static_assert(ForceModel<EihRelativity>);
}  // namespace orbital_core::forces
