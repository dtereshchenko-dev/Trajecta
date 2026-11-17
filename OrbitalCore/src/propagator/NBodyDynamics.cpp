#include "propagator/NBodyDynamics.h"
#include <cstdint>

#include "propagator/ForceAudit.h"

#include "orbital_core/TimeScale.h"

#include "orbital_core/Constants.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/forces/EihRelativity.h"
#include "orbital_core/forces/EarthJ2.h"
#include "orbital_core/forces/ForceField.h"
#include "orbital_core/forces/SunJ2.h"
#include "orbital_core/forces/SunGravity.h"
#include "orbital_core/nongrav/NonGravForceWrapper.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include "orbital_core/perturbations/AsteroidPerturbations.h"
#include "orbital_core/perturbations/PlanetaryPerturbations.h"

#include <array>
#include <cmath>

namespace forces = orbital_core::forces;

namespace {

	constexpr std::array<double, 8> kMeanOrbitRadiiAu{0.38709927, 0.72333566, 1.00000261,  1.52371034,
	                                                  5.20288799, 9.53667594, 19.18916464, 30.06992276};

	constexpr std::array<double, 8> kInteriorBodyGm{
	    forces::detail::GM_MERCURY, forces::detail::GM_VENUS,  forces::detail::GM_EARTH,  forces::detail::GM_MARS,
	    forces::detail::GM_JUPITER, forces::detail::GM_SATURN, forces::detail::GM_URANUS, forces::detail::GM_NEPTUNE};

	[[nodiscard]] double smoothInteriorMassScale(double heliocentricDistanceAu, std::uint16_t bodyMask) noexcept
	{
		constexpr double kBlendFraction = 0.2;
		double scale = 1.0;
		for (std::size_t index = 0; index < kMeanOrbitRadiiAu.size(); ++index) {
			const double radius = kMeanOrbitRadiiAu[index];
			if (heliocentricDistanceAu <= radius)
				break;
			if ((bodyMask & static_cast<std::uint16_t>(1u << index)) != 0)
				continue;

			const double gmRatio = kInteriorBodyGm[index] / Constants::GM_SUN;
			if (heliocentricDistanceAu > radius * (1.0 + kBlendFraction)) {
				scale += gmRatio;
			} else {
				scale += gmRatio * ((heliocentricDistanceAu / radius - 1.0) / kBlendFraction);
			}
		}
		return scale;
	}

}  // namespace

void NBodyDynamics::setPlanetEnabled(EphemerisPlanet planet, bool enabled)
{
	const auto mask = orbital_core::planetaryBodyMask(planet);
	if (!mask)
		return;
	if (enabled)
		m_bodyMask |= mask;
	else
		m_bodyMask &= static_cast<std::uint16_t>(~mask);
}

bool NBodyDynamics::planetEnabled(EphemerisPlanet planet) const
{
	const auto mask = orbital_core::planetaryBodyMask(planet);
	if (!(m_bodyMask & mask))
		return false;
	for (const auto& body : forces::detail::kPerturbers)
		if (body.planet == planet && body.naifId == m_targetNaifId)
			return false;
	return true;
}

void NBodyDynamics::setAllPlanetsEnabled(bool enabled) { m_bodyMask = enabled ? FULL_BODY_MASK : 0; }

std::uint16_t NBodyDynamics::bodyMask() const noexcept { return m_bodyMask; }

void NBodyDynamics::setBodyMask(std::uint16_t mask) noexcept { m_bodyMask = mask; }

void NBodyDynamics::setUseSb441(bool enabled) noexcept { m_useSb441 = enabled; }

bool NBodyDynamics::useSb441() const noexcept { return m_useSb441; }

void NBodyDynamics::setInactiveInteriorMassPolicy(orbital_core::InactiveInteriorMassPolicy policy) noexcept
{
	m_inactiveInteriorMassPolicy = policy;
}

orbital_core::InactiveInteriorMassPolicy NBodyDynamics::inactiveInteriorMassPolicy() const noexcept
{
	return m_inactiveInteriorMassPolicy;
}

double NBodyDynamics::gravitationalParameter(EphemerisPlanet planet, double fallback) const noexcept
{
	if (!m_source)
		m_source = dynamic_cast<DE440Source*>(getDefaultEphemerisSource());
	if (m_source && m_source->isAvailable()) {
		const double gm = m_source->gravitationalParameter(planet);
		if (gm > 0.0 && std::isfinite(gm))
			return gm;
	}
	return fallback;
}

double NBodyDynamics::centralNewtonianGm(double heliocentricDistanceAu) const noexcept
{
	const double sunGm = gravitationalParameter(EphemerisPlanet::Sun, Constants::GM_SUN);
	if (m_inactiveInteriorMassPolicy != orbital_core::InactiveInteriorMassPolicy::SmoothCentralMass) {
		return sunGm;
	}
	return sunGm * smoothInteriorMassScale(heliocentricDistanceAu, m_bodyMask);
}

Vec3D NBodyDynamics::accelerationAtJd(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng) const
{
	return accelerationAtMjd(r, v, orbital_core::time::mjdFromJd(jd), ng);
}

Eigen::Matrix3d NBodyDynamics::accelerationPartialRMjd(const Vec3D& r, const Vec3D& v, double mjd,
                                                       const NonGravModel& ng) const
{
	Vec3D unusedAcceleration;
	Eigen::Matrix3d dAdR;
	Eigen::Matrix3d unusedV;
	evaluate(r, v, mjd, 0.0, ng, true, unusedAcceleration, dAdR, unusedV);
	return dAdR;
}

Eigen::Matrix3d NBodyDynamics::accelerationPartialVMjd(const Vec3D& r, const Vec3D& v, double mjd,
                                                       const NonGravModel& ng) const
{
	Vec3D unusedAcceleration;
	Eigen::Matrix3d unusedR;
	Eigen::Matrix3d dAdV;
	evaluate(r, v, mjd, 0.0, ng, true, unusedAcceleration, unusedR, dAdV);
	return dAdV;
}

Eigen::Matrix3d NBodyDynamics::accelerationPartialR(const Vec3D& r, const Vec3D& v, double jd,
                                                    const NonGravModel& ng) const
{
	return accelerationPartialRMjd(r, v, orbital_core::time::mjdFromJd(jd), ng);
}

Eigen::Matrix3d NBodyDynamics::accelerationPartialV(const Vec3D& r, const Vec3D& v, double jd,
                                                    const NonGravModel& ng) const
{
	return accelerationPartialVMjd(r, v, orbital_core::time::mjdFromJd(jd), ng);
}

void NBodyDynamics::accelerationAndPartialsJd(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng,
                                              Vec3D& acceleration, Eigen::Matrix3d& dAdR, Eigen::Matrix3d& dAdV) const
{
	accelerationAndPartialsMjd(r, v, orbital_core::time::mjdFromJd(jd), ng, acceleration, dAdR, dAdV);
}

namespace {

	template <typename Pipeline>
	void evaluateAudited(const Pipeline& pipeline, forces::ForceContext& context, Vec3D& acceleration,
	                     Eigen::Matrix3d& dAdR, Eigen::Matrix3d& dAdV);

}

Vec3D NBodyDynamics::accelerationAtMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng) const
{
	Vec3D acceleration;
	Eigen::Matrix3d unusedR, unusedV;
	evaluate(r, v, mjd, 0.0, ng, false, acceleration, unusedR, unusedV);
	return acceleration;
}

void NBodyDynamics::accelerationAndPartialsMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng,
                                               Vec3D& acceleration, Eigen::Matrix3d& dAdR, Eigen::Matrix3d& dAdV) const
{
	evaluate(r, v, mjd, 0.0, ng, true, acceleration, dAdR, dAdV);
}

void NBodyDynamics::evaluate(const Vec3D& r, const Vec3D& v, double epochMjd, double offsetDays, const NonGravModel& ng,
                             bool computePartials, Vec3D& acceleration, Eigen::Matrix3d& dAdR,
                             Eigen::Matrix3d& dAdV) const
{
	acceleration = Vec3D(0.0, 0.0, 0.0);
	dAdR = Eigen::Matrix3d::Zero();
	dAdV = Eigen::Matrix3d::Zero();
	if (r.magnitude() < 1.0e-12)
		return;

	if (!m_source)
		m_source = dynamic_cast<DE440Source*>(getDefaultEphemerisSource());
	if (!m_asteroidSource)
		m_asteroidSource = getDefaultAsteroidSource();

	forces::ForceContext context{};
	context.R = r;
	context.V = v;
	context.setEpochMjd(epochMjd, offsetDays);
	context.centralNewtonianGm = centralNewtonianGm(r.magnitude());
	context.sunGm = gravitationalParameter(EphemerisPlanet::Sun, Constants::GM_SUN);
	for (int index = 0; index < 8; ++index) {
		const auto& body = forces::detail::kPerturbers[index];
		context.planetGms[index] = gravitationalParameter(body.planet, body.gm);
	}
	context.moonGm = gravitationalParameter(EphemerisPlanet::Moon, forces::detail::GM_MOON);
	context.earthGm = gravitationalParameter(EphemerisPlanet::Earth, forces::detail::GM_EARTH);
	context.plutoGm = gravitationalParameter(EphemerisPlanet::Pluto, forces::detail::GM_PLUTO);
	context.ephemerisSource = m_source;
	context.asteroidSource = m_asteroidSource;
	context.bodyMask = m_bodyMask;
	if (m_targetNaifId)
		for (const auto& body : forces::detail::kPerturbers)
			if (body.naifId == m_targetNaifId)
				context.bodyMask &= static_cast<std::uint16_t>(~orbital_core::planetaryBodyMask(body.planet));
	context.targetNaifId = m_targetNaifId;
	context.useSb441 = m_useSb441;
	context.forceModel = m_forceModel;
	context.computePartials = computePartials;
	context.nonGrav = &ng;

	static constexpr forces::ForceField<forces::SunGravity, forces::PlanetaryPerturbations, forces::EihRelativity,
	                                    forces::SunJ2, forces::AsteroidPerturbations, forces::EarthJ2,
	                                    forces::NonGravForceWrapper>
	    kForcePipeline;

	if (forces::audit::enabled()) {
		evaluateAudited(kForcePipeline, context, acceleration, dAdR, dAdV);
		return;
	}
	kForcePipeline.evaluate(context, acceleration, dAdR, dAdV);
}

namespace {

	template <typename Pipeline>
	void evaluateAudited(const Pipeline& pipeline, forces::ForceContext& context, Vec3D& acceleration,
	                     Eigen::Matrix3d& dAdR, Eigen::Matrix3d& dAdV)
	{
		namespace audit = forces::audit;

		const Vec3D rEarth = context.ephemerisSource && context.ephemerisSource->isAvailable()
		                         ? context.ephemerisSource->getPlanetPositionOffsetMjd(
		                               EphemerisPlanet::Earth, context.epochMjd, context.offsetDays)
		                         : Vec3D(0.0, 0.0, 0.0);
		const double geocentricAu = (context.R - rEarth).magnitude();
		const bool inEarthJ2Region = geocentricAu < Constants::EARTH_LOCAL_DYNAMICS_AU;

		if (!audit::claim(inEarthJ2Region)) {
			pipeline.evaluate(context, acceleration, dAdR, dAdV);
			return;
		}

		++audit::evaluationIndex();
		audit::write("state", [&](auto& out) {
			out << " mjdTdb=" << context.mjd << " jdTdb=" << context.jd << " x=" << context.R.x << " y=" << context.R.y
			    << " z=" << context.R.z << " vx=" << context.V.x << " vy=" << context.V.y << " vz=" << context.V.z
			    << " rGeoAu=" << geocentricAu << " earthJ2Region=" << (inEarthJ2Region ? 1 : 0)
			    << " partials=" << (context.computePartials ? 1 : 0);
		});

		Vec3D previous(0.0, 0.0, 0.0);
		pipeline.forEach([&](const auto& force) {
			force.evaluate(context, acceleration, dAdR, dAdV);
			const Vec3D differenced = acceleration - previous;
			previous = acceleration;

			Vec3D isolated(0.0, 0.0, 0.0);
			Eigen::Matrix3d scratchR = Eigen::Matrix3d::Zero();
			Eigen::Matrix3d scratchV = Eigen::Matrix3d::Zero();
			force.evaluate(context, isolated, scratchR, scratchV);

			audit::write("force", [&](auto& out) {
				out << " name=" << force.name() << " ax=" << isolated.x << " ay=" << isolated.y << " az=" << isolated.z
				    << " norm=" << isolated.magnitude() << " differencedNorm=" << differenced.magnitude();
			});
		});

		if (context.havePlanets) {
			static constexpr const char* kPerturberNames[8] = {"Mercury", "Venus",  "Earth",  "Mars",
			                                                   "Jupiter", "Saturn", "Uranus", "Neptune"};
			for (int index = 0; index < 8; ++index) {
				const Vec3D& p = context.planetPositions[index];
				audit::write("body", [&](auto& out) {
					out << " name=" << kPerturberNames[index] << " gm=" << context.planetGms[index] << " x=" << p.x
					    << " y=" << p.y << " z=" << p.z;
				});
			}
			audit::write("body", [&](auto& out) {
				out << " name=Moon gm=" << context.moonGm << " x=" << context.moonPosition.x
				    << " y=" << context.moonPosition.y << " z=" << context.moonPosition.z;
			});
		}

		audit::write("total", [&](auto& out) {
			out << " ax=" << acceleration.x << " ay=" << acceleration.y << " az=" << acceleration.z
			    << " norm=" << acceleration.magnitude();
		});
	}

}  // namespace

#ifdef _MSC_VER
#pragma float_control(precise, on, push)
#endif
Vec3D NBodyDynamics::accelerationAtOffsetMjd(const Vec3D& r, const Vec3D& v, double epochMjd, double offsetDays,
                                             const NonGravModel& ng) const
{
	Vec3D acceleration;
	Eigen::Matrix3d unusedR, unusedV;
	evaluate(r, v, epochMjd, offsetDays, ng, false, acceleration, unusedR, unusedV);
	return acceleration;
}

void NBodyDynamics::accelerationAndPartialsOffsetMjd(const Vec3D& r, const Vec3D& v, double epochMjd, double offsetDays,
                                                     const NonGravModel& ng, Vec3D& acceleration, Eigen::Matrix3d& dAdR,
                                                     Eigen::Matrix3d& dAdV) const
{
	evaluate(r, v, epochMjd, offsetDays, ng, true, acceleration, dAdR, dAdV);
}
#ifdef _MSC_VER
#pragma float_control(pop)
#endif
