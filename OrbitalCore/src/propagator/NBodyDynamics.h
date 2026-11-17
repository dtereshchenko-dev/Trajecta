#pragma once

#include "integrators/IIntegrator.h"
#include "orbital_core/ForceModelPolicy.h"
#include "orbital_core/ForceModelSettings.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/PlanetaryBodyMask.h"

#include <cstdint>

class DE440Source;
class SpkReader;

static constexpr std::uint16_t ASTEROID_BODY_MASK = 0x0700;
static constexpr std::uint16_t MOON_BODY_MASK = orbital_core::planetaryBodyMask(EphemerisPlanet::Moon);

static constexpr std::uint16_t PLANET_BODY_MASK =
    orbital_core::planetaryBodyMask(EphemerisPlanet::Mercury) |
    orbital_core::planetaryBodyMask(EphemerisPlanet::Venus) | orbital_core::planetaryBodyMask(EphemerisPlanet::Earth) |
    orbital_core::planetaryBodyMask(EphemerisPlanet::Mars) | orbital_core::planetaryBodyMask(EphemerisPlanet::Jupiter) |
    orbital_core::planetaryBodyMask(EphemerisPlanet::Saturn) |
    orbital_core::planetaryBodyMask(EphemerisPlanet::Uranus) |
    orbital_core::planetaryBodyMask(EphemerisPlanet::Neptune) | orbital_core::planetaryBodyMask(EphemerisPlanet::Pluto);

static constexpr std::uint16_t FULL_BODY_MASK = PLANET_BODY_MASK | ASTEROID_BODY_MASK | MOON_BODY_MASK;

class NBodyDynamics final : public IIntegrationDynamics {
   public:
	[[nodiscard]] Vec3D accelerationAtMjd(const Vec3D& r, const Vec3D& v, double mjd,
	                                      const NonGravModel& ng) const override;

	void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng,
	                                Vec3D& acceleration, Eigen::Matrix3d& dAdR, Eigen::Matrix3d& dAdV) const override;

	[[nodiscard]] Vec3D accelerationAtOffsetMjd(const Vec3D& r, const Vec3D& v, double epochMjd, double offsetDays,
	                                            const NonGravModel& ng) const override;

	void accelerationAndPartialsOffsetMjd(const Vec3D& r, const Vec3D& v, double epochMjd, double offsetDays,
	                                      const NonGravModel& ng, Vec3D& acceleration, Eigen::Matrix3d& dAdR,
	                                      Eigen::Matrix3d& dAdV) const override;

	[[nodiscard]] Eigen::Matrix3d accelerationPartialRMjd(const Vec3D& r, const Vec3D& v, double mjd,
	                                                      const NonGravModel& ng) const;

	[[nodiscard]] Eigen::Matrix3d accelerationPartialVMjd(const Vec3D& r, const Vec3D& v, double mjd,
	                                                      const NonGravModel& ng) const;

	[[nodiscard]] Vec3D accelerationAtJd(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng) const;

	void accelerationAndPartialsJd(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng,
	                               Vec3D& acceleration, Eigen::Matrix3d& dAdR, Eigen::Matrix3d& dAdV) const;

	[[nodiscard]] Eigen::Matrix3d accelerationPartialR(const Vec3D& r, const Vec3D& v, double jd,
	                                                   const NonGravModel& ng) const;

	[[nodiscard]] Eigen::Matrix3d accelerationPartialV(const Vec3D& r, const Vec3D& v, double jd,
	                                                   const NonGravModel& ng) const;

	void setPlanetEnabled(EphemerisPlanet planet, bool enabled);

	[[nodiscard]] bool planetEnabled(EphemerisPlanet planet) const;
	void setAllPlanetsEnabled(bool enabled);

	[[nodiscard]] std::uint16_t bodyMask() const noexcept;
	void setBodyMask(std::uint16_t mask) noexcept;

	void setUseSb441(bool enabled) noexcept;
	[[nodiscard]] bool useSb441() const noexcept;

	void setForceModel(const orbital_core::ForceModelSettings& settings) noexcept { m_forceModel = settings; }

	[[nodiscard]] const orbital_core::ForceModelSettings& forceModel() const noexcept { return m_forceModel; }

	void setTargetNaifId(int id) noexcept { m_targetNaifId = id; }

	[[nodiscard]] int targetNaifId() const noexcept { return m_targetNaifId; }

	void setInactiveInteriorMassPolicy(orbital_core::InactiveInteriorMassPolicy policy) noexcept;

	[[nodiscard]] orbital_core::InactiveInteriorMassPolicy inactiveInteriorMassPolicy() const noexcept;

   private:
	void evaluate(const Vec3D& r, const Vec3D& v, double epochMjd, double offsetDays, const NonGravModel& ng,
	              bool computePartials, Vec3D& acceleration, Eigen::Matrix3d& dAdR, Eigen::Matrix3d& dAdV) const;

	[[nodiscard]] double gravitationalParameter(EphemerisPlanet planet, double fallback) const noexcept;
	[[nodiscard]] double centralNewtonianGm(double heliocentricDistanceAu) const noexcept;

	std::uint16_t m_bodyMask{FULL_BODY_MASK};
	bool m_useSb441{false};
	orbital_core::ForceModelSettings m_forceModel;
	int m_targetNaifId{0};

	orbital_core::InactiveInteriorMassPolicy m_inactiveInteriorMassPolicy{
	    orbital_core::InactiveInteriorMassPolicy::ExplicitOnly};
	mutable DE440Source* m_source{nullptr};
	mutable SpkReader* m_asteroidSource{nullptr};
};
