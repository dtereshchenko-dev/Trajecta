#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/Constants.h"
#include <memory>
#include <span>
#include <vector>
#include <cstdint>
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/ForceModelPolicy.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/IntegratorKind.h"
#include "integrators/IIntegrator.h"
#include "propagator/NBodyDynamics.h"
#include <Eigen/Dense>

class ORBITAL_CORE_API NBodyPropagator {
   public:
	struct StateVel {
		Vec3D R;
		Vec3D V;
	};

	explicit NBodyPropagator(IntegratorKind integrator = IntegratorKind::RA15);
	NBodyPropagator(const NBodyPropagator& other);
	NBodyPropagator& operator=(const NBodyPropagator& other);
	NBodyPropagator(NBodyPropagator&&) noexcept = default;
	NBodyPropagator& operator=(NBodyPropagator&&) noexcept = default;

	[[nodiscard]] NBodyPropagator forEncounterMonitoring(
	    IEphemerisSource& eph, int npoint = Constants::ENCOUNTER_STEP_POINTS,
	    double earthEncounterRadiusAu = Constants::EARTH_TARGET_PLANE_AU) const;

	[[nodiscard]] int encounterStepPointCount() const noexcept { return m_encounterPoints; }

	[[nodiscard]] double earthEncounterRadiusAu() const noexcept { return m_earthEncounterRadiusAu; }

	Vec3D propagatePosition(const Vec3D& R0, const Vec3D& V0, double epoch, double targetJD) const;

	std::vector<Vec3D> propagateAll(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                const std::vector<double>& targetJDs) const;

	Vec3D propagatePosition(const Vec3D& R0, const Vec3D& V0, double epoch, double targetJD,
	                        const NonGravModel& ng) const;

	std::vector<Vec3D> propagateAll(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                const std::vector<double>& targetJDs, const NonGravModel& ng) const;


	std::vector<StateVel> propagateAllWithVelocity(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                               const std::vector<double>& targetJDs) const;

	std::vector<StateVel> propagateAllWithVelocity(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                               const std::vector<double>& targetJDs, const NonGravModel& ng,
	                                               IntegrationReport* report = nullptr) const;

	struct PropagationWithSurfaceStop {
		Vec3D R{0.0, 0.0, 0.0};
		Vec3D V{0.0, 0.0, 0.0};
		double jd{0.0};
		bool hitSurface{false};
		double mjd{0.0};
	};

	enum class SurfaceShape { ReferenceEllipsoid, Sphere };

	PropagationWithSurfaceStop propagateUntilSurface(const Vec3D& R0, const Vec3D& V0, double epoch, double maxTargetJd,
	                                                 EphemerisPlanet planet, double rPlanetKm, IEphemerisSource& eph,
	                                                 SurfaceShape shape = SurfaceShape::ReferenceEllipsoid,
	                                                 double altitudeKm = 0.0) const;

	PropagationWithSurfaceStop propagateUntilSurface(const Vec3D& R0, const Vec3D& V0, double epoch, double maxTargetJd,
	                                                 EphemerisPlanet planet, double rPlanetKm, IEphemerisSource& eph,
	                                                 const NonGravModel& ng,
	                                                 SurfaceShape shape = SurfaceShape::ReferenceEllipsoid,
	                                                 double altitudeKm = 0.0) const;

	PropagationWithSurfaceStop propagateUntilSurface(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                                 std::span<const double> targetJds, std::vector<StateVel>& samples,
	                                                 EphemerisPlanet planet, double rPlanetKm, IEphemerisSource& eph,
	                                                 const NonGravModel& ng,
	                                                 SurfaceShape shape = SurfaceShape::ReferenceEllipsoid,
	                                                 double altitudeKm = 0.0) const;

	PropagationWithSurfaceStop propagateUntilSurfaceMjd(const Vec3D& R0, const Vec3D& V0, double epochMjd,
	                                                    std::span<const double> targetMjds,
	                                                    std::vector<StateVel>& samples, EphemerisPlanet planet,
	                                                    double rPlanetKm, IEphemerisSource& eph, const NonGravModel& ng,
	                                                    SurfaceShape shape, double altitudeKm,
	                                                    std::stop_token cancel = {}) const;

	Eigen::Matrix3d accelerationPartialR(const Vec3D& r, const Vec3D& v, double jd) const;
	Eigen::Matrix3d accelerationPartialV(const Vec3D& r, const Vec3D& v, double jd) const;

	Vec3D accelerationAt(const Vec3D& r, const Vec3D& v, double jd) const;
	Vec3D accelerationAt(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng) const;

	void accelerationAndPartials(const Vec3D& r, const Vec3D& v, double jd, const NonGravModel& ng, Vec3D& a_out,
	                             Eigen::Matrix3d& dAdR_out, Eigen::Matrix3d& dAdV_out) const;

	struct PropagationResult {
		Vec3D R;
		Vec3D V;
		Eigen::Matrix<double, 6, 6> Phi;
	};

	std::vector<PropagationResult> propagateAllWithStm(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                                   const std::vector<double>& targetJDs) const;

	std::vector<PropagationResult> propagateAllWithStm(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                                   const std::vector<double>& targetJDs, const NonGravModel& ng,
	                                                   IntegrationReport* report = nullptr) const;

	std::vector<PropagationResult> propagateAllWithStmChained(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                                          const std::vector<double>& targetJDs,
	                                                          const NonGravModel& ng,
	                                                          IntegrationReport* report = nullptr) const;

	struct PropagationResultFull {
		Vec3D R;
		Vec3D V;
		Eigen::Matrix<double, 6, 6> Phi;
		Eigen::Matrix<double, 6, 4> SensA;
	};

	struct PropagationWithStmSurfaceStop {
		PropagationResultFull state;
		double jd{};
		bool hitSurface{};
		double mjd{};
	};

	PropagationWithStmSurfaceStop propagateWithStmToSurface(const Vec3D& r0, const Vec3D& v0, double epoch,
	                                                        std::span<const double> checkpointJds,
	                                                        EphemerisPlanet planet, double radiusKm,
	                                                        IEphemerisSource& eph, const NonGravModel& ng,
	                                                        IntegrationReport* report = nullptr) const;

	Vec3D accelerationAtMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng) const;

	std::vector<PropagationResultFull> propagateAllWithStmAndSens(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                                              const std::vector<double>& targetJDs) const;

	std::vector<PropagationResultFull> propagateAllWithStmAndSens(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                                              const std::vector<double>& targetJDs,
	                                                              const NonGravModel& ng,
	                                                              IntegrationReport* report = nullptr) const;

	// Chained local STM with accumulated non-gravitational sensitivities:
	// S(t,t0)=Phi(t,ti)*S(ti,t0)+S(t,ti).
	std::vector<PropagationResultFull> propagateAllWithStmAndSensChained(const Vec3D& R0, const Vec3D& V0, double epoch,
	                                                                     const std::vector<double>& targetJDs,
	                                                                     const NonGravModel& ng,
	                                                                     IntegrationReport* report = nullptr) const;

	void setStepSize(double h) { m_stepSize = h; }

	std::vector<StateVel> propagateAllWithVelocityMjd(const Vec3D& R0, const Vec3D& V0, double epochMjd,
	                                                  const std::vector<double>& targetMjds, const NonGravModel& ng,
	                                                  std::stop_token cancel = {}) const;

	std::vector<PropagationResult> propagateAllWithStmMjd(const Vec3D& R0, const Vec3D& V0, double epochMjd,
	                                                      const std::vector<double>& targetMjds, const NonGravModel& ng,
	                                                      IntegrationReport* report = nullptr) const;

	std::vector<PropagationResultFull> propagateAllWithStmAndSensMjd(const Vec3D& R0, const Vec3D& V0, double epochMjd,
	                                                                 const std::vector<double>& targetMjds,
	                                                                 const NonGravModel& ng,
	                                                                 IntegrationReport* report = nullptr,
	                                                                 std::stop_token cancel = {}) const;

	double stepSize() const { return m_stepSize; }

	[[nodiscard]] const char* integratorName() const noexcept;
	[[nodiscard]] IntegratorKind integratorKind() const noexcept;
	void setIntegratorKind(IntegratorKind kind);

	void setPlanetEnabled(EphemerisPlanet planet, bool enabled);
	bool planetEnabled(EphemerisPlanet planet) const;
	void setAllPlanetsEnabled(bool enabled);

	uint16_t planetMask() const noexcept;
	void setPlanetMask(uint16_t mask);

	void setUseSb441(bool enabled);
	bool useSb441() const noexcept;
	void setTargetNaifId(int id) noexcept;
	[[nodiscard]] int targetNaifId() const noexcept;

	void setInactiveInteriorMassPolicy(orbital_core::InactiveInteriorMassPolicy policy) noexcept;
	[[nodiscard]] orbital_core::InactiveInteriorMassPolicy inactiveInteriorMassPolicy() const noexcept;

	NonGravModel& nonGrav();
	const NonGravModel& nonGrav() const;

	void setNonGravEnabled(bool enabled);
	bool nonGravEnabled() const;
	double A1() const noexcept;
	double A2() const noexcept;
	double A3() const noexcept;

	[[nodiscard]] double fixedAccuracyParameter() const noexcept;

   private:
	[[nodiscard]] IntegrationOptions integrationOptions() const;

	double m_stepSize = 10.0;

	IIntegrator_p m_integrator;
	IntegratorKind m_integratorKind{IntegratorKind::IAS15};
	NBodyDynamics m_dynamics;
	NonGravModel m_nonGrav;
	IEphemerisSource* m_encounterEphemeris{nullptr};
	int m_encounterPoints{0};
	double m_earthEncounterRadiusAu{Constants::EARTH_TARGET_PLANE_AU};
};

using NBodyPropagator_p = std::shared_ptr<NBodyPropagator>;
