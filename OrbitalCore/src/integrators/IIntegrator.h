#pragma once

#include "orbital_core/TimeScale.h"
#include "orbital_core/Vec3D.h"
#include "propagator/VariationalEqs.h"

#include <Eigen/Dense>

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <stop_token>

class NonGravModel;

class IIntegrationDynamics {
   public:
	virtual ~IIntegrationDynamics() = default;

	// Parameter forcing must use the same physical frame as the force model,
	// including when an integration adapter translates the state origin.
	[[nodiscard]] virtual Eigen::Matrix<double, 3, 4> nonGravPartialsMjd(const Vec3D& r, const Vec3D& v, double mjd,
	                                                                     const NonGravModel& ng) const
	{
		return VariationalEqs::buildNongravSensForcing(ng, r, v, mjd);
	}

	[[nodiscard]] virtual Eigen::Matrix<double, 3, 4> nonGravPartialsOffsetMjd(const Vec3D& r, const Vec3D& v,
	                                                                           double epochMjd, double offsetDays,
	                                                                           const NonGravModel& ng) const
	{
		return nonGravPartialsMjd(r, v, epochMjd + offsetDays, ng);
	}

	// Time is MJD TDB throughout; see orbital_core/TimeScale.h for why the
	// integration stack must not be handed a JD.
	[[nodiscard]] virtual Vec3D accelerationAtMjd(const Vec3D& r, const Vec3D& v, double mjd,
	                                              const NonGravModel& ng) const = 0;

	virtual void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D& v, double mjd, const NonGravModel& ng,
	                                        Vec3D& acceleration, Eigen::Matrix3d& dAdR,
	                                        Eigen::Matrix3d& dAdV) const = 0;

	// Ephemeris-backed implementations override these to preserve both parts
	// through the reader. Single-epoch dynamics retain the default adapter.
	[[nodiscard]] virtual Vec3D accelerationAtOffsetMjd(const Vec3D& r, const Vec3D& v, double epochMjd,
	                                                    double offsetDays, const NonGravModel& ng) const
	{
		return accelerationAtMjd(r, v, epochMjd + offsetDays, ng);
	}

	virtual void accelerationAndPartialsOffsetMjd(const Vec3D& r, const Vec3D& v, double epochMjd, double offsetDays,
	                                              const NonGravModel& ng, Vec3D& acceleration, Eigen::Matrix3d& dAdR,
	                                              Eigen::Matrix3d& dAdV) const
	{
		accelerationAndPartialsMjd(r, v, epochMjd + offsetDays, ng, acceleration, dAdR, dAdV);
	}
};

struct IntegratedState {
	Vec3D R;
	Vec3D V;
};

struct IntegratedStateWithStm {
	Vec3D R;
	Vec3D V;
	Eigen::Matrix<double, 6, 6> Phi{Eigen::Matrix<double, 6, 6>::Identity()};
};

struct IntegratedStateFull {
	Vec3D R;
	Vec3D V;
	Eigen::Matrix<double, 6, 6> Phi{Eigen::Matrix<double, 6, 6>::Identity()};
	Eigen::Matrix<double, 6, 4> Sens{Eigen::Matrix<double, 6, 4>::Zero()};
};

struct DenseIntegratedStep {
	double beginMjd{};
	double endMjd{};
	std::function<IntegratedStateFull(double)> evaluate;
};

struct IntegrationOptions {
	double initialStepDays{10.0};
	// A signed fixed step for the next accepted-step interval; nullopt
	// restores adaptive control. The callback is copied for each session.
	std::function<std::optional<double>(double, double, const IntegratedState&)> nextEncounterStep;
	std::function<bool(double, const IntegratedStateFull&)> variationalObserver;
	// Reference cloend/propin: restart the polynomial and local STM on exit.
	bool restartOnEncounterExit{false};
	// Prefer this observer when the accepted epoch must retain the small
	// offset separately. It supersedes variationalObserver when both are set.
	std::function<bool(double, double, const IntegratedStateFull&)> variationalObserverOffset;
	// Owns the accepted-step polynomial, including any propagated variations.
	std::function<void(DenseIntegratedStep)> denseStepObserver;
	// Preserve the force epoch for encounter control. When provided, this
	// callback takes step, base epoch and offset, superseding nextEncounterStep.
	std::function<std::optional<double>(double, double, double, const IntegratedState&)> nextEncounterStepOffset;
	// Checked before integration and between corrector stages, including rejected
	// attempts. Cancellation returns a report, never a partially successful batch.
	std::stop_token cancel;
};

struct AcceptedStateStep {
	double beginMjd{0.0};
	double endMjd{0.0};
	IntegratedState begin;
	IntegratedState end;
	// Optional dense output in the same frame as begin/end, at an MJD TDB
	// inside this accepted step. Borrowed from the active integrator session:
	// call only during the observer; do not retain this callback or the step.
	// An empty callback requires the caller's usual independent propagation.
	std::function<IntegratedState(double)> interpolateMjd;
};

using AcceptedStateStepObserver = std::function<bool(const AcceptedStateStep&)>;

struct IntegrationReport {
	IntegratedState finalState;
	double reachedMjd{0.0};
	std::size_t completedOutputs{0};
	std::size_t acceptedSteps{0};
	std::size_t rejectedSteps{0};
	std::size_t rhsEvaluations{0};
	std::size_t correctorIterations{0};
	std::size_t encounterSteps{0};
	std::size_t encounterRestarts{0};
	std::size_t unconvergedEncounterSteps{0};
	bool stoppedByObserver{false};
	bool cancelled{false};
};

// Integrators own the complete adaptive session. targetMjds must be monotonic
// and all targets must lie in one direction from epoch. Continuation state is
// preserved between requested outputs. Every time here is MJD TDB.
class IIntegrator {
   public:
	virtual ~IIntegrator() = default;

	virtual IntegrationReport integrateState(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
	                                         const IntegratedState& initial, double epochMjd,
	                                         std::span<const double> targetMjds, std::span<IntegratedState> outputs,
	                                         const IntegrationOptions& options,
	                                         const AcceptedStateStepObserver& observer = {}) const = 0;

	virtual IntegrationReport integrateStateWithStm(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
	                                                const IntegratedStateWithStm& initial, double epochMjd,
	                                                std::span<const double> targetMjds,
	                                                std::span<IntegratedStateWithStm> outputs,
	                                                const IntegrationOptions& options) const = 0;

	virtual IntegrationReport integrateStateWithStmAndSens(const IIntegrationDynamics& dynamics, const NonGravModel& ng,
	                                                       const IntegratedStateFull& initial, double epochMjd,
	                                                       std::span<const double> targetMjds,
	                                                       std::span<IntegratedStateFull> outputs,
	                                                       const IntegrationOptions& options) const = 0;

	[[nodiscard]] virtual const char* name() const noexcept = 0;
};

using IIntegrator_p = std::unique_ptr<IIntegrator>;
