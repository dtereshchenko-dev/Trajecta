#pragma once

#include "hazard/EncounterPlane.h"
#include "hazard/TargetPlaneTrace.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/hazard/CloseApproach.h"

#include <Eigen/Core>

#include <optional>
#include <span>
#include <vector>

class IEphemerisSource;
class NBodyPropagator;
class NonGravModel;
struct IntegrationReport;

namespace orbital_core::hazard {

	struct CloseApproachTargetPlaneDifferentialInput {
		const CloseApproachEvent* event{};
		Eigen::VectorXd parameters;
		FitParameterCovariance parameterCovariance;
		Eigen::VectorXd parameterStepPerSigma;
		const NonGravModel* nominalNonGrav{};
		const NBodyPropagator* propagator{};
		IEphemerisSource* eph{};
		CloseApproachConfig closeApproach;
		EncounterPlaneKind planeKind{EncounterPlaneKind::ScaledBPlane};
		double epochTdb{};
		double sigma{};
		IntegrationReport* integrationReport{};
	};

	struct CloseApproachTargetPlaneTrace {
		CloseApproachEvent event;
		TargetPlaneTracePoint trace;
		Eigen::MatrixXd dCoordsDParams;
		BPlaneBasis basis;
		double gamma{1.0};
		double vInfKmS{};
	};

	// Outside encounter monitoring, intermediate epochs split propagation.
	// Monitoring instead rebases at encounter exits. The requested event
	// remains the sole output point.
	// External events project the fixed-time STM onto r_rel.dot(v_rel) = 0
	// using relative acceleration, including the target's motion. A degenerate
	// minimum has no regular event differential and returns nullopt.
	// Interior events use a refined incoming sphere crossing and the derivative
	// of each trajectory's own two-body pericenter. Both event-time terms are
	// included without changing the propagator's fixed-time STM. The surface
	// force transition estimates the target's acceleration from ephemeris
	// velocities at MJD TDB +/- 0.001 day; no force model is modified.
	[[nodiscard]] ORBITAL_CORE_API std::optional<CloseApproachTargetPlaneTrace> makeCloseApproachTargetPlaneTracePoint(
	    const CloseApproachTargetPlaneDifferentialInput& input, std::span<const double> intermediateCheckpointJds);

	// Evaluates several encounters of the same VA in one STM propagation.
	// Results stay aligned with events; an invalid encounter yields nullopt.
	[[nodiscard]] ORBITAL_CORE_API std::vector<std::optional<CloseApproachTargetPlaneTrace>>
	makeCloseApproachTargetPlaneTracePoints(const CloseApproachTargetPlaneDifferentialInput& input,
	                                        std::span<const CloseApproachEvent> events);

}  // namespace orbital_core::hazard
