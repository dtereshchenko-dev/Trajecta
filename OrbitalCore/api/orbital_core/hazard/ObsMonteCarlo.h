#pragma once

#include "orbital_core/FitConfig.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/hazard/CloseApproach.h"
#include "orbital_core/hazard/Planet.h"

#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <string>
#include <vector>

class IEphemerisSource;
class NBodyPropagator;
class NonGravModel;
class ObsNode;

// Impact Monte Carlo: covariance draws or observation perturbation/refits.
// Each clone is integrated with surface stopping, never classified from a
// two-body pericenter estimate. Events must be forward from the clone epoch.
// Hits belong to the selected event's +/-1 day window. A previous collision
// with that target is not counted again at a later event. This is target-event
// sampling, not a global first-collision search over all Solar System bodies.
namespace orbital_core::hazard {

	struct ORBITAL_CORE_API ObsMcConfig {
		int numClones{1000};
		double noiseSigmaMult{1.0};
		std::uint64_t rngSeed{42};
		double altitudeKm{0.0};
	};

	struct ORBITAL_CORE_API ImpactCloneInfo {
		int cloneIndex{};
		int eventIndex{};
		double tcaJd{};
		Planet planet{Planet::Earth};
		double missKm{};


		Vec3D rBodyFixedKm{};
		double latitudeDeg{};
		double longitudeDeg{};
		bool bodyFixedValid{};

		double velocityKmS{};
		double epochMjdTdb{};
		double epochMjdUtc{};
		Vec3D vBodyFixedKmS{};
	};

	struct ORBITAL_CORE_API ObsMcResult {
		std::string earthFrameModel;
		double altitudeKm{};

		// ---- bookkeeping ----
		int numAttempted{};
		int numConverged{};
		int numFailedFit{};
		bool cancelled{};

		// ---- per-event statistics ----
		// hitsPerEvent[k] / numConverged = empirical impact probability for
		// events[k].  σ via binomial: √(p(1-p)/N).
		std::vector<int> hitsPerEvent;
		std::vector<double> minMissPerEventKm;
		std::vector<double> medianMissPerEventKm;

		int numAnyImpactClones{};
		double aggregateImpactProb{};
		double aggregateImpactSigma{};

		std::vector<ImpactCloneInfo> corridor;

		double centroidLatDeg{std::numeric_limits<double>::quiet_NaN()};
		double centroidLonDeg{std::numeric_limits<double>::quiet_NaN()};
		double centroidTcaJd{std::numeric_limits<double>::quiet_NaN()};

		// < 50% converged → arc too short for obs-MC to be meaningful.
		[[nodiscard]] bool isReliable() const noexcept { return numAttempted > 0 && numConverged * 2 >= numAttempted; }
	};

	// Run obs-MC over `events`; OpenMP-parallel.  `progressCb(done,total)`
	// returns false to cancel; pass `{}` for no callback.  `eph=nullptr`
	// resolves the global default.
	[[nodiscard]] ORBITAL_CORE_API ObsMcResult runObservationMonteCarlo(
	    const Orbit& nominalOrbit, std::span<const ObsNode> observations, const FitConfig& fitCfg,
	    const NBodyPropagator& propagator, std::span<const CloseApproachEvent> events, const ObsMcConfig& cfg = {},
	    IEphemerisSource* eph = nullptr, std::function<bool(int /*done*/, int /*total*/)> progressCb = {});

	// Fast corridor sampler: draw clones directly from the fit-epoch
	// cartesian covariance and propagate them non-linearly to the selected
	// close-approach events. This is the UI Risk Corridor path; it does not
	// refit observations per clone. Supply the same final objective-scaled
	// fit-epoch COV as LOV. The sampler uses it as supplied (also permits an
	// explicitly chosen formal COV for diagnostics), with no hidden rescov.
	[[nodiscard]] ORBITAL_CORE_API ObsMcResult runStateCovarianceMonteCarlo(
	    const Orbit& nominalOrbit, const StateCov6x6& stateCovariance, const NBodyPropagator& propagator,
	    std::span<const CloseApproachEvent> events, const ObsMcConfig& cfg = {}, IEphemerisSource* eph = nullptr,
	    std::function<bool(int /*done*/, int /*total*/)> progressCb = {});

	[[nodiscard]] ORBITAL_CORE_API ObsMcResult runParameterCovarianceMonteCarlo(
	    const Orbit& nominalOrbit, const FitParameterCovariance& parameterCovariance,
	    const NonGravModel& nominalNonGrav, const NBodyPropagator& propagator,
	    std::span<const CloseApproachEvent> events, const ObsMcConfig& cfg = {}, IEphemerisSource* eph = nullptr,
	    std::function<bool(int /*done*/, int /*total*/)> progressCb = {});
}  // namespace orbital_core::hazard
