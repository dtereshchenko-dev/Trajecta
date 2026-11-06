#include "fit/FitEpochPolicy.h"
#include <memory>
#include <span>

#include "orbital_core/TimeScale.h"

#include "TwoBodyPropagation.h"
#include "NBodyPropagator.h"
#include "logger/Logger.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"

#include <cmath>

namespace orbital_core::fit {
	namespace {

		[[nodiscard]] NonGravModel makePropagationNonGrav(const FitConfig& cfg, const NonGravModel* warmStartNonGrav,
		                                                  bool disableNonGrav) noexcept
		{
			if (disableNonGrav || !cfg.nonGrav.enabled) {
				NonGravModel ng;
				ng.setEnabled(false);
				ng.setModel(NgModel::None);
				return ng;
			}

			if (warmStartNonGrav && warmStartNonGrav->enabled())
				return *warmStartNonGrav;

			NonGravModel ng;
			ng.setEnabled(true);
			ng.setModel(static_cast<NgModel>(cfg.nonGrav.model));

			const auto requestedMask = static_cast<NgParameterMask>(cfg.nonGrav.solveMask);
			ng.setFitParameterMask(requestedMask == NgParameterMask::None ? NonGravModel::defaultFitParameterMaskFor(
			                                                                    static_cast<NgModel>(cfg.nonGrav.model))
			                                                              : requestedMask);

			if (cfg.nonGrav.useInitialValues) {
				ng.setParams(cfg.nonGrav.initialA1, cfg.nonGrav.initialA2, cfg.nonGrav.initialA3);
				ng.setDT(cfg.nonGrav.initialDT);
			}

			return ng;
		}

		[[nodiscard]] NBodyPropagator makePropagationContext(const FitConfig& cfg, const NonGravModel& ng)
		{
			NBodyPropagator prop;
			prop.setStepSize(cfg.nbody.stepSize);
			prop.setPlanetMask(cfg.nbody.enabled ? cfg.nbody.planetMask : 0);
			prop.setUseSb441(cfg.nbody.useSb441);
			prop.setTargetNaifId(cfg.nbody.targetNaifId);
			prop.setInactiveInteriorMassPolicy(cfg.nbody.inactiveInteriorMassPolicy);
			prop.nonGrav() = ng;
			if (auto* ephemeris = getDefaultEphemerisSource())
				prop = prop.forEncounterMonitoring(*ephemeris);
			return prop;
		}

		[[nodiscard]] Orbit_p copyOrbit(const Orbit_p& orbit)
		{
			return orbit ? std::make_shared<Orbit>(*orbit) : nullptr;
		}

	}  // namespace

	double defaultOutputEpochJdTdb(std::span<const ObsNode> observations)
	{
		if (observations.empty())
			return 0.0;

		const double jdFirstTdb = observations.front().epochJdTdb();
		const double jdLastTdb = observations.back().epochJdTdb();

		double epoch = std::floor(jdLastTdb) + 0.5;

		if (epoch > jdLastTdb && epoch - jdLastTdb > jdFirstTdb - (epoch - 1.0))
			epoch -= 1.0;

		return epoch;
	}

	double arcMidpointEpochTdbMjd(std::span<const ObsNode> observations)
	{
		if (observations.empty())
			return 0.0;

		const double firstTdb = observations.front().epochMjdTdb();
		const double lastTdb = observations.back().epochMjdTdb();
		return 0.5 * (firstTdb + lastTdb);
	}

	double arcMidpointEpochTdb(std::span<const ObsNode> observations)
	{
		const double mjd = arcMidpointEpochTdbMjd(observations);
		return mjd > 0.0 ? orbital_core::time::jdFromMjd(mjd) : 0.0;
	}

	Orbit_p rebaseInitialOrbitToFitEpoch(const Orbit_p& initial, std::span<const ObsNode> observations,
	                                     const FitConfig& cfg, const NonGravModel* warmStartNonGrav,
	                                     bool disableNonGrav)
	{
		if (!initial || observations.empty())
			return copyOrbit(initial);

		const double targetEpochMjd = arcMidpointEpochTdbMjd(observations);
		if (targetEpochMjd <= 0.0 || !std::isfinite(targetEpochMjd))
			return copyOrbit(initial);

		const double sourceEpochMjd = initial->epochMjd();
		if (std::abs(targetEpochMjd - sourceEpochMjd) <= 1.0e-9)
			return copyOrbit(initial);

		const NonGravModel ng = makePropagationNonGrav(cfg, warmStartNonGrav, disableNonGrav);

		Orbit_p rebased;
		if (cfg.nbody.enabled || ng.enabled()) {
			NBodyPropagator prop = makePropagationContext(cfg, ng);
			const Vec3D vReal = initial->getV() * Constants::GAUSS_K;
			const auto states =
			    prop.propagateAllWithVelocityMjd(initial->getR(), vReal, sourceEpochMjd, {targetEpochMjd}, ng);
			if (!states.empty()) {
				const auto& state = states.front();
				rebased = make_orbit_from_states_vectors(StateVectors(state.R, state.V * (1.0 / Constants::GAUSS_K)),
				                                         orbital_core::time::jdFromMjd(targetEpochMjd));
			}
		} else {
			const StateVectors state = Calculations::propagateKeplerState(
			    StateVectors(initial->getR(), initial->getV()), targetEpochMjd - sourceEpochMjd);
			if (state.isValid())
				rebased = make_orbit_from_states_vectors(state, orbital_core::time::jdFromMjd(targetEpochMjd));
		}
		if (rebased)
			rebased->setEpochMjd(targetEpochMjd);

		if (!rebased) {
			LOG_INFO << "[FitEpoch] rebase failed; keeping initial epoch MJD " << sourceEpochMjd;
			return copyOrbit(initial);
		}

		LOG_INFO << "[FitEpoch] rebased initial state epoch MJD " << sourceEpochMjd << " -> " << targetEpochMjd;
		return rebased;
	}

}  // namespace orbital_core::fit
