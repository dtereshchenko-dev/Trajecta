#include "orbital_core/hazard/SolutionPropagation.h"
#include "orbital_core/SpkReader.h"
#include "NBodyPropagator.h"
#include <cmath>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace orbital_core::hazard {
	namespace {
		class WorkerLimit {
		   public:
			explicit WorkerLimit(unsigned workers)
			{
#ifdef _OPENMP
				m_previous = omp_get_max_threads();
				omp_set_num_threads(static_cast<int>(workers));
#else
				if (workers != 1)
					throw std::invalid_argument("OrbitalCore was built without OpenMP; set workers=1");
#endif
			}

			~WorkerLimit()
			{
#ifdef _OPENMP
				omp_set_num_threads(m_previous);
#endif
			}

		   private:
#ifdef _OPENMP
			int m_previous{};
#endif
		};
	}  // namespace

	SolutionPropagationResult propagateSolution(const OrbitSolution& solution, const SolutionPropagationConfig& config,
	                                            std::stop_token cancel, std::function<bool(int, int)> progress)
	{
		const auto& ca = config.encounter;
		if (config.workers == 0 || config.workers > 1024 || !config.forces.enabled ||
		    (config.integrator != IntegratorKind::RA15 && config.integrator != IntegratorKind::IAS15) ||
		    !std::isfinite(config.forces.stepSize) || config.forces.stepSize <= 0 || !std::isfinite(ca.horizonYears) ||
		    ca.horizonYears <= 0 || ca.lookbackYears != 0 || !std::isfinite(ca.sampleStepDay) ||
		    ca.sampleStepDay <= 0 || !std::isfinite(ca.encounterRadiusAu) || ca.encounterRadiusAu <= 0 ||
		    !std::isfinite(ca.thresholdAu) || ca.thresholdAu <= 0 || !std::isfinite(config.sampling.altitudeKm) ||
		    config.sampling.altitudeKm < 0 || config.sampling.numClones <= 0)
			throw std::invalid_argument("Invalid solution propagation settings");
		if (solution.covariance.rows() != 6 + solution.nonGrav.fitParameterCount() ||
		    solution.covariance.cols() != solution.covariance.rows() || !solution.covariance.allFinite())
			throw std::invalid_argument("Solution covariance and non-grav layout differ");
		SolutionPropagationResult result;
		result.mc.altitudeKm = config.sampling.altitudeKm;
		if (cancel.stop_requested()) {
			result.mc.cancelled = true;
			return result;
		}
		auto* eph = getDefaultEphemerisSource();
		if (!eph || !eph->isAvailable())
			throw std::runtime_error("Planetary ephemerides are unavailable");
		if (config.forces.useSb441 && (!getDefaultAsteroidSource() || !getDefaultAsteroidSource()->isAvailable()))
			throw std::runtime_error("SB441 was requested but is unavailable");
		WorkerLimit workers(config.workers);
		NBodyPropagator prop(config.integrator);
		prop.setStepSize(config.forces.stepSize);
		prop.setPlanetMask(config.forces.planetMask);
		prop.setUseSb441(config.forces.useSb441);
		prop.setTargetNaifId(config.forces.targetNaifId);
		prop.setInactiveInteriorMassPolicy(config.forces.inactiveInteriorMassPolicy);
		prop.nonGrav() = solution.nonGrav;
		prop = prop.forEncounterMonitoring(*eph, Constants::ENCOUNTER_STEP_POINTS, ca.encounterRadiusAu);
		result.events = findCloseApproaches(solution.nominal, prop, ca, eph);
		if (cancel.stop_requested()) {
			result.mc.cancelled = true;
			return result;
		}
		if (result.events.empty())
			return result;
		auto onProgress = [&](int done, int total) {
			return !cancel.stop_requested() && (!progress || progress(done, total));
		};
		if (solution.nonGrav.enabled()) {
			result.mc = runParameterCovarianceMonteCarlo(solution.nominal, solution.covariance, solution.nonGrav, prop,
			                                             result.events, config.sampling, eph, onProgress);
		} else {
			const StateCov6x6 cov = solution.covariance;
			result.mc = runStateCovarianceMonteCarlo(solution.nominal, cov, prop, result.events, config.sampling, eph,
			                                         onProgress);
		}
		return result;
	}

}  // namespace orbital_core::hazard
