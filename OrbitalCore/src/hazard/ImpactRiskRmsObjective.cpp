#include "orbital_core/hazard/ImpactRiskRmsObjective.h"
#include <span>
#include <utility>

#include "OrbitOptimizer.h"
#include "fit/FitObservations.h"
#include "fit/OptimizerConfigFactory.h"
#include "logger/Logger.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <vector>

namespace orbital_core::hazard {

	namespace {
		[[nodiscard]] OrbitOptimizerConfig makeObjectiveConfig(const FitConfig& fitConfig,
		                                                       const NonGravModel& nominalNonGrav) noexcept
		{
			OrbitOptimizerConfig cfg = fit::makeOrbitOptimizerConfig(fitConfig);
			cfg.computeUncertainty = false;
			cfg.populateResiduals = false;
			cfg.logDiagnostics = false;
			cfg.logCacheStats = false;

			if (nominalNonGrav.enabled() && nominalNonGrav.numFitParams() > 0) {
				cfg.nonGrav.enabled = true;
				cfg.nonGrav.model = static_cast<orbital_core::NonGravModelKind>(nominalNonGrav.model());
				cfg.nonGrav.solveMask = static_cast<std::uint8_t>(nominalNonGrav.fitParameterMask());
				cfg.nonGrav.useInitialValues = true;
				cfg.nonGrav.initialA1 = nominalNonGrav.A1();
				cfg.nonGrav.initialA2 = nominalNonGrav.A2();
				cfg.nonGrav.initialA3 = nominalNonGrav.A3();
				cfg.nonGrav.initialDT = nominalNonGrav.DT();
			}

			return cfg;
		}
	}  // namespace

	ImpactRiskRmsScaling makeImpactRiskRmsScaling(Orbit_p fitEpochOrbit, std::span<const ObsNode> fitObservations,
	                                              std::span<const std::size_t> keptIndices, const FitConfig& fitConfig,
	                                              const NonGravModel& nominalNonGrav, double objectiveScale,
	                                              int numSamples)
	{
		ImpactRiskRmsScaling scaling;
		scaling.numSamples = std::max(0, numSamples);
		scaling.objectiveScale = (std::isfinite(objectiveScale) && objectiveScale > 0.0) ? objectiveScale : 1.0;

		const std::size_t keptCount = keptIndices.empty() ? fitObservations.size() : keptIndices.size();
		if (!fitEpochOrbit || keptCount < 3 || scaling.numSamples == 0)
			return scaling;

		scaling.gammaCache = std::make_shared<ImpactRiskRmsScaling::GammaNormalization>();

		try {
			auto orbitCopy = std::make_shared<Orbit>(*fitEpochOrbit);
			const auto config = makeObjectiveConfig(fitConfig, nominalNonGrav);
			auto observations = fit::prepareFitObservations(fitObservations, config);
			if (!keptIndices.empty()) {
				std::vector<std::uint8_t> included(observations.size(), 0);
				for (const std::size_t index : keptIndices)
					included.at(index) = 1;
				observations = fit::selectObservations(observations, included);
			}
			auto optimizer = std::make_shared<OrbitOptimizer>(orbitCopy, std::move(observations), config);

			auto optimizerMutex = std::make_shared<std::mutex>();

			scaling.objective = [optimizer, optimizerMutex](const Eigen::VectorXd& x) -> double {
				if (x.size() < 6 || !x.array().isFinite().all())
					return std::numeric_limits<double>::quiet_NaN();

				std::vector<double> params(static_cast<std::size_t>(x.size()));
				for (Eigen::Index i = 0; i < x.size(); ++i)
					params[static_cast<std::size_t>(i)] = x(i);

				const std::lock_guard<std::mutex> lock(*optimizerMutex);
				return optimizer->weightedObservationChiSquarePhysical(params);
			};

			scaling.linearization = [optimizer = std::move(optimizer), optimizerMutex](
			                            const Eigen::VectorXd& x) -> ImpactRiskRmsScaling::ObservationLinearization {
				ImpactRiskRmsScaling::ObservationLinearization out;
				if (x.size() < 6 || !x.array().isFinite().all())
					return out;

				std::vector<double> params(static_cast<std::size_t>(x.size()));
				for (Eigen::Index i = 0; i < x.size(); ++i)
					params[static_cast<std::size_t>(i)] = x(i);

				const std::lock_guard<std::mutex> lock(*optimizerMutex);
				const auto lin = optimizer->weightedObservationLinearizationPhysical(params);
				out.residuals = lin.residuals;
				out.jacobian = lin.jacobian;
				out.sumSquares = lin.sumSquares;
				out.valid = lin.valid;
				return out;
			};
		}

		catch (const std::exception& e) {
			LOG_INFO << "[ImpactRisk] RMS scaling objective unavailable: " << e.what();
			scaling.objective = {};
		}

		catch (...) {
			LOG_INFO << "[ImpactRisk] RMS scaling objective unavailable";
			scaling.objective = {};
		}

		return scaling;
	}

}  // namespace orbital_core::hazard
