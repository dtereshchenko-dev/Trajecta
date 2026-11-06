#include "orbital_core/FitPipeline.h"
#include "orbital_core/FitSession.h"

#include "InitialOrbit.h"
#include "IODExtender.h"
#include "fit/FitEpochPolicy.h"
#include "fit/FitCancellation.h"
#include "fit/FitPipelineDetail.h"
#include "fit/FitObservations.h"
#include "fit/FitSessionRunner.h"
#include "fit/OptimizerConfigFactory.h"
#include "fit/ObservationSelection.h"
#include "logger/Logger.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <unordered_set>
#include <vector>

namespace orbital_core {
	namespace {

		OrbitOptimizerConfig makeExtensionConfig(const FitConfig& config, PerturberSelectionPolicy perturberPolicy)
		{
			OrbitOptimizerConfig result;
			result.useWeights = config.useWeights;
			result.obsWeighting = config.obsWeighting;
			result.gaiaPhotocenterRadiusKm = config.gaiaPhotocenterRadiusKm;
			result.nbody.enabled = config.nbody.enabled;
			result.nbody.planetMask = config.nbody.planetMask;
			result.nbody.inactiveInteriorMassPolicy = config.nbody.inactiveInteriorMassPolicy;
			result.perturberSelectionMode = perturberPolicy == PerturberSelectionPolicy::Configured
			                                    ? Calculations::PerturberSelectionMode::ManualExact
			                                    : Calculations::PerturberSelectionMode::AutoByArcLength;
			result.nbody.useSb441 = config.nbody.useSb441;
			result.nbody.targetNaifId = config.nbody.targetNaifId;
			result.nbody.stepSize = config.nbody.stepSize;
			result.nonGrav.enabled = false;
			result.weightProfile = config.weightProfile;
			return result;
		}

		class PipelineRun {
		   public:
			PipelineRun(FitConfig config, FitPipelineOptions options, FitPipelineInput input,
			            const FitRunControl& control, const FitPipeline::InitialOrbitCallback& initialOrbitReady)
			    : options_(options), input_(input), control_(control), initialOrbitReady_(initialOrbitReady)
			{
				output_.effectiveConfig = std::move(config);
				if (!output_.effectiveConfig.gaiaPhotocenterRadiusKm)
					output_.effectiveConfig.gaiaPhotocenterRadiusKm = input_.observations.gaiaPhotocenterRadiusKm;
				if (input_.observations.targetNaifId) {
					auto& target = output_.effectiveConfig.nbody.targetNaifId;
					if (target && target != input_.observations.targetNaifId)
						throw std::invalid_argument("Configured target conflicts with observation identity");
					target = input_.observations.targetNaifId;
				}
			}

			FitPipelineResult run()
			{
				fit::checkCancellation(control_.cancel);
				validateInput();
				const auto manual =
				    fit::resolveObservationSelection(input_.observations.all, {}, input_.observationOverrides);
				std::unordered_set<ObservationId> excluded;
				for (std::size_t i = 0; i < input_.observations.all.size(); ++i) {
					if (!manual.included[i]) {
						excluded.insert(input_.observations.all[i].id);
						continue;
					}
					if (input_.observations.all[i].radar)
						continue;
					opticalIndices_.push_back(i);
					optical_.push_back(input_.observations.all[i]);
				}
				if (input_.observationOverrides.empty()) {
					for (const auto& o : input_.observations.iod)
						if (!excluded.contains(o.id))
							iod_.push_back(o);
				} else {
					// Rebuild after manual selection: the original IOD representative of
					// a duplicate may be excluded while another copy remains selected.
					iod_ = optical_;
					const auto end = std::unique(iod_.begin(), iod_.end(), [](const ObsNode& a, const ObsNode& b) {
						return a.JD == b.JD && a.submittedRaRad() == b.submittedRaRad() &&
						       a.submittedDecRad() == b.submittedDecRad() && a.magString == b.magString &&
						       a.obsCode == b.obsCode;
					});
					iod_.erase(end, iod_.end());
				}
				if (optical_.size() < 3 || iod_.size() < 3)
					throw std::invalid_argument("Initial orbit requires at least 3 non-excluded optical observations");
				fit::reportProgress(control_, "Starting IOD");
				const auto iod = Calculations::bestInitialOrbit(iod_, optical_, output_.effectiveConfig.weightProfile);
				fit::checkCancellation(control_.cancel);
				if (!iod.state.isValid())
					throw std::runtime_error("All IOD methods (Gauss/Herget/Vaisala) failed");

				output_.initial = {make_orbit_from_states_vectors(iod.state, iod.epoch), iod.method ? iod.method : "",
				                   iod.state.first.magnitude(), iod.state.second.magnitude()};

				if (!output_.initial.orbit)
					throw std::runtime_error("Failed to build orbit from IOD");

				if (initialOrbitReady_)
					initialOrbitReady_(output_.initial);

				fit::checkCancellation(control_.cancel);

				auto seed = prepareSeed(iod);
				fit::checkCancellation(control_.cancel);

				auto& config = output_.effectiveConfig;
				const auto finalSelection =
				    fit::resolveObservationSelection(input_.observations.all, mask_, input_.observationOverrides);
				output_.fitInputObservations = std::ranges::count(finalSelection.included, std::uint8_t{1});

				fit::reportProgress(control_, "Optimizing " + std::to_string(output_.fitInputObservations) + " obs...");

				output_.fit =
				    fit::runPreparedFitSession(config,
				                               {.initialOrbit = seed,
				                                .observations = input_.observations.all,
				                                .initialIncluded = mask_,
				                                .outputEpochJdTdb = input_.outputEpochJdTdb.value_or(outputEpochJdTdb_),
				                                .observationOverrides = input_.observationOverrides},
				                               prepared_, control_);
				fit::checkCancellation(control_.cancel);

				output_.observations = input_.observations;
				output_.selection = FitSelectionMode::Configured;
				output_.observationOverrides.assign(input_.observationOverrides.begin(),
				                                    input_.observationOverrides.end());
				fit::checkCancellation(control_.cancel);

				return std::move(output_);
			}

		   private:
			void validateInput() const
			{
				if (input_.observations.all.size() < 3 || input_.observations.iod.size() < 3)
					throw std::invalid_argument("FitPipeline: need at least 3 observations in each pool");

				if (input_.outputEpochJdTdb &&
				    (!std::isfinite(*input_.outputEpochJdTdb) || *input_.outputEpochJdTdb <= 0.0))
					throw std::invalid_argument("FitPipeline: output epoch must be a finite positive JD TDB");
			}

			std::optional<Calculations::ExtensionResult> extendInitialOrbit(
			    const Calculations::InitialOrbitResult& iod) const
			{
				const auto extensionConfig =
				    makeExtensionConfig(output_.effectiveConfig, options_.perturberSelectionPolicy);
				fit::reportProgress(control_, "Extending IOD winner arc...");
				const auto subarc = fit::mapIodSubarc(iod_, optical_, {iod.subarcStart, iod.subarcEnd});
				if (!subarc)
					return std::nullopt;

				std::vector<obs::PreparedObservation> opticalPrepared;
				opticalPrepared.reserve(opticalIndices_.size());
				for (const auto index : opticalIndices_)
					opticalPrepared.push_back(prepared_.at(index));
				auto extension = Calculations::attemptExtensions(output_.initial.orbit, optical_, subarc->first,
				                                                 subarc->last, extensionConfig, opticalPrepared);
				if (extension.orbit)
					return extension;
				return std::nullopt;
			}

			void applyExtension(Calculations::ExtensionResult& extension)
			{
				auto& config = output_.effectiveConfig;
				const bool useSelection = config.sigmaClip.enabled && options_.observationSelectionPolicy ==
				                                                          ObservationSelectionPolicy::ArcExtension;
				const bool restricts = extension.keptObs.size() < optical_.size();
				const bool spansAll = Calculations::spansInputRange(extension, optical_.size());
				const bool keepMask = useSelection && (!spansAll || restricts);
				if (options_.perturberSelectionPolicy == PerturberSelectionPolicy::AutomaticFromArc) {
					config.nbody.planetMask = extension.effectivePerturberMask;
					config.nbody.enabled = config.nbody.enabled && (extension.effectivePerturberMask != 0);
				}
				outputEpochJdTdb_ = fit::defaultOutputEpochJdTdb(
				    (useSelection && !extension.keptObs.empty()) ? extension.keptObs : input_.observations.all);
				if (keepMask) {
					mask_.assign(input_.observations.all.size(), 1);
					for (std::size_t i = 0; i < opticalIndices_.size(); ++i)
						mask_[opticalIndices_[i]] = extension.includedMask.at(i);
				}
			}

			Orbit_p prepareSeed(const Calculations::InitialOrbitResult& iod)
			{
				auto& config = output_.effectiveConfig;
				config.nbody.inactiveInteriorMassPolicy =
				    options_.perturberSelectionPolicy == PerturberSelectionPolicy::AutomaticFromArc
				        ? InactiveInteriorMassPolicy::SmoothCentralMass
				        : InactiveInteriorMassPolicy::ExplicitOnly;
				outputEpochJdTdb_ = fit::defaultOutputEpochJdTdb(input_.observations.all);
				prepared_ = fit::prepareFitObservations(input_.observations.all, fit::makeOrbitOptimizerConfig(config));
				auto extension = extendInitialOrbit(iod);
				if (!extension)
					return output_.initial.orbit;

				output_.extension = ArcExtensionInfo{
				    extension->keptObs.size(), extension->configuredPerturberMask, extension->effectivePerturberMask,
				    extension->automaticPerturberMask, extension->blockedPerturberMask};
				LOG_INFO << "[ATT-EXT] pipeline seed: selected=" << extension->keptObs.size() << " arcIdx=["
				         << extension->arcFirst << ".." << extension->arcLast << "]"
				         << " spanIdx=[" << extension->firstIncluded << ".." << extension->lastIncluded
				         << "] finalRMS=" << extension->finalRMS << "\""
				         << " perturbers effective=0x" << std::hex << extension->effectivePerturberMask << " auto=0x"
				         << extension->automaticPerturberMask << " blocked=0x" << extension->blockedPerturberMask
				         << std::dec;
				applyExtension(*extension);
				return extension->orbit;
			}

			FitPipelineOptions options_;
			FitPipelineInput input_;
			const FitRunControl& control_;
			const FitPipeline::InitialOrbitCallback& initialOrbitReady_;
			FitPipelineResult output_;
			double outputEpochJdTdb_{};
			std::vector<std::uint8_t> mask_;
			std::vector<ObsNode> optical_;
			std::vector<ObsNode> iod_;
			std::vector<std::size_t> opticalIndices_;
			std::vector<obs::PreparedObservation> prepared_;
		};

	}  // namespace

	FitPipeline::FitPipeline(FitConfig config, FitPipelineOptions options)
	    : config_(std::move(config)), options_(options)
	{
	}

	FitPipelineResult FitPipeline::run(FitPipelineInput input, const FitRunControl& control,
	                                   const InitialOrbitCallback& initialOrbitReady) const
	try {
		return PipelineRun{config_, options_, input, control, initialOrbitReady}.run();
	} catch (const fit::FitCancelled&) {
		return {.fit = {.cancelled = true}};
	}

}  // namespace orbital_core
