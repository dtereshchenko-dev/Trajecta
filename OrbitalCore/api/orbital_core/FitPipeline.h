#pragma once

#include "orbital_core/FitConfig.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/FitRunControl.h"
#include "orbital_core/FitPipelineOptions.h"
#include "orbital_core/FitSession.h"
#include "orbital_core/ObservationSet.h"
#include "orbital_core/hazard/BoundaryCrossing.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <ranges>
#include <span>
#include <string>

namespace orbital_core {

	// Borrowed, JD-sorted observation set, prepared as requested by the caller.
	// The owner must remain alive and unchanged throughout the synchronous run.
	struct FitPipelineInput {
		const ObservationSet& observations;
		// Applied only to the final fit, after selecting/extending the seed arc.
		// JD TDB; absent selects the existing rounded output-epoch policy.
		std::optional<double> outputEpochJdTdb;
		// Manual exclusions do not seed IOD/arc extension. All manual choices
		// lock the final fit selection, while automatic selection starts fresh.
		std::span<const ObservationOverride> observationOverrides{};
	};

	struct InitialOrbitInfo {
		Orbit_p orbit;
		std::string method;
		double distanceAu{};
		double speedAuPerDay{};
	};

	struct ArcExtensionInfo {
		std::size_t observationCount{};
		std::uint16_t configuredMask{};
		std::uint16_t effectiveMask{};
		std::uint16_t automaticMask{};
		std::uint16_t blockedMask{};
	};

	struct NominalEarthEntry {
		double altitudeKm{};
		hazard::BoundaryCrossing crossing;
	};

	struct FitPipelineResult {
		FitResult fit;
		FitConfig effectiveConfig;
		InitialOrbitInfo initial;
		std::optional<ArcExtensionInfo> extension;
		std::size_t fitInputObservations{};
		// Owning snapshot of both input pools and their preparation/provenance.
		// The final fit's keptObsIndices always address observations.all.
		// Populated on completion only; cancellation leaves this snapshot empty.
		ObservationSet observations;
		FitSelectionMode selection{FitSelectionMode::Configured};
		std::vector<ObservationOverride> observationOverrides;
		std::optional<NominalEarthEntry> nominalEarthEntry;
		std::string nominalEarthEntryError;

		// Derived from the display-epoch orbit and the submitted designation.
		// A failed fit cannot establish an orbital group, but a designation
		// can still identify a comet or an interstellar object.
		[[nodiscard]] OrbitalClass orbitalClass() const noexcept
		{
			if (!fit.orbit)
				return OrbitalClass::Unknown;
			const auto& designation = observations.classificationDesignation.empty() && observations.format != "ADES"
			                              ? observations.objectName
			                              : observations.classificationDesignation;
			const auto classification = orbitalClassification(*fit.orbit, designation);
			if (fit.converged || classification == OrbitalClass::Interstellar)
				return classification;
			if (classification == OrbitalClass::Comet || classification == OrbitalClass::HyperbolicComet ||
			    classification == OrbitalClass::ParabolicComet)
				return OrbitalClass::Comet;
			return OrbitalClass::Unknown;
		}

		[[nodiscard]] CometOrbitClasses cometClasses() const noexcept
		{
			if (!fit.orbit || !fit.converged)
				return {};
			const auto& designation = observations.classificationDesignation.empty() && observations.format != "ADES"
			                              ? observations.objectName
			                              : observations.classificationDesignation;
			return cometOrbitClasses(*fit.orbit, designation);
		}

		[[nodiscard]] bool usesCometPhotometry() const noexcept
		{
			const auto type = orbitalClass();
			if (type == OrbitalClass::Comet || type == OrbitalClass::HyperbolicComet ||
			    type == OrbitalClass::ParabolicComet)
				return true;
			if (type != OrbitalClass::Interstellar)
				return false;
			for (const auto& row : fit.residuals.data())
				if (row.nuclearMagnitude.has_value())
					return true;
			return false;
		}

		// Read-only view in fit/residual order; an empty index list means all rows.
		// Borrows this result, invalidated by its destruction, move or modification.
		// No observation copies; malformed indices throw on access.
		[[nodiscard]] auto fittedObservations() const&
		{
			const auto count = fit.keptObsIndices.empty() ? observations.all.size() : fit.keptObsIndices.size();
			return std::views::iota(std::size_t{0}, count) |
			       std::views::transform([this](std::size_t row) -> const ObsNode& {
				       const auto index = fit.keptObsIndices.empty() ? row : fit.keptObsIndices.at(row);
				       return observations.all.at(index);
			       });
		}

		auto fittedObservations() const&& = delete;
	};

	// Shared IOD -> arc extension -> FitSession orchestration for UI and CLI.
	// Loading, debiasing, export and hazard analysis are separate operations.
	// Use FitSession directly when an initial orbit is already available.
	class ORBITAL_CORE_API FitPipeline {
	   public:
		using InitialOrbitCallback = std::function<void(const InitialOrbitInfo&)>;

		explicit FitPipeline(FitConfig config, FitPipelineOptions options = {});

		// Cancellation returns fit.cancelled=true, without partial results/metadata.
		// IOD and individual arc-extension attempts finish before checking stop.
		// initialOrbitReady runs synchronously before extension; exceptions propagate.
		[[nodiscard]] FitPipelineResult run(FitPipelineInput input, const FitRunControl& control = {},
		                                    const InitialOrbitCallback& initialOrbitReady = {}) const;

	   private:
		FitConfig config_;
		FitPipelineOptions options_;
	};

}  // namespace orbital_core
