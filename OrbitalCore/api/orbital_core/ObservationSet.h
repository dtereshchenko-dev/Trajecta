#pragma once

#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/ObservationDiagnostics.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace orbital_core {
	enum class GaiaRadiusSource { Unknown, Sbdb, FitEstimate, None };

	struct ObservationPreparationOptions {
		bool debias{false};
		AstrometricWeightProfile weightProfile{AstrometricWeightProfile::Vfcc17};
		std::filesystem::path runtimeDirectory;
	};

	struct ObservationPreparation {
		bool debiasRequested{false};
		AstrometricWeightProfile weightProfile{AstrometricWeightProfile::Vfcc17};
		std::filesystem::path biasTable;
		std::size_t correctedAll{};
		std::size_t correctedIod{};
	};

	// Owns the complete input and the independently deduplicated IOD pool.
	// Both are JD-sorted. Selection/clipping must not discard entries from all.
	// After debias, ObsNode::submittedCoordinates retains original angles.
	struct ObservationSet {
		std::vector<ObsNode> all;
		std::vector<ObsNode> iod;
		std::string objectName;
		std::string classificationDesignation;
		int targetNaifId{0};
		std::filesystem::path source;
		std::string format;
		std::optional<ObservationPreparation> preparation;
		std::vector<ObservationTimeGroup> coincidentObservations;
		std::vector<ObservationVisibilityDiagnostic> visibilityDiagnostics;
		std::optional<double> gaiaPhotocenterRadiusKm;
		GaiaRadiusSource gaiaRadiusSource{GaiaRadiusSource::Unknown};
	};

	struct GaiaFprObservation {
		double epochMjdUtc{};
		double raDeg{}, decDeg{};
		double sigmaRaArcsec{}, sigmaDecArcsec{}, correlation{};
		Vec3D observerOffsetEquatorialAu;
	};

	// Detects ADES/MPC from content, loads observer geometry and preserves the
	// existing observation ordering. No debias or fit is performed here.
	// Throws std::invalid_argument for an empty path, std::runtime_error for
	// read/parse/catalog failures or fewer than three usable IOD observations.
	[[nodiscard]] ORBITAL_CORE_API ObservationSet loadObservations(const std::filesystem::path& path);

	// Add an explicitly identified FPR dataset, retaining original MPC/ADES rows.
	// Radius is optional; its absence disables only the size-dependent correction.
	[[nodiscard]] ORBITAL_CORE_API ObservationSet
	withGaiaFprObservations(const ObservationSet& source, std::span<const GaiaFprObservation> gaia,
	                        std::optional<double> photocenterRadiusKm = {});

	// Prepare exactly once, after any diagnostic coordinate overrides.
	// Applies debias to both pools in the same order as the historical clients;
	// weighting itself remains in the fit. Repeated preparation is rejected.
	// Reload the source to change preparation options on the same observations.
	ORBITAL_CORE_API void prepareObservations(ObservationSet& observations,
	                                          const ObservationPreparationOptions& options = {});

	// Prepare an independent copy from submitted coordinates, even when the
	// source is already debiased. Preserves IDs, geometry and both full pools;
	// never reads the source file or changes the source snapshot.
	[[nodiscard]] ORBITAL_CORE_API ObservationSet
	prepareObservationCopy(const ObservationSet& source, const ObservationPreparationOptions& options = {});

}  // namespace orbital_core
