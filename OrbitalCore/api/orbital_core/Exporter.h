#pragma once

#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/FitResult.h"
#include "orbital_core/FitConfig.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/OrbitalCoreExport.h"

#include <filesystem>
#include <initializer_list>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace orbital_core {

	enum class ExportProduct {
		FitSummary,
		// Ordinary ellipses retain EQU/KEP output. The cometary coordinate
		// policy writes CAR first, then a separate COT OEF at the CAR epoch.
		// COT uses (q,e,i,node,peri,true anomaly); angular values and both
		// covariance axes are exported in degrees. Unavailable covariance
		// is marked explicitly, never replaced by zeros or KEP/EQU fields.
		OrbitElements,
		StateVectors,
		Covariance,
		Residuals,
		Observations,
		All
	};

	struct ExportContext {
		std::string objectName;
		const FitResult* fitResult{nullptr};
		std::span<const ObsNode> observations{};
		AstrometricWeightProfile weightProfile{AstrometricWeightProfile::Vfcc17};
		bool useWeights{true};
		bool cometPhotometry{false};
		PhotometryOptions photometry;
		std::filesystem::path rootDirectory{};
		std::optional<double> gaiaPhotocenterRadiusKm;
	};

	struct ExportReport {
		std::filesystem::path directory;
		// Outputs written by this save call, not an inventory of directory.
		// Unrequested files from earlier exports can remain in the directory.
		std::vector<std::filesystem::path> files;
	};

	class ORBITAL_CORE_API Exporter final {
	   public:
		[[nodiscard]] static ExportReport save(const ExportContext& context,
		                                       std::initializer_list<ExportProduct> products);

		[[nodiscard]] static ExportReport save(const ExportContext& context, std::span<const ExportProduct> products);
	};

}  // namespace orbital_core
