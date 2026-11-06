#pragma once

#include "FitConfig.h"
#include "FitResult.h"
#include "FitRunControl.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <cstdint>
#include <optional>
#include <span>

namespace orbital_core {
	enum class ObservationDisposition { Auto, Include, Exclude };

	struct ObservationOverride {
		ObservationId id;
		ObservationDisposition disposition{ObservationDisposition::Auto};
		bool operator==(const ObservationOverride&) const = default;
	};

	enum class FitSelectionMode { Configured, Fixed };

	enum class FitUncertaintyMode { Compute, Skip };

	struct FitRunOptions {
		FitSelectionMode selection{FitSelectionMode::Configured};
		FitUncertaintyMode uncertainty{FitUncertaintyMode::Compute};
		bool evaluateOnly{false};
	};

	struct FitRequest {
		Orbit_p initialOrbit;
		std::span<const ObsNode> observations;
		std::span<const std::uint8_t> initialIncluded{};
		const NonGravModel* warmStartNonGrav{nullptr};
		FitRunOptions options{};
		std::optional<double> outputEpochJdTdb;
		std::span<const ObservationOverride> observationOverrides{};
	};

	class ORBITAL_CORE_API FitSession {
	   public:
		explicit FitSession(FitConfig cfg);
		[[nodiscard]] FitResult optimize(const FitRequest& request, const FitRunControl& control = {}) const;

	   private:
		FitConfig config_;
	};
}  // namespace orbital_core
