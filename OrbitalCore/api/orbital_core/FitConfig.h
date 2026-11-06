#pragma once
#include "orbital_core/AstrometricWeights.h"
#include "orbital_core/ForceModelPolicy.h"
#include <cstddef>
#include <cstdint>
#include <optional>

namespace orbital_core {
	enum class NonGravModelKind : int {
		None = 0,
		MarsdenH2O = 1,
		MarsdenCO = 2,
		InverseSquare = 3,
		Yarkovsky = 4,
	};

	struct NBodySettings {
		bool enabled = false;
		double stepSize = 10.0;
		uint16_t planetMask = 0x18FF;  // Mercury..Neptune, Moon (bit 11), Pluto (bit 12).
		bool useSb441 = false;
		// NAIF identifier of the body being fitted; zero means unknown.
		int targetNaifId = 0;
		InactiveInteriorMassPolicy inactiveInteriorMassPolicy{InactiveInteriorMassPolicy::ExplicitOnly};
	};

	struct NonGravSettings {
		NonGravModelKind model{NonGravModelKind::None};
		bool enabled{false};
		std::uint8_t solveMask{0};  // 0 => model default.
		bool useInitialValues{false};
		double initialA1{0.0};
		double initialA2{0.0};
		double initialA3{0.0};
		double initialDT{0.0};
	};

	struct SigmaClipSettings {
		bool enabled{false};
		// Reject above this radius (post-fit chi-square). Recovery uses
		// 9.21/10 of the squared threshold, retaining a hysteresis gap.
		double thresholdNSigma{3.0};
		// Maximum selection/refit passes per force-model stage.
		int maxPasses{30};
		// Optional absolute floor; 0 uses count/span and positive-DOF safeguards.
		std::size_t minObsKept{0};
	};

	struct PhotometryOptions {
		double asteroidG{0.15};
		double cometK1{10.0};
		double cometK2{10.0};
		double cometNuclearPhaseCoefficient{0.0};  // mag/degree
	};

	struct FitConfig {
		bool useWeights{false};
		NBodySettings nbody;
		NonGravSettings nonGrav;
		SigmaClipSettings sigmaClip;
		AstrometricWeightProfile weightProfile{AstrometricWeightProfile::Vfcc17};
		ObservationWeightingOptions obsWeighting;
		PhotometryOptions photometry;
		// Physical radius for Gaia FPR photocenter correction; absent means that
		// the size-dependent correction is not available.
		std::optional<double> gaiaPhotocenterRadiusKm;
	};
}  // namespace orbital_core
