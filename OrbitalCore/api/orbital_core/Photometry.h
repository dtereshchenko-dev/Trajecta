#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Residuals.h"
#include "orbital_core/FitConfig.h"

#include <optional>
#include <span>
#include <string>
#include <string_view>

namespace orbital_core {
	struct FitResult;
}

namespace orbital_core::photometry {
	inline constexpr double kDefaultAsteroidSlope = PhotometryOptions{}.asteroidG;

	struct Solution {
		std::optional<double> absoluteMagnitude;
		double rms{-1.0};
		double weightSum{0.0};
		std::size_t count{0};
		double slopeG{kDefaultAsteroidSlope};
	};

	struct Evaluation {
		double sigma{0.0};
		std::optional<double> residual;
		bool selected{false};
	};

	struct CometComponent {
		std::optional<double> absoluteMagnitude;
		double slope{10.0};
		double phaseCoefficient{0.0};
		std::size_t count{0};
	};

	struct CometSolution {
		CometComponent total;
		CometComponent nuclear;
		std::size_t unknownTypeCount{0};
	};

	ORBITAL_CORE_API double resolveSigma(std::string_view magnitudeText, double reportedSigma = 0.0) noexcept;

	// Borrows the full observation geometry when available. Astrometric rejection
	// does not reject photometry; older/partial results fall back to fitted rows.
	ORBITAL_CORE_API const Residuals& observations(const FitResult& result) noexcept;

	ORBITAL_CORE_API Solution fit(std::span<const ::Residual> observations, double slopeG = kDefaultAsteroidSlope);
	// K is an assumed, fixed coefficient, not an estimated uncertainty.
	ORBITAL_CORE_API CometSolution fitComet(std::span<const ::Residual> observations, double slope = 10.0);
	ORBITAL_CORE_API CometSolution fitComet(std::span<const ::Residual> observations, const PhotometryOptions& options);
	ORBITAL_CORE_API Evaluation evaluateComet(const ::Residual& observation, const CometSolution& solution,
	                                          bool includedInFit = true);
	ORBITAL_CORE_API std::optional<double> cometMagnitude(double absoluteMagnitude, double slope, double heliocentricAu,
	                                                      double observerAu, double phaseRad = 0.0,
	                                                      double phaseCoefficient = 0.0) noexcept;
	ORBITAL_CORE_API Evaluation evaluate(const ::Residual& observation, const Solution& solution,
	                                     bool includedInFit = true);

	ORBITAL_CORE_API std::string formatMagnitude(double magnitude, char band, std::string_view magnitudeText);

	// Predicted V magnitude from the H-G law used to reduce observed magnitudes to H:
	// V = H + 5 log10(r Delta) - 2.5 log10((1 - G) Phi1 + G Phi2). Distances in AU,
	// phase angle in radians. Empty for non-finite or non-positive geometry.
	ORBITAL_CORE_API std::optional<double> visualMagnitude(double absoluteMagnitude, double slopeG,
	                                                       double heliocentricAu, double observerAu,
	                                                       double phaseRad) noexcept;

}  // namespace orbital_core::photometry
