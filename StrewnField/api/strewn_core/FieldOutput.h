#pragma once

#include "strewn_core/Export.h"
#include <array>
#include <filesystem>
#include <optional>
#include <span>
#include <stop_token>
#include <string>
#include <vector>

namespace strewn {

	enum class KdeMethod { Fixed, Adaptive };
	enum class Basemap { None, EsriImagery };
	enum class ArtifactStatus { Skipped, Complete, Empty, Failed, Cancelled };

	struct MapOptions {
		std::string title{"Strewn field"};
		KdeMethod method{KdeMethod::Adaptive};
		double cellM{25};
		// Increase the requested cell size only when grid/FFT memory limits require it.
		bool autoCellSize{true};
		double fixedSigmaM{60};
		double adaptiveAlpha{0.5};
		std::array<double, 2> adaptiveLimits{0.2, 5.0};
		unsigned adaptiveClasses{16};
		// Latitude/longitude center, north/east half-size. Absent center derives coverage.
		std::optional<std::array<double, 2>> centerDeg;
		std::array<double, 2> halfSizeKm{14, 20};
		// Legacy PNG settings are accepted for archive compatibility, but unused.
		unsigned imageWidth{1600};
		Basemap basemap{Basemap::None};
		bool offline{false};
		std::filesystem::path fontFile;
	};

	struct FieldPoint {
		double latitudeDeg{}, longitudeDeg{}, massKg{}, multiplicity{1}, runWeight{1};
	};

	struct FieldDensity {
		std::vector<double> probability;
		std::array<double, 3> covarianceM2{};  // east/east, east/north, north/north.
		double effectiveSamples{};
		double totalWeight{};
		bool fixedKernelFallback{};
		std::array<double, 3> hpdThresholds{};  // 50%, 80%, 99%.
	};

	struct DensityFields {
		std::size_t rows{}, columns{};
		double cellM{};
		double firstLatitudeDeg{}, firstLongitudeDeg{}, latitudeStepDeg{}, longitudeStepDeg{};
		// Cell centers, south-to-north rows; longitude is the fastest index.
		FieldDensity count, mass;
	};

	struct ArtifactReport {
		ArtifactStatus status{ArtifactStatus::Skipped};
		std::vector<std::filesystem::path> files;
		std::vector<std::string> warnings;
		std::string error;
		double elapsedS{};
	};

	STREWN_CORE_EXPORT void validate(const MapOptions& options);
	STREWN_CORE_EXPORT MapOptions loadMapOptions(const std::filesystem::path& path);
	STREWN_CORE_EXPORT DensityFields estimateField(std::span<const FieldPoint> points, const MapOptions& options,
	                                               std::stop_token cancel = {});
	// Output must be a new directory distinct from the numerical archive.
	// No orbital propagation. Display traces replay and verify every saved realization.
	// cacheDirectory stores only imagery; numerical archives are never modified.
	STREWN_CORE_EXPORT ArtifactReport renderSavedRun(const std::filesystem::path& runDirectory,
	                                                 const std::filesystem::path& outputDirectory,
	                                                 const MapOptions& options, std::stop_token cancel = {},
	                                                 const std::filesystem::path& cacheDirectory = {});

}  // namespace strewn
