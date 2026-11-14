#pragma once

#include "orbital_core/OrbitalCoreExport.h"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <variant>

namespace orbital_core {

	struct Iau1980PoleOffsets {
		double dpsiRad = 0.0;
		double depsRad = 0.0;
	};

	struct Iau2000PoleOffsets {
		double dxRad = 0.0;
		double dyRad = 0.0;
	};

	struct EarthOrientationSample {
		double xpRad = 0.0;
		double ypRad = 0.0;
		double ttMinusUt1Seconds = 0.0;
		std::variant<Iau1980PoleOffsets, Iau2000PoleOffsets> poleOffsets;
		bool valid = false;
		bool extrapolated = false;
		bool poleOffsetsAvailable = false;
	};

	enum class EarthOrientationProfile {
		IAU1980,  // C04/IAU1980, native 80-column records and consistency corrections.
		FullC04,  // Same interpolation/corrections, without native record truncation.
	};

	class ORBITAL_CORE_API EarthOrientation {
	   public:
		explicit EarthOrientation(const std::filesystem::path& dataDirectory,
		                          EarthOrientationProfile profile = EarthOrientationProfile::IAU1980,
		                          const std::filesystem::path& orientationDirectory = {});
		EarthOrientation(const EarthOrientation&) = delete;
		EarthOrientation& operator=(const EarthOrientation&) = delete;
		~EarthOrientation();

		[[nodiscard]] EarthOrientationSample lookupMjdTt(double mjdTt) const noexcept;
		[[nodiscard]] bool coversObservedMjdTt(double mjdTt) const noexcept;
		[[nodiscard]] bool ttToUt1Mjd(double mjdTt, double& mjdUt1) const noexcept;
		[[nodiscard]] bool ut1ToTtSeconds(double mjdUt1, double& seconds) const noexcept;
		[[nodiscard]] EarthOrientationProfile profile() const noexcept;
		[[nodiscard]] bool isAvailable() const noexcept;
		[[nodiscard]] std::size_t size() const noexcept;

	   private:
		struct Impl;
		std::unique_ptr<Impl> m_impl;
	};

	[[nodiscard]] ORBITAL_CORE_API const EarthOrientation& earthOrientation();

}  // namespace orbital_core
