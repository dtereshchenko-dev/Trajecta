#pragma once

#include "EarthOrientation.h"
#include <filesystem>

#include <array>
#include <vector>

namespace orbital_core::earth_orientation_detail {

	// ET-UT.dat: day, month, year, TT-UT1 seconds after the dashed header separator.
	// Immutable, strictly increasing samples; no extrapolation or invented fallback.
	class DeltaTData {
	   public:
		explicit DeltaTData(const std::filesystem::path& file);
		[[nodiscard]] bool offsetSeconds(double mjdUt1, double& seconds) const noexcept;
		[[nodiscard]] bool ttToUt1Mjd(double mjdTt, double& mjdUt1) const noexcept;

		[[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }

	   private:
		struct Record {
			double mjdUt1, seconds;
		};

		std::vector<Record> m_records;
	};

	class C04Data {
	   public:
		explicit C04Data(const std::filesystem::path& directory, EarthOrientationProfile profile);
		[[nodiscard]] EarthOrientationSample lookupMjdTt(double mjdTt) const noexcept;

		[[nodiscard]] bool isAvailable() const noexcept { return !m_segments.empty(); }

		[[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }

	   private:
		using Values = std::array<double, 5>;

		struct Record {
			double mjdTt;
			Values values;
		};

		struct Segment {
			std::array<std::array<double, 16>, 5> coefficients;
		};

		std::vector<Record> m_records;
		std::vector<Segment> m_segments;
		Values m_correction{};
		Values m_trend{};

		double m_correctionEpoch{0.0};
		double m_begin{0.0};
		double m_end{0.0};
		double m_ut1Begin{0.0};
		double m_ut1End{0.0};
		double m_ut1Slope{0.0};

		void load(const std::filesystem::path& directory, EarthOrientationProfile profile);
		[[nodiscard]] Values interpolate(double mjdTt) const noexcept;
	};

	class FinalsData {
	   public:
		explicit FinalsData(const std::filesystem::path& file);
		[[nodiscard]] EarthOrientationSample lookupMjdTt(double mjdTt) const noexcept;
		[[nodiscard]] bool ttToUt1Mjd(double mjdTt, double& mjdUt1) const noexcept;

		[[nodiscard]] bool ut1ToTtSeconds(double mjdUt1, double& seconds) const noexcept
		{
			return clockOffsetAtDailyEpoch(mjdUt1, seconds);
		}

		[[nodiscard]] bool isAvailable() const noexcept { return !m_records.empty(); }

		[[nodiscard]] std::size_t size() const noexcept { return m_records.size(); }

		[[nodiscard]] double firstMjdUtc() const noexcept { return m_records.empty() ? 0 : m_records.front().mjdUtc; }

		[[nodiscard]] double lastMjdUtc() const noexcept { return m_records.empty() ? 0 : m_records.back().mjdUtc; }

	   private:
		struct Record {
			bool haveOffsets;
			double mjdUtc;
			double mjdTt;
			double xpRad;
			double ypRad;
			double ttMinusUt1Seconds;
			Iau2000PoleOffsets offsets;
		};

		std::vector<Record> m_records;
		[[nodiscard]] bool clockOffsetAtDailyEpoch(double mjd, double& ttMinusUt1) const noexcept;
	};

}  // namespace orbital_core::earth_orientation_detail
