#include "EopData.h"
#include <filesystem>
#include <string>
#include "logger/Logger.h"
#include "orbital_core/Types.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <string_view>

namespace orbital_core::earth_orientation_detail {
	namespace {

		bool number(std::string_view line, std::size_t start, std::size_t length, double& value)
		{
			if (start + length > line.size())
				return false;
			auto field = line.substr(start, length);
			const auto first = field.find_first_not_of(' ');
			if (first == std::string_view::npos)
				return false;
			field.remove_prefix(first);
			field = field.substr(0, field.find_last_not_of(' ') + 1);
			const auto result = std::from_chars(field.data(), field.data() + field.size(), value);
			return result.ec == std::errc{} && result.ptr == field.data() + field.size() && std::isfinite(value);
		}

		bool available(char flag) { return flag == 'I' || flag == 'P'; }

		EarthOrientationSample fallback(double mjdTt, bool extrapolated)
		{
			EarthOrientationSample result;
			result.poleOffsets = Iau2000PoleOffsets{};
			result.extrapolated = extrapolated;
			if (std::isfinite(mjdTt)) {
				double mjdUtc = mjdTt;
				for (int i = 0; i < 3; ++i) {
					result.ttMinusUt1Seconds = Convertors::utcToTdtSeconds(mjdUtc + 2400000.5);
					mjdUtc = mjdTt - result.ttMinusUt1Seconds / 86400.0;
				}
			}
			return result;
		}

	}  // namespace

	FinalsData::FinalsData(const std::filesystem::path& file)
	{
		std::ifstream input(file);
		if (!input) {
			LOG_WARN << "[EarthOrientation] finals2000A missing: " << file.string()
			         << "; zero observed pole corrections, UT1=UTC";
			return;
		}

		std::size_t skipped = 0, predicted = 0, offsets = 0;
		std::string line;
		while (std::getline(input, line)) {
			double mjdUtc, xp, yp, ut1MinusUtc, dx = 0.0, dy = 0.0;
			if (line.size() < 68 || !available(line[16]) || !available(line[57]) || !number(line, 7, 8, mjdUtc) ||
			    !number(line, 18, 9, xp) || !number(line, 37, 9, yp) || !number(line, 58, 10, ut1MinusUtc)) {
				++skipped;
				continue;
			}
			const bool haveOffsets =
			    line.size() >= 125 && available(line[95]) && number(line, 97, 9, dx) && number(line, 116, 9, dy);
			const double ttMinusUtc = Convertors::utcToTdtSeconds(mjdUtc + 2400000.5);
			const double mjdTt = mjdUtc + ttMinusUtc / 86400.0;
			if (!m_records.empty() && mjdTt <= m_records.back().mjdTt) {
				m_records.clear();
				LOG_WARN << "[EarthOrientation] non-increasing finals2000A epochs: " << file.string();
				return;
			}
			constexpr double masToRad = Constants::ARCSEC_RAD / 1000.0;
			m_records.push_back(
			    Record{.haveOffsets = haveOffsets,
			           .mjdUtc = mjdUtc,
			           .mjdTt = mjdTt,
			           .xpRad = xp * Constants::ARCSEC_RAD,
			           .ypRad = yp * Constants::ARCSEC_RAD,
			           .ttMinusUt1Seconds = ttMinusUtc - ut1MinusUtc,
			           .offsets = {haveOffsets ? dx * masToRad : 0.0, haveOffsets ? dy * masToRad : 0.0}});
			predicted += line[16] == 'P' || line[57] == 'P';
			offsets += haveOffsets;
		}
		LOG_INFO << "[EarthOrientation] finals2000A records=" << m_records.size() << " predicted=" << predicted
		         << " celestialOffsets=" << offsets << " skipped=" << skipped << "; source=" << file.string();
	}

	EarthOrientationSample FinalsData::lookupMjdTt(double mjdTt) const noexcept
	{
		if (!std::isfinite(mjdTt) || m_records.empty())
			return fallback(mjdTt, false);
		if (mjdTt < m_records.front().mjdTt || mjdTt > m_records.back().mjdTt)
			return fallback(mjdTt, true);

		const auto upper = std::upper_bound(m_records.begin(), m_records.end(), mjdTt,
		                                    [](double epoch, const Record& record) { return epoch < record.mjdTt; });
		const auto& first = *std::prev(upper);
		if (first.mjdTt == mjdTt)
			return {first.xpRad, first.ypRad, first.ttMinusUt1Seconds, first.offsets, true, false, first.haveOffsets};

		const auto& second = *upper;
		const double span = second.mjdTt - first.mjdTt;
		if (span > 1.5)
			return fallback(mjdTt, false);
		const double fraction = (mjdTt - first.mjdTt) / span;
		const auto interpolate = [fraction](double a, double b) { return a + fraction * (b - a); };
		const bool haveOffsets = first.haveOffsets && second.haveOffsets;
		Iau2000PoleOffsets offsets;
		if (haveOffsets) {
			offsets.dxRad = interpolate(first.offsets.dxRad, second.offsets.dxRad);
			offsets.dyRad = interpolate(first.offsets.dyRad, second.offsets.dyRad);
		}
		return {interpolate(first.xpRad, second.xpRad),
		        interpolate(first.ypRad, second.ypRad),
		        interpolate(first.ttMinusUt1Seconds, second.ttMinusUt1Seconds),
		        offsets,
		        true,
		        false,
		        haveOffsets};
	}

	bool FinalsData::clockOffsetAtDailyEpoch(double mjd, double& ttMinusUt1) const noexcept
	{
		if (!std::isfinite(mjd) || m_records.empty() || mjd < m_records.front().mjdUtc || mjd > m_records.back().mjdUtc)
			return false;

		const auto upper = std::upper_bound(m_records.begin(), m_records.end(), mjd,
		                                    [](double epoch, const Record& record) { return epoch < record.mjdUtc; });
		const auto& first = *std::prev(upper);
		const auto& second = upper == m_records.end() ? first : *upper;
		const double span = second.mjdUtc - first.mjdUtc;
		if (span > 1.5 && mjd != first.mjdUtc)
			return false;
		const double fraction = span > 0.0 ? (mjd - first.mjdUtc) / span : 0.0;
		ttMinusUt1 = first.ttMinusUt1Seconds + fraction * (second.ttMinusUt1Seconds - first.ttMinusUt1Seconds);
		return true;
	}

	bool FinalsData::ttToUt1Mjd(double mjdTt, double& mjdUt1) const noexcept
	{
		if (!std::isfinite(mjdTt))
			return false;

		constexpr int maxIterations = 5;
		constexpr double toleranceSeconds = 1.0e-10;
		double estimate = mjdTt;
		for (int iteration = 0; iteration < maxIterations; ++iteration) {
			double offset;
			if (!clockOffsetAtDailyEpoch(estimate, offset)) {
				if (iteration != 0)
					return false;
				const auto sample = lookupMjdTt(mjdTt);
				if (!sample.valid)
					return false;
				offset = sample.ttMinusUt1Seconds;
			}
			const double next = mjdTt - offset / 86400.0;
			const double correctionSeconds = (next - estimate) * 86400.0;
			estimate = next;
			if (std::abs(correctionSeconds) <= toleranceSeconds)
				break;
		}
		mjdUt1 = estimate;
		return true;
	}

}  // namespace orbital_core::earth_orientation_detail
