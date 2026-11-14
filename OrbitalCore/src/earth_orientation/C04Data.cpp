#include "EopData.h"
#include <array>
#include <filesystem>
#include <string>
#include <vector>
#include "orbital_core/Types.h"
#include "logger/Logger.h"
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace orbital_core::earth_orientation_detail {
	namespace {
		constexpr int degree = 15;
		constexpr int unused = 3;
		constexpr int overlap = 3;
		constexpr int stride = 6;
		constexpr int leading = 18;
		constexpr double besselianYear = 365.242198781;
		constexpr double arcsecToRad = Constants::DEG_TO_RAD / 3600.0;

		double ttAtUtcDay(int mjdUtc) noexcept
		{
			return mjdUtc + Convertors::utcToTdtSeconds(mjdUtc + 2400000.5) / 86400.0;
		}

		std::array<double, 16> chebyshev(double x) noexcept
		{
			std::array<double, 16> p{};
			p[0] = 1.0;
			p[1] = x;
			for (int i = 2; i <= degree; ++i)
				p[i] = (2.0 * x) * p[i - 1] - p[i - 2];

			return p;
		}

		double smooth(double x, bool transition) noexcept
		{
			if (x <= 0.0)
				return 0.0;

			if (x >= 1.0)
				return 1.0;

			const double u = 2.0 * x - 1.0;
			const double u2 = u * u;
			if (transition)
				return 0.5 + u * ((3.0 / 16.0 * u2 - 5.0 / 8.0) * u2 + 15.0 / 16.0);

			constexpr double c[] = {3003.0 / 2048, -3003.0 / 1024, 9009.0 / 2048, -2145.0 / 512,
			                        5005.0 / 2048, -819.0 / 1024,  231.0 / 2048};
			double y = c[6];
			for (int i = 5; i >= 0; --i)
				y = y * u2 + c[i];

			return 0.5 + u * y;
		}

		double field(const std::string& line, std::size_t start, std::size_t length)
		{
			std::istringstream stream(line.substr(start, length));
			stream.imbue(std::locale::classic());

			double value;
			if (!(stream >> value) || !std::isfinite(value))
				throw std::runtime_error("invalid C04 field");

			return value;
		}
	}  // namespace

	C04Data::C04Data(const std::filesystem::path& directory, EarthOrientationProfile profile)
	{
		try {
			load(directory, profile);
			LOG_INFO << "[EarthOrientation] C04/IAU1980, " << m_records.size() << " records incl. extensions; TT MJD "
			         << m_begin << ".." << m_end
			         << "; ccor=eopc04.cnc; extrapolation=linear UT1 + ccor, zero observed pole"
			         << "; record=" << (profile == EarthOrientationProfile::IAU1980 ? "IAU1980" : "full-C04")
			         << "; source=" << directory.string();
		} catch (const std::exception& e) {
			m_records.clear();
			m_segments.clear();
			LOG_WARN << "[EarthOrientation] " << e.what() << "; source=" << directory.string()
			         << "; J2 falls back to NOMINAL pole";
		}
	}

	void C04Data::load(const std::filesystem::path& directory, EarthOrientationProfile profile)
	{
		std::ifstream data(directory / "eopc04.62-now");
		if (!data)
			throw std::runtime_error("C04 file is missing");

		std::string line;
		std::vector<Record> observed;
		int firstDay = 0, lastDay = 0;
		while (std::getline(data, line)) {
			if (line.size() < 87 || line[0] < '0' || line[0] > '9')
				continue;

			const int mjd = static_cast<int>(field(line, 12, 7));
			if (mjd <= 41334)
				continue;

			if (!observed.empty() && mjd != lastDay + 1)
				throw std::runtime_error("non-daily C04 sequence");

			if (observed.empty())
				firstDay = mjd;

			lastDay = mjd;
			if (profile == EarthOrientationProfile::IAU1980) {
				line.resize(80);
				line.resize(87, ' ');
			}

			const double ttMinusUtc = Convertors::utcToTdtSeconds(mjd + 2400000.5);
			observed.push_back({mjd + ttMinusUtc / 86400.0,
			                    {field(line, 19, 11), field(line, 30, 11), ttMinusUtc - field(line, 41, 12),
			                     field(line, 65, 11), field(line, 76, 11)}});
		}

		if (observed.size() < 32)
			throw std::runtime_error("C04 series is too short");

		std::ifstream corrections(directory / "eopc04.cnc");
		if (!corrections)
			throw std::runtime_error("C04 consistency file is missing");

		corrections.imbue(std::locale::classic());
		while (std::getline(corrections, line) && !line.starts_with("-----")) {
		}

		for (auto& value : m_correction)
			corrections >> value;

		for (auto& value : m_trend)
			corrections >> value;

		if (!(corrections >> m_correctionEpoch) || !std::isfinite(m_correctionEpoch))
			throw std::runtime_error("invalid C04 consistency file");

		for (int i = 0; i < 5; ++i) {
			if (!std::isfinite(m_correction[i]) || !std::isfinite(m_trend[i]))
				throw std::runtime_error("nonfinite C04 consistency correction");

			const double scale = i == 2 ? 1e4 : 1e3;
			m_correction[i] /= scale;
			m_trend[i] /= scale;
		}

		const int count = static_cast<int>(observed.size());
		const int intervals = (count - degree) / stride;
		const int endOffset = degree - overlap - unused + (intervals - 1) * stride;
		const int trailing = std::max(2, lastDay - firstDay - endOffset + 2 * stride);
		const int extendedIntervals = (count + trailing - degree) / stride;
		const int beginDay = firstDay + unused + overlap - leading;
		const int endDay = firstDay + degree - overlap - unused + (extendedIntervals - 1) * stride;

		m_begin = ttAtUtcDay(beginDay);
		m_end = ttAtUtcDay(endDay);

		const auto& first = observed.front();
		const auto& last = observed.back();
		m_observedBegin = first.mjdTt;
		m_observedEnd = last.mjdTt;

		m_ut1Slope = (last.values[2] - first.values[2]) / (last.mjdTt - first.mjdTt);
		m_ut1Begin = first.values[2] + m_ut1Slope * (m_begin - first.mjdTt);
		m_ut1End = last.values[2] + m_ut1Slope * (m_end - last.mjdTt);
		m_records.reserve(leading + count + trailing);

		for (int i = 0; i < leading; ++i) {
			const double t = ttAtUtcDay(firstDay - leading + i);
			Values v{};
			if (profile == EarthOrientationProfile::FullC04) {
				const double weight = smooth((t - m_begin) / (first.mjdTt - m_begin), true);
				for (int j = 0; j < 5; ++j)
					v[j] = weight * first.values[j];
			}

			v[2] = m_ut1Begin + m_ut1Slope * (t - m_begin);
			m_records.push_back({t, v});
		}

		m_records.insert(m_records.end(), observed.begin(), observed.end());
		for (int i = 1; i <= trailing; ++i) {
			const double t = ttAtUtcDay(lastDay + i);
			const double weight = smooth((t - m_end) / (last.mjdTt - m_end), true);

			Values v{};
			for (int j = 0; j < 5; ++j)
				v[j] = weight * last.values[j];

			v[2] = m_ut1End + m_ut1Slope * (t - m_end);
			m_records.push_back({t, v});
		}

		const int segmentCount = (static_cast<int>(m_records.size()) - degree) / stride;
		m_segments.resize(segmentCount);

		for (int segment = 0; segment < segmentCount; ++segment) {
			const int start = segment * stride;
			const double t0 = m_records[start].mjdTt;
			const double factor = 2.0 / (m_records[start + degree].mjdTt - t0);

			Eigen::Matrix<double, 16, 16> design;
			Eigen::Matrix<double, 16, 5> values;

			for (int row = 0; row <= degree; ++row) {
				const auto p = chebyshev(factor * (m_records[start + row].mjdTt - t0) - 1.0);
				for (int col = 0; col <= degree; ++col)
					design(row, col) = p[col];
				for (int col = 0; col < 5; ++col)
					values(row, col) = m_records[start + row].values[col];
			}

			const Eigen::Matrix<double, 16, 5> coefficients = design.partialPivLu().solve(values);
			if (!coefficients.allFinite())
				throw std::runtime_error("invalid C04 polynomial coefficients");

			for (int field = 0; field < 5; ++field)
				for (int power = 0; power <= degree; ++power)
					m_segments[segment].coefficients[field][power] = coefficients(power, field);
		}
	}

	C04Data::Values C04Data::interpolate(double t) const noexcept
	{
		const auto upper = std::upper_bound(m_records.begin(), m_records.end(), t,
		                                    [](double epoch, const Record& record) { return epoch < record.mjdTt; });

		const int interval = static_cast<int>(upper - m_records.begin()) - 1;
		const int primary = (interval - unused) / stride;

		if (primary < 0 || primary >= static_cast<int>(m_segments.size())) {
			Values values{};
			values[2] =
			    t < m_begin + 1e-6 ? m_ut1Begin + m_ut1Slope * (t - m_begin) : m_ut1End + m_ut1Slope * (t - m_end);
			return values;
		}

		const auto evaluate = [&](int segment) {
			const int start = segment * stride;
			const double t0 = m_records[start].mjdTt;
			const auto p = chebyshev((2.0 / (m_records[start + degree].mjdTt - t0)) * (t - t0) - 1.0);

			Values values{};
			for (int field = 0; field < 5; ++field)
				for (int power = 0; power <= degree; ++power)
					values[field] += m_segments[segment].coefficients[field][power] * p[power];

			return values;
		};

		Values result = evaluate(primary);
		const int start = primary * stride;
		if (interval - start + 1 <= unused + overlap) {
			if (primary == 0)
				return result;
			const Values previous = evaluate(primary - 1);
			const double t1 = m_records[start + unused].mjdTt;
			const double t2 = m_records[start + unused + overlap].mjdTt;
			const double weight = smooth((t - t1) / (t2 - t1), false);
			for (int i = 0; i < 5; ++i)
				result[i] = weight * result[i] + (1.0 - weight) * previous[i];
		}

		return result;
	}

	bool C04Data::observedTtMinusUt1(double mjdTt, double& seconds) const noexcept
	{
		if (!coversObserved(mjdTt))
			return false;
		const auto upper = std::lower_bound(m_records.begin(), m_records.end(), mjdTt,
		                                    [](const Record& record, double epoch) { return record.mjdTt < epoch; });
		if (upper == m_records.end())
			return false;
		if (upper->mjdTt == mjdTt || upper == m_records.begin()) {
			seconds = upper->values[2];
			return true;
		}
		const auto& before = *(upper - 1);
		const double fraction = (mjdTt - before.mjdTt) / (upper->mjdTt - before.mjdTt);
		seconds = before.values[2] + fraction * (upper->values[2] - before.values[2]);
		return std::isfinite(seconds);
	}

	EarthOrientationSample C04Data::lookupMjdTt(double mjdTt) const noexcept
	{
		if (!isAvailable() || !std::isfinite(mjdTt))
			return {};

		const bool extrapolated = mjdTt < m_begin || mjdTt > m_end;
		Values values{};
		if (extrapolated) {
			values[2] =
			    mjdTt < m_begin ? m_ut1Begin + m_ut1Slope * (mjdTt - m_begin) : m_ut1End + m_ut1Slope * (mjdTt - m_end);
		}

		else {
			values = interpolate(mjdTt);
		}

		const double years = (1900.0 + (mjdTt - 15019.81352) / besselianYear) - m_correctionEpoch;
		for (int i = 0; i < 5; ++i)
			values[i] = (values[i] + m_correction[i]) + m_trend[i] * years;

		return {values[0] * arcsecToRad,
		        values[1] * arcsecToRad,
		        values[2],
		        Iau1980PoleOffsets{values[3] * arcsecToRad, values[4] * arcsecToRad},
		        true,
		        extrapolated,
		        true};
	}

}  // namespace orbital_core::earth_orientation_detail
