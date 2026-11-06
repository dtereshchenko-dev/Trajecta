#include "EopData.h"
#include "logger/Logger.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <string>

namespace orbital_core::earth_orientation_detail {

	DeltaTData::DeltaTData(const std::filesystem::path& file)
	{
		std::ifstream input(file);
		if (!input) {
			LOG_WARN << "[EarthOrientation] historical DeltaT unavailable: " << file.string();
			return;
		}
		bool header = true;
		std::size_t lineNumber = 0;
		std::string line;
		while (std::getline(input, line)) {
			++lineNumber;
			if (header) {
				header = !line.starts_with("----------");
				continue;
			}
			if (line.find_first_not_of(" \t\r") == std::string::npos)
				continue;
			std::istringstream row(line);
			row.imbue(std::locale::classic());
			int d, m, y;
			double seconds;
			using namespace std::chrono;
			if (row >> d >> m >> y >> seconds) {
				const year_month_day date{year{y}, month{static_cast<unsigned>(m)}, day{static_cast<unsigned>(d)}};
				row >> std::ws;
				const double epoch = static_cast<double>(sys_days(date).time_since_epoch().count()) + 40587;
				if (date.ok() && std::isfinite(seconds) && row.eof() &&
				    (m_records.empty() ||
				     (epoch > m_records.back().mjdUt1 &&
				      epoch + seconds / 86400.0 > m_records.back().mjdUt1 + m_records.back().seconds / 86400.0))) {
					m_records.push_back({epoch, seconds});
					continue;
				}
			}
			m_records.clear();
			LOG_WARN << "[EarthOrientation] invalid historical DeltaT table: " << file.string() << ":" << lineNumber;
			return;
		}
		if (header || m_records.size() < 2)
			m_records.clear();
	}

	bool DeltaTData::offsetSeconds(double mjdUt1, double& seconds) const noexcept
	{
		if (!std::isfinite(mjdUt1) || m_records.empty() || mjdUt1 < m_records.front().mjdUt1 ||
		    mjdUt1 > m_records.back().mjdUt1)
			return false;
		const auto next = std::upper_bound(m_records.begin(), m_records.end(), mjdUt1,
		                                   [](double t, const Record& r) { return t < r.mjdUt1; });
		const auto& first = *std::prev(next);
		if (next == m_records.end() || mjdUt1 == first.mjdUt1) {
			seconds = first.seconds;
			return true;
		}
		const double weight = (next->mjdUt1 - mjdUt1) / (next->mjdUt1 - first.mjdUt1);
		seconds = weight * first.seconds + (1.0 - weight) * next->seconds;
		return true;
	}

	bool DeltaTData::ttToUt1Mjd(double mjdTt, double& mjdUt1) const noexcept
	{
		const auto tt = [](const Record& r) { return r.mjdUt1 + r.seconds / 86400.0; };
		if (!std::isfinite(mjdTt) || m_records.empty() || mjdTt < tt(m_records.front()) || mjdTt > tt(m_records.back()))
			return false;
		const auto next = std::upper_bound(m_records.begin(), m_records.end(), mjdTt,
		                                   [&](double t, const Record& r) { return t < tt(r); });
		const auto& first = *std::prev(next);
		double seconds = first.seconds;
		if (next != m_records.end()) {
			const double fraction = (mjdTt - tt(first)) / (tt(*next) - tt(first));
			seconds += fraction * (next->seconds - first.seconds);
		}
		mjdUt1 = mjdTt - seconds / 86400.0;
		return true;
	}
}  // namespace orbital_core::earth_orientation_detail
