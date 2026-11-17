#include "EarthPckPoleCache.h"

#include "orbital_core/Constants.h"

#include <SpiceUsr.h>
#include <algorithm>
#include <cmath>
#include <utility>

namespace orbital_core {
	namespace {
		constexpr SpiceInt kEarthFrame = 3000;
		constexpr SpiceInt kEclipticJ2000 = 17;
		constexpr SpiceInt kChebyshevAngles = 2;
		constexpr double kObliquityCos = 0.91748206206918181;
		constexpr double kObliquitySin = 0.39777715593191371;

		double evaluateChebyshev(const double* coefficients, std::size_t count, double t) noexcept
		{
			double b1 = 0.0, b2 = 0.0;
			for (std::size_t i = count - 1; i > 0; --i) {
				const double b = 2.0 * t * b1 - b2 + coefficients[i];
				b2 = b1;
				b1 = b;
			}
			return t * b1 - b2 + coefficients[0];
		}
	}  // namespace

	bool EarthPckPoleCache::addFile(const std::filesystem::path& path)
	{
		if (!m_valid)
			return false;

		SpiceInt handle = 0;
		dafopr_c(path.string().c_str(), &handle);
		if (failed_c()) {
			reset_c();
			m_valid = false;
			return false;
		}

		bool success = true;
		try {
			dafbfs_c(handle);
			SpiceBoolean found = SPICEFALSE;
			while (!failed_c()) {
				daffna_c(&found);
				if (failed_c() || !found)
					break;
				double descriptor[5]{}, bounds[2]{};
				SpiceInt fields[5]{};
				dafgs_c(descriptor);
				dafus_c(descriptor, 2, 5, bounds, fields);
				if (failed_c())
					break;
				if (fields[0] != kEarthFrame)
					continue;

				Segment segment;
				segment.start = bounds[0];
				segment.end = bounds[1];
				if (!std::isfinite(segment.start) || !std::isfinite(segment.end) || segment.start > segment.end) {
					success = false;
					break;
				}

				// Retain unsupported segments in priority order, so they fall back to CSPICE.
				if (fields[1] != kEclipticJ2000 || fields[2] != kChebyshevAngles) {
					m_segments.push_back(std::move(segment));
					continue;
				}

				if (fields[3] < 1 || fields[4] < fields[3] || fields[4] - fields[3] < 8) {
					success = false;
					break;
				}
				const auto length = static_cast<std::size_t>(fields[4] - fields[3]) + 1;
				const auto fileBytes = std::filesystem::file_size(path);
				if (length > fileBytes / sizeof(double)) {
					success = false;
					break;
				}
				segment.data.resize(length);
				dafgda_c(handle, fields[3], fields[4], segment.data.data());
				if (failed_c())
					break;
				segment.initial = segment.data[length - 4];
				segment.interval = segment.data[length - 3];
				const double recordSize = segment.data[length - 2];
				const double recordCount = segment.data[length - 1];
				if (!std::isfinite(segment.initial) || !std::isfinite(segment.interval) || segment.interval <= 0 ||
				    !std::isfinite(recordSize) || !std::isfinite(recordCount) || recordSize < 5 || recordCount < 1 ||
				    recordSize > static_cast<double>(length) || recordCount > static_cast<double>(length) ||
				    std::floor(recordSize) != recordSize || std::floor(recordCount) != recordCount) {
					success = false;
					break;
				}
				segment.recordSize = static_cast<std::size_t>(recordSize);
				segment.recordCount = static_cast<std::size_t>(recordCount);
				if ((segment.recordSize - 2) % 3 != 0 || segment.recordCount != (length - 4) / segment.recordSize ||
				    (length - 4) % segment.recordSize != 0) {
					success = false;
					break;
				}
				m_segments.push_back(std::move(segment));
			}
			if (failed_c()) {
				reset_c();
				success = false;
			}
		} catch (...) {
			success = false;
		}
		if (failed_c()) {
			reset_c();
			success = false;
		}
		dafcls_c(handle);
		if (failed_c()) {
			reset_c();
			success = false;
		}
		if (!success) {
			m_valid = false;
			m_segments.clear();
		}
		return success;
	}

	bool EarthPckPoleCache::pole(double mjdTdb, Vec3D& result) const noexcept
	{
		if (!m_valid || !std::isfinite(mjdTdb))
			return false;
		const double et = (mjdTdb - 51544.5) * Constants::SECONDS_PER_DAY;
		if (!std::isfinite(et))
			return false;
		for (auto it = m_segments.rbegin(); it != m_segments.rend(); ++it) {
			const auto& segment = *it;
			if (et < segment.start || et > segment.end)
				continue;
			if (segment.data.empty())
				return false;
			const double record = (et - segment.initial) / segment.interval;
			if (!std::isfinite(record))
				return false;
			const auto index = record <= 0                                          ? std::size_t{0}
			                   : record >= static_cast<double>(segment.recordCount) ? segment.recordCount - 1
			                                                                        : static_cast<std::size_t>(record);
			const double* data = segment.data.data() + index * segment.recordSize;
			if (!std::isfinite(data[0]) || !std::isfinite(data[1]) || data[1] <= 0)
				return false;
			const double t = (et - data[0]) / data[1];
			const std::size_t count = (segment.recordSize - 2) / 3;
			const double phi = evaluateChebyshev(data + 2, count, t);
			const double delta = evaluateChebyshev(data + 2 + count, count, t);
			if (!std::isfinite(phi) || !std::isfinite(delta))
				return false;
			const double x = std::sin(phi) * std::sin(delta);
			const double y = -std::cos(phi) * std::sin(delta);
			const double z = std::cos(delta);
			result = {x, kObliquityCos * y - kObliquitySin * z, kObliquitySin * y + kObliquityCos * z};
			return true;
		}
		return false;
	}
}  // namespace orbital_core
