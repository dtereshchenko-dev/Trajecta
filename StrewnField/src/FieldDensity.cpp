#include "strewn_core/FieldOutput.h"
#include "VectorMath.h"
#include <unsupported/Eigen/FFT>
#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <limits>
#include <stdexcept>

namespace strewn {
	namespace {

		using Covariance = std::array<double, 3>;
		constexpr double radiusM = 6367e3;

		void cancelled(std::stop_token token)
		{
			if (token.stop_requested())
				throw std::runtime_error("Field estimation cancelled");
		}

		struct Samples {
			std::vector<double> x, y, weight;
			double total{}, squared{}, meanX{}, meanY{};
			Covariance covariance{};
		};

		Samples samples(std::span<const FieldPoint> points, bool mass, double lat0, double lon0)
		{
			Samples s;
			for (const auto& p : points) {
				const double w = p.runWeight * p.multiplicity * (mass ? p.massKg : 1);
				const double x = (p.longitudeDeg - lon0) * detail::deg * radiusM * std::cos(lat0 * detail::deg);
				const double y = (p.latitudeDeg - lat0) * detail::deg * radiusM;
				s.x.push_back(x);
				s.y.push_back(y);
				s.weight.push_back(w);
				s.total += w;
				s.squared += w * w;
				s.meanX += w * x;
				s.meanY += w * y;
			}
			if (!std::isfinite(s.total) || !std::isfinite(s.squared) || s.total <= 0 || s.squared <= 0)
				throw std::invalid_argument("Invalid total field weight");
			s.meanX /= s.total;
			s.meanY /= s.total;
			const double denominator = s.total - s.squared / s.total;
			if (denominator > 0) {
				for (std::size_t i = 0; i < s.x.size(); ++i) {
					const double dx = s.x[i] - s.meanX, dy = s.y[i] - s.meanY;
					s.covariance[0] += s.weight[i] * dx * dx;
					s.covariance[1] += s.weight[i] * dx * dy;
					s.covariance[2] += s.weight[i] * dy * dy;
				}
				const double scott = std::pow(s.total * s.total / s.squared, -1.0 / 3.0) / denominator;
				for (auto& v : s.covariance)
					v *= scott;
			}
			return s;
		}

		std::size_t fftSize(std::size_t n)
		{
			std::size_t size = 1;
			while (size < n && size < 16384)
				size *= 2;
			if (size < n)
				throw std::invalid_argument("KDE convolution dimension exceeds limit");
			return size;
		}

		void fft2(std::vector<std::complex<double>>& data, std::size_t nx, std::size_t ny, bool inverse,
		          std::stop_token cancel)
		{
			Eigen::FFT<double> fft;
			std::vector<std::complex<double>> in, out;
			in.resize(nx);
			out.resize(nx);
			for (std::size_t y = 0; y < ny; ++y) {
				if (y % 16 == 0)
					cancelled(cancel);
				std::copy_n(data.data() + y * nx, nx, in.data());
				if (inverse)
					fft.inv(out, in);
				else
					fft.fwd(out, in);
				std::copy(out.begin(), out.end(), data.begin() + y * nx);
			}
			in.resize(ny);
			out.resize(ny);
			for (std::size_t x = 0; x < nx; ++x) {
				if (x % 16 == 0)
					cancelled(cancel);
				for (std::size_t y = 0; y < ny; ++y)
					in[y] = data[y * nx + x];
				if (inverse)
					fft.inv(out, in);
				else
					fft.fwd(out, in);
				for (std::size_t y = 0; y < ny; ++y)
					data[y * nx + x] = out[y];
			}
		}

		std::vector<double> smooth(const std::vector<double>& hist, const DensityFields& grid, Covariance cov,
		                           double cellM, std::stop_token cancel)
		{
			for (auto& v : cov)
				v /= cellM * cellM;
			if (!std::isfinite(cov[0]) || !std::isfinite(cov[1]) || !std::isfinite(cov[2]) || cov[0] <= 0 ||
			    cov[2] <= 0 || cov[0] > 4e6 || cov[2] > 4e6)
				throw std::invalid_argument("KDE kernel exceeds convolution limits");
			const auto rx = static_cast<std::size_t>(std::ceil(4 * std::sqrt(cov[0]))) + 1;
			const auto ry = static_cast<std::size_t>(std::ceil(4 * std::sqrt(cov[2]))) + 1;
			std::size_t x0 = grid.columns, y0 = grid.rows, x1 = 0, y1 = 0;
			for (std::size_t y = 0; y < grid.rows; ++y)
				for (std::size_t x = 0; x < grid.columns; ++x)
					if (hist[y * grid.columns + x] > 0) {
						x0 = std::min(x0, x);
						x1 = std::max(x1, x);
						y0 = std::min(y0, y);
						y1 = std::max(y1, y);
					}
			std::vector<double> result(hist.size());
			if (x0 == grid.columns)
				return result;
			const auto width = x1 - x0 + 1, height = y1 - y0 + 1;
			const auto nx = fftSize(width + 2 * rx), ny = fftSize(height + 2 * ry);
			if (nx * ny > 8388608)
				throw std::invalid_argument("KDE convolution exceeds memory limit");
			std::vector<std::complex<double>> h(nx * ny), kernel(nx * ny);
			for (std::size_t y = 0; y < height; ++y)
				for (std::size_t x = 0; x < width; ++x)
					h[y * nx + x] = hist[(y + y0) * grid.columns + x + x0];
			const double determinant = cov[0] * cov[2] - cov[1] * cov[1];
			double sum = 0;
			for (std::size_t y = 0; y <= 2 * ry; ++y)
				for (std::size_t x = 0; x <= 2 * rx; ++x) {
					const double dx = static_cast<double>(x) - rx, dy = static_cast<double>(y) - ry;
					const double value =
					    std::exp(-.5 * (cov[2] * dx * dx - 2 * cov[1] * dx * dy + cov[0] * dy * dy) / determinant);
					kernel[y * nx + x] = value;
					sum += value;
				}
			for (auto& v : kernel)
				v /= sum;
			fft2(h, nx, ny, false, cancel);
			fft2(kernel, nx, ny, false, cancel);
			for (std::size_t i = 0; i < h.size(); ++i)
				h[i] *= kernel[i];
			kernel.clear();
			kernel.shrink_to_fit();
			fft2(h, nx, ny, true, cancel);
			const auto west = x0 > rx ? x0 - rx : 0, south = y0 > ry ? y0 - ry : 0;
			const auto east = std::min(grid.columns, x1 + rx + 1), north = std::min(grid.rows, y1 + ry + 1);
			for (std::size_t y = south; y < north; ++y)
				for (std::size_t x = west; x < east; ++x)
					result[y * grid.columns + x] = std::max(0.0, h[(y + ry - y0) * nx + x + rx - x0].real());
			return result;
		}

		FieldDensity density(const Samples& s, const DensityFields& grid, const MapOptions& options,
		                     std::stop_token cancel)
		{
			FieldDensity result;
			result.totalWeight = s.total;
			result.effectiveSamples = s.total * s.total / s.squared;
			auto cov = s.covariance;
			const double det = cov[0] * cov[2] - cov[1] * cov[1];
			result.fixedKernelFallback = options.method != KdeMethod::Fixed &&
			                             (!std::isfinite(det) || det <= 1e-12 * std::max(1.0, cov[0] * cov[2]));
			if (options.method == KdeMethod::Fixed || result.fixedKernelFallback)
				cov = {options.fixedSigmaM * options.fixedSigmaM, 0, options.fixedSigmaM * options.fixedSigmaM};
			result.covarianceM2 = cov;
			std::vector<std::size_t> cells(s.x.size());
			std::vector<double> hist(grid.rows * grid.columns);
			for (std::size_t i = 0; i < s.x.size(); ++i) {
				const double x = s.x[i] / options.cellM + grid.columns / 2.0 + .5;
				const double y = s.y[i] / options.cellM + grid.rows / 2.0 + .5;
				if (x < 0 || y < 0 || x > grid.columns || y > grid.rows)
					throw std::invalid_argument("Impact outside output grid; enlarge map extent");
				cells[i] = std::min(static_cast<std::size_t>(y), grid.rows - 1) * grid.columns +
				           std::min(static_cast<std::size_t>(x), grid.columns - 1);
				hist[cells[i]] += s.weight[i];
			}
			auto values = smooth(hist, grid, cov, options.cellM, cancel);
			if (options.method == KdeMethod::Adaptive && !result.fixedKernelFallback) {
				std::vector<double> pilot(s.x.size());
				double minPositive = std::numeric_limits<double>::max();
				for (std::size_t i = 0; i < pilot.size(); ++i) {
					pilot[i] = values[cells[i]];
					if (pilot[i] > 0)
						minPositive = std::min(minPositive, pilot[i]);
				}
				double logMean = 0;
				for (std::size_t i = 0; i < pilot.size(); ++i) {
					pilot[i] = std::max(minPositive, pilot[i]);
					logMean += s.weight[i] * std::log(pilot[i]);
				}
				logMean /= s.total;
				const double lo = options.adaptiveLimits[0], hi = options.adaptiveLimits[1];
				const double logStep = std::log(hi / lo) / options.adaptiveClasses;
				std::vector<unsigned> classes(pilot.size());
				for (std::size_t i = 0; i < pilot.size(); ++i) {
					const double lambda =
					    std::clamp(std::exp(-options.adaptiveAlpha * (std::log(pilot[i]) - logMean)), lo, hi);
					classes[i] =
					    std::min(static_cast<unsigned>(std::max(0.0, std::floor(std::log(lambda / lo) / logStep))),
					             options.adaptiveClasses - 1);
				}
				std::fill(values.begin(), values.end(), 0);
				for (unsigned c = 0; c < options.adaptiveClasses; ++c) {
					cancelled(cancel);
					std::fill(hist.begin(), hist.end(), 0);
					bool used = false;
					for (std::size_t i = 0; i < cells.size(); ++i)
						if (classes[i] == c) {
							hist[cells[i]] += s.weight[i];
							used = true;
						}
					if (!used)
						continue;
					const double factor = lo * std::exp((c + .5) * logStep);
					auto scaled = cov;
					for (auto& v : scaled)
						v *= factor * factor;
					auto layer = smooth(hist, grid, scaled, options.cellM, cancel);
					for (std::size_t i = 0; i < values.size(); ++i)
						values[i] += layer[i];
				}
			}
			const double total = std::accumulate(values.begin(), values.end(), 0.0);
			if (!std::isfinite(total) || total <= 0)
				throw std::runtime_error("Empty output density");
			for (auto& v : values)
				v /= total;
			auto sorted = values;
			std::sort(sorted.begin(), sorted.end(), std::greater<double>());
			const std::array<double, 3> levels{.5, .8, .99};
			double cumulative = 0;
			std::size_t next = 0;
			for (double v : sorted) {
				cumulative += v;
				while (next < levels.size() && cumulative >= levels[next])
					result.hpdThresholds[next++] = v;
				if (next == levels.size())
					break;
			}
			double edge = 0;
			for (std::size_t x = 0; x < grid.columns; ++x)
				edge = std::max({edge, values[x], values[(grid.rows - 1) * grid.columns + x]});
			for (std::size_t y = 0; y < grid.rows; ++y)
				edge = std::max({edge, values[y * grid.columns], values[y * grid.columns + grid.columns - 1]});
			if (edge >= result.hpdThresholds[2])
				throw std::invalid_argument("99% HPD touches output grid edge; enlarge map extent");
			result.probability = std::move(values);
			return result;
		}

	}  // namespace

	void validate(const MapOptions& o)
	{
		if (o.method != KdeMethod::Fixed && o.method != KdeMethod::Adaptive)
			throw std::invalid_argument("Unknown KDE method");
		if (o.basemap != Basemap::None && o.basemap != Basemap::EsriImagery)
			throw std::invalid_argument("Unknown basemap");
		if (!std::isfinite(o.cellM) || o.cellM <= 0 || !std::isfinite(o.fixedSigmaM) || o.fixedSigmaM <= 0 ||
		    !std::isfinite(o.adaptiveAlpha) || o.adaptiveAlpha < 0 || o.adaptiveAlpha > 1 ||
		    !std::isfinite(o.adaptiveLimits[0]) || !std::isfinite(o.adaptiveLimits[1]) || o.adaptiveLimits[0] <= 0 ||
		    o.adaptiveLimits[1] <= o.adaptiveLimits[0] || o.adaptiveLimits[1] > 10 || o.adaptiveClasses < 1 ||
		    o.adaptiveClasses > 32 || o.imageWidth < 800 || o.imageWidth > 4096)
			throw std::invalid_argument("Invalid map settings");
		if (o.centerDeg && (!std::isfinite((*o.centerDeg)[0]) || !std::isfinite((*o.centerDeg)[1]) ||
		                    std::abs((*o.centerDeg)[0]) >= 85 || std::abs((*o.centerDeg)[1]) > 180))
			throw std::invalid_argument("Invalid map center");
		for (double v : o.halfSizeKm)
			if (!std::isfinite(v) || v <= 0 || v > 1000)
				throw std::invalid_argument("Invalid map extent");
	}

	DensityFields estimateField(std::span<const FieldPoint> points, const MapOptions& options, std::stop_token cancel)
	{
		validate(options);
		cancelled(cancel);
		if (points.empty())
			throw std::invalid_argument("Cannot estimate an empty field");
		double south = 90, north = -90, west = 180, east = -180;
		for (const auto& p : points) {
			if (!std::isfinite(p.latitudeDeg) || std::abs(p.latitudeDeg) >= 85 || !std::isfinite(p.longitudeDeg) ||
			    std::abs(p.longitudeDeg) > 180 || !std::isfinite(p.massKg) || p.massKg <= 0 ||
			    !std::isfinite(p.multiplicity) || p.multiplicity <= 0 || !std::isfinite(p.runWeight) ||
			    p.runWeight <= 0)
				throw std::invalid_argument("Invalid weighted impact");
			south = std::min(south, p.latitudeDeg);
			north = std::max(north, p.latitudeDeg);
			west = std::min(west, p.longitudeDeg);
			east = std::max(east, p.longitudeDeg);
		}
		if (east - west >= 180)
			throw std::invalid_argument("Dateline fields are not supported");
		const auto center = options.centerDeg.value_or(std::array<double, 2>{(south + north) / 2, (west + east) / 2});
		const auto count = samples(points, false, center[0], center[1]);
		const auto mass = samples(points, true, center[0], center[1]);
		auto half = options.halfSizeKm;
		if (!options.centerDeg) {
			const double factor = options.method == KdeMethod::Adaptive ? options.adaptiveLimits[1] : 1;
			const double sx =
			    std::max({std::sqrt(count.covariance[0]), std::sqrt(mass.covariance[0]), options.fixedSigmaM});
			const double sy =
			    std::max({std::sqrt(count.covariance[2]), std::sqrt(mass.covariance[2]), options.fixedSigmaM});
			half = {((north - south) * detail::deg * radiusM / 2 + 6 * sy * factor) / 1000,
			        ((east - west) * detail::deg * radiusM * std::cos(center[0] * detail::deg) / 2 + 6 * sx * factor) /
			            1000};
		}
		auto resolved = options;
		const auto fits = [&](double cell) {
			const double nx = std::floor(2000 * half[1] / cell), ny = std::floor(2000 * half[0] / cell);
			if (!std::isfinite(nx) || !std::isfinite(ny) || nx < 4 || ny < 4)
				return false;
			if (nx * ny > 4000000)
				return false;
			const double factor = options.method == KdeMethod::Adaptive ? std::max(1.0, options.adaptiveLimits[1]) : 1;
			const double sx =
			    options.method == KdeMethod::Fixed
			        ? options.fixedSigmaM
			        : std::max({std::sqrt(count.covariance[0]), std::sqrt(mass.covariance[0]), options.fixedSigmaM});
			const double sy =
			    options.method == KdeMethod::Fixed
			        ? options.fixedSigmaM
			        : std::max({std::sqrt(count.covariance[2]), std::sqrt(mass.covariance[2]), options.fixedSigmaM});
			const double width =
			    std::ceil((east - west) * detail::deg * radiusM * std::cos(center[0] * detail::deg) / cell) + 2;
			const double height = std::ceil((north - south) * detail::deg * radiusM / cell) + 2;
			const double fx = width + 2 * (std::ceil(4 * sx * factor / cell) + 1);
			const double fy = height + 2 * (std::ceil(4 * sy * factor / cell) + 1);
			if (!std::isfinite(fx) || !std::isfinite(fy) || fx > 16384 || fy > 16384)
				return false;
			return fftSize(static_cast<std::size_t>(fx)) * fftSize(static_cast<std::size_t>(fy)) <= 8388608;
		};
		if (resolved.autoCellSize) {
			unsigned attempts = 0;
			while (!fits(resolved.cellM)) {
				if (++attempts > 64 || !std::isfinite(half[0]) || !std::isfinite(half[1]) ||
				    2000 * std::min(half[0], half[1]) / resolved.cellM < 8)
					throw std::invalid_argument("Cannot fit map grid and KDE into memory limits");
				resolved.cellM *= 2;
			}
		}
		const double ny = std::floor(2000 * half[0] / resolved.cellM), nx = std::floor(2000 * half[1] / resolved.cellM);
		if (!std::isfinite(nx) || !std::isfinite(ny) || nx < 4 || ny < 4 || nx * ny > 4000000)
			throw std::invalid_argument("Output grid must contain 4..4000000 cells with each dimension >=4");
		DensityFields grid;
		grid.cellM = resolved.cellM;
		grid.rows = static_cast<std::size_t>(ny);
		grid.columns = static_cast<std::size_t>(nx);
		grid.latitudeStepDeg = resolved.cellM / (radiusM * detail::deg);
		grid.longitudeStepDeg = grid.latitudeStepDeg / std::cos(center[0] * detail::deg);
		grid.firstLatitudeDeg = center[0] - ny / 2 * grid.latitudeStepDeg;
		grid.firstLongitudeDeg = center[1] - nx / 2 * grid.longitudeStepDeg;
		if (grid.firstLatitudeDeg <= -85 || grid.firstLatitudeDeg + ny * grid.latitudeStepDeg >= 85 ||
		    grid.firstLongitudeDeg <= -180 || grid.firstLongitudeDeg + nx * grid.longitudeStepDeg >= 180)
			throw std::invalid_argument("Map extent crosses a pole or the dateline");
		grid.count = density(count, grid, resolved, cancel);
		grid.mass = density(mass, grid, resolved, cancel);
		return grid;
	}

}  // namespace strewn
