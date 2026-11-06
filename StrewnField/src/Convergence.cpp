#include "Convergence.h"
#include "VectorMath.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace strewn::detail {

	void validateConvergence(const ConvergenceConfig& c)
	{
		if (!std::isfinite(c.threshold) || c.threshold <= 0 || c.threshold > 1 || c.batchSize == 0 ||
		    c.consecutiveChecks == 0 || !std::isfinite(c.centerDeg[0]) || std::abs(c.centerDeg[0]) >= 89 ||
		    !std::isfinite(c.centerDeg[1]) || std::abs(c.centerDeg[1]) > 180 || !std::isfinite(c.cellM) ||
		    c.cellM <= 0 || !std::isfinite(c.kernelSigmaM) || c.kernelSigmaM <= 0 || c.kernelSigmaM / c.cellM > 64)
			throw std::invalid_argument("Invalid JS convergence parameters");
		double cells = 1;
		for (const auto half : c.halfSizeKm) {
			const double n = std::floor(2000 * half / c.cellM);
			if (!std::isfinite(n) || n < 2 || n > 4000000)
				throw std::invalid_argument("Invalid JS grid extent");
			cells *= n;
		}
		if (cells > 4000000)
			throw std::invalid_argument("JS grid exceeds 4000000 cells");
	}

	double jsDivergence(std::span<const double> p, std::span<const double> q)
	{
		if (p.size() != q.size())
			throw std::invalid_argument("JS map dimensions differ");
		double sum = 0;
		for (std::size_t i = 0; i < p.size(); ++i) {
			const double mean = 0.5 * (p[i] + q[i]);
			if (p[i] > 0)
				sum += p[i] * std::log2(p[i] / mean);
			if (q[i] > 0)
				sum += q[i] * std::log2(q[i] / mean);
		}
		return std::clamp(0.5 * sum, 0.0, 1.0);
	}

	ConvergenceMonitor::ConvergenceMonitor(ConvergenceConfig config) : m_config(config)
	{
		validateConvergence(config);
		m_ny = static_cast<std::size_t>(2000 * config.halfSizeKm[0] / config.cellM);
		m_nx = static_cast<std::size_t>(2000 * config.halfSizeKm[1] / config.cellM);
		m_xmin = m_nx;
		m_ymin = m_ny;
		m_dy = config.cellM / (6367e3 * deg);
		m_dx = m_dy / std::cos(config.centerDeg[0] * deg);
		for (auto& h : m_histogram)
			h.resize(m_nx * m_ny);
		const double sigma = config.kernelSigmaM / config.cellM;
		const auto radius = static_cast<std::size_t>(std::ceil(8 * sigma));
		m_kernel.resize(radius + 1);
		for (std::size_t i = 0; i <= radius; ++i)
			m_kernel[i] = std::exp(-0.5 * std::pow(static_cast<double>(i) / sigma, 2));
	}

	std::vector<double> ConvergenceMonitor::density(const std::vector<double>& h) const
	{
		std::vector<double> scratch(h.size()), result(h.size());
		const auto radius = m_kernel.size() - 1;
		const auto x0 = m_xmin > radius ? m_xmin - radius : 0;
		const auto x1 = std::min(m_xmax + radius + 1, m_nx);
		const auto y0 = m_ymin > radius ? m_ymin - radius : 0;
		const auto y1 = std::min(m_ymax + radius + 1, m_ny);
		for (std::size_t y = m_ymin; y <= m_ymax; ++y) {
			for (std::size_t x = m_xmin; x <= m_xmax; ++x) {
				const double weight = h[y * m_nx + x];
				if (weight == 0)
					continue;
				scratch[y * m_nx + x] += weight;
				for (std::size_t d = 1; d <= radius; ++d) {
					if (d <= x)
						scratch[y * m_nx + x - d] += weight * m_kernel[d];
					if (x + d < m_nx)
						scratch[y * m_nx + x + d] += weight * m_kernel[d];
				}
			}
		}
		for (std::size_t y = y0; y < y1; ++y) {
			const auto first = y > radius ? y - radius : 0;
			const auto last = std::min(y + radius + 1, m_ny);
			for (std::size_t row = first; row < last; ++row) {
				const double kernel = m_kernel[row > y ? row - y : y - row];
				for (std::size_t x = x0; x < x1; ++x)
					result[y * m_nx + x] += kernel * scratch[row * m_nx + x];
			}
		}
		const double total = std::accumulate(result.begin(), result.end(), 0.0);
		if (!std::isfinite(total) || total <= 0)
			throw std::runtime_error("Invalid total weight in JS density");
		for (double& value : result)
			value /= total;
		return result;
	}

	void ConvergenceMonitor::update(const RunResult& run, RunSummary& summary)
	{
		summary.convergenceBlocked |= run.status != RunStatus::Complete;
		for (const auto& impact : run.impacts) {
			const auto& f = impact.fragment;
			const double x = (impact.longitudeDeg - m_config.centerDeg[1]) / m_dx + m_nx / 2.0 + 0.5;
			const double y = (impact.latitudeDeg - m_config.centerDeg[0]) / m_dy + m_ny / 2.0 + 0.5;
			const double mass = f.multiplicity * f.state.massKg;
			if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(mass) || mass <= 0 ||
			    !std::isfinite(f.multiplicity) || f.multiplicity <= 0)
				throw std::runtime_error("Invalid impact in JS convergence monitor");
			if (x < 0 || y < 0 || x > m_nx || y > m_ny) {
				++summary.convergenceOutOfGridImpacts;
				summary.convergenceBlocked = true;
				continue;
			}
			const auto ix = std::min(static_cast<std::size_t>(x), m_nx - 1);
			const auto iy = std::min(static_cast<std::size_t>(y), m_ny - 1);
			m_histogram[0][iy * m_nx + ix] += f.multiplicity;
			m_histogram[1][iy * m_nx + ix] += mass;
			m_xmin = std::min(m_xmin, ix);
			m_xmax = std::max(m_xmax, ix);
			m_ymin = std::min(m_ymin, iy);
			m_ymax = std::max(m_ymax, iy);
		}
		if (summary.completed % m_config.batchSize != 0)
			return;
		ConvergenceCheck check;
		check.runs = summary.completed;
		if (!summary.convergenceBlocked && m_xmin < m_nx) {
			std::array<std::vector<double>, 2> maps{density(m_histogram[0]), density(m_histogram[1])};
			if (!m_previous[0].empty()) {
				check.jsCount = jsDivergence(m_previous[0], maps[0]);
				check.jsMass = jsDivergence(m_previous[1], maps[1]);
				m_stableChecks = std::max(*check.jsCount, *check.jsMass) < m_config.threshold ? m_stableChecks + 1 : 0;
			}
			m_previous = std::move(maps);
		} else {
			m_stableChecks = 0;
		}
		check.stableChecks = m_stableChecks;
		summary.convergenceHistory.push_back(check);
		summary.converged = m_stableChecks >= m_config.consecutiveChecks;
	}

}  // namespace strewn::detail
