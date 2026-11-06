#pragma once

#include "strewn_core/Types.h"
#include <span>

namespace strewn::detail {

	void validateConvergence(const ConvergenceConfig& config);
	double jsDivergence(std::span<const double> p, std::span<const double> q);

	class ConvergenceMonitor {
	   public:
		explicit ConvergenceMonitor(ConvergenceConfig config);
		void update(const RunResult& result, RunSummary& summary);

	   private:
		std::vector<double> density(const std::vector<double>& histogram) const;

		ConvergenceConfig m_config;
		std::size_t m_nx{}, m_ny{};
		std::size_t m_xmin{}, m_xmax{}, m_ymin{}, m_ymax{};
		double m_dx{}, m_dy{};
		std::vector<double> m_kernel;
		std::array<std::vector<double>, 2> m_histogram, m_previous;
		std::size_t m_stableChecks{};
	};

}  // namespace strewn::detail
