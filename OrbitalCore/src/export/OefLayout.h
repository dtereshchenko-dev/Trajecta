#pragma once

#include "orbital_core/nongrav/NonGravModel.h"
#include <Eigen/Core>
#include <stdexcept>
#include <vector>

namespace orbital_core::oef {

	// OEF model slots are not NgModel values. Yarkovsky uses slot 2 and
	// 1e-10 AU/day^2; comet coefficients use physical AU/day^2, DT uses days.
	struct Layout {
		int model{};
		std::vector<double> values;
		std::vector<int> solved;
		Eigen::VectorXd scales;
	};

	// Reserved Core-only layouts must never be advertised as standard OEF2.0.
	inline bool requiresExtension(const NonGravModel& ng)
	{
		return ng.enabled() && ng.model() != NgModel::None && ng.model() != NgModel::Yarkovsky &&
		       ng.model() != NgModel::MarsdenH2O;
	}

	inline Layout layout(const NonGravModel& ng)
	{
		Layout out;
		out.scales = Eigen::VectorXd::Ones(6 + ng.numFitParams());
		if (!ng.enabled() || ng.model() == NgModel::None)
			return out;
		switch (ng.model()) {
			case NgModel::Yarkovsky:
				out.model = 1;
				out.values = {0.0, ng.A2() * 1e10};
				for (int i = 0; i < ng.numFitParams(); ++i) {
					out.solved.push_back(2);
					out.scales(6 + i) = 1e10;
				}
				break;
			case NgModel::MarsdenH2O:
				out.model = 2;
				out.values = {ng.A1(), ng.A2(), ng.A3(), ng.DT()};
				for (int i = 0; i < ng.numFitParams(); ++i)
					out.solved.push_back(1 + static_cast<int>(ng.fitParameterAt(i)));
				break;
			case NgModel::MarsdenCO:
			case NgModel::InverseSquare:
				out.model = ng.model() == NgModel::MarsdenCO ? 101 : 103;
				out.values = {ng.A1(), ng.A2(), ng.A3(), ng.DT()};
				for (int i = 0; i < ng.numFitParams(); ++i)
					out.solved.push_back(1 + static_cast<int>(ng.fitParameterAt(i)));
				break;
			default:
				throw std::invalid_argument("Unknown non-gravitational model");
		}
		return out;
	}

	inline Eigen::MatrixXd scaledCovariance(const Eigen::MatrixXd& covariance, const NonGravModel& ng)
	{
		const auto spec = layout(ng);
		if (covariance.rows() != spec.scales.size() || covariance.cols() != spec.scales.size())
			throw std::invalid_argument("OEF covariance does not include the complete fitted parameter layout");
		return spec.scales.asDiagonal() * covariance * spec.scales.asDiagonal();
	}

}  // namespace orbital_core::oef
