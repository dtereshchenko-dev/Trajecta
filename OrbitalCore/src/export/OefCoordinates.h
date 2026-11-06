#pragma once
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/Types.h"
#include <stdexcept>

namespace orbital_core::oef {

	// Analytic KEP <-> COT derivative at a fixed epoch. This export/import boundary
	// does not use the finite-difference LOV chart and does not change LOV sampling.
	inline cometary::Matrix6 keplerianToCometaryJacobian(const cometary::Vector6& cot)
	{
		const double e = cot(1), nu = cot(5);
		if (!(e >= 0 && e < 1 && cot(0) > 0))
			throw std::invalid_argument("KEP covariance requires a nonsingular elliptic orbit");
		const double beta2 = 1 - e * e;
		const double w = 1 + e * std::cos(nu);
		auto jac = cometary::Matrix6::Identity().eval();
		jac(0, 0) = 1 - e;
		jac(0, 1) = -cot(0) / (1 - e);
		jac(5, 1) = std::sin(nu) * (2 + e * std::cos(nu)) / beta2;
		jac(5, 5) = w * w / std::pow(beta2, 1.5);
		return jac;
	}

	inline cometary::Vector6 cotFromKeplerian(const cometary::Vector6& kep)
	{
		const double e = kep(1);
		if (!(kep(0) > 0 && e >= 0 && e < 1))
			throw std::invalid_argument("Unsupported KEP orbit");
		const double E = Calculations::kepler_solve_e(Helpers::wrapTwoPi(kep(5)), e);
		auto cot = kep;
		cot(0) *= 1 - e;
		cot(5) = std::atan2(std::sqrt(1 - e * e) * std::sin(E), std::cos(E) - e);
		return cot;
	}

	inline Eigen::MatrixXd keplerianToPhysicalJacobian(const Eigen::VectorXd& kep)
	{
		const auto cot = cotFromKeplerian(kep.head<6>());
		cometary::Matrix6 carCot;
		if (!cometary::toCartesian(cot, &carCot))
			throw std::invalid_argument("Singular KEP covariance chart");
		Eigen::MatrixXd jac = Eigen::MatrixXd::Identity(kep.size(), kep.size());
		jac.topLeftCorner<6, 6>() = carCot * keplerianToCometaryJacobian(cot);
		return jac;
	}

	inline Eigen::MatrixXd physicalToKeplerianJacobian(const Eigen::VectorXd& kep)
	{
		const auto cot = cotFromKeplerian(kep.head<6>());
		const auto cotCar = cometary::cartesianToElementsJacobian(cot);
		if (!cotCar)
			throw std::invalid_argument("Singular KEP covariance chart");
		const auto forward = keplerianToCometaryJacobian(cot);
		auto inverse = cometary::Matrix6::Identity().eval();
		inverse(0, 0) = 1.0 / forward(0, 0);
		inverse(0, 1) = -forward(0, 1) / forward(0, 0);
		inverse(5, 5) = 1.0 / forward(5, 5);
		inverse(5, 1) = -forward(5, 1) / forward(5, 5);
		Eigen::MatrixXd jac = Eigen::MatrixXd::Identity(kep.size(), kep.size());
		jac.topLeftCorner<6, 6>() = inverse * *cotCar;
		return jac;
	}
}  // namespace orbital_core::oef
