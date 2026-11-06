#pragma once

#include "NBodyPropagator.h"
#include "obs/PreparedObservation.h"
#include <Eigen/Core>
#include <functional>
#include <span>

namespace orbital_core::hazard::mov {

	using Vector6 = Eigen::Matrix<double, 6, 1>;
	using Matrix6 = Eigen::Matrix<double, 6, 6>;

	struct Attributable {
		double epoch{}, curvatureSnr{};
		Eigen::Vector4d angles;
		obs::ObserverData observer;
	};

	Attributable fitAttributable(std::span<const obs::PreparedObservation> observations);

	// ATT is equatorial; Cartesian state/derivative are heliocentric ecliptic
	// AU and AU/day. The reference observer is fixed during the four-dimensional
	// fit, but depends on rho through reception time when differentiating rho.
	struct Chart {
		double epoch{};
		std::function<NBodyPropagator::StateVel(double)> observer;
		Vector6 state(const Vector6& att) const;
		Matrix6 derivative(const Vector6& att) const;
		Vector6 inverse(const Vector6& state) const;
	};

	// Heliocentric bound-orbit condition a <= 100 AU, using the preliminary
	// attributable. Connected range components come from the sixth-degree
	// energy polynomial, with spurious squared roots rejected against energy.
	struct AdmissibleRegion {
		Eigen::Vector3d position, velocity, direction, angularVelocity;
		std::vector<std::pair<double, double>> ranges;
		double rateMin{}, rateMax{};
		double discriminant(double range) const;
		bool contains(double range, double rate) const;
	};

	AdmissibleRegion admissibleRegion(const Chart& chart, const Eigen::Vector4d& angles);

	// Formal observation metric, without posterior covariance inflation. The
	// tangent dA/d(rho,rhodot) = -(B_A^T B_A)^-1 B_A^T B_rho is coordinate-dependent;
	// units are the ATT units above. Sampling Jacobians are applied separately.
	double manifoldJacobian(const Eigen::MatrixXd& weightedDesign);

}  // namespace orbital_core::hazard::mov
