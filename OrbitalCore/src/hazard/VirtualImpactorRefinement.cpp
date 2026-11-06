#include "hazard/VirtualImpactorRefinement.h"

#include "hazard/ObservationNormalMetric.h"
#include "orbital_core/Constants.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>

namespace orbital_core::hazard {
	namespace {
		[[nodiscard]] bool normalizeAbove(Eigen::Ref<Eigen::VectorXd> vector, double tolerance) noexcept
		{
			const double norm = vector.norm();
			if (!std::isfinite(norm) || !(norm > tolerance))
				return false;
			vector /= norm;
			return vector.array().isFinite().all();
		}

		[[nodiscard]] bool buildTargetPlaneAdaptedBasis(const Eigen::MatrixXd& dCoordsDParams, Eigen::MatrixXd& basis,
		                                                const char*& failureStage) noexcept
		{
			const Eigen::Index n = dCoordsDParams.cols();
			if (dCoordsDParams.rows() != 2 || n < 2) {
				failureStage = "adapted-basis-shape";
				return false;
			}

			const Eigen::VectorXd first = dCoordsDParams.row(0).transpose();
			const Eigen::VectorXd second = dCoordsDParams.row(1).transpose();
			const double firstNorm = first.norm();
			const double secondNorm = second.norm();
			const double tolerance = 1.0e-12 * std::min(firstNorm, secondNorm);
			if (!std::isfinite(tolerance) || tolerance == 0.0) {
				failureStage = "adapted-basis-rank-zero";
				return false;
			}

			basis = Eigen::MatrixXd::Zero(n, n);
			basis.col(0) = first;
			if (!normalizeAbove(basis.col(0), tolerance)) {
				failureStage = "adapted-basis-first-small";
				return false;
			}

			Eigen::VectorXd secondUnit = second;
			if (!normalizeAbove(secondUnit, tolerance)) {
				failureStage = "adapted-basis-second-small";
				return false;
			}
			basis.col(1) = secondUnit - basis.col(0).dot(secondUnit) * basis.col(0);
			if (!normalizeAbove(basis.col(1), tolerance)) {
				failureStage = "adapted-basis-rank-one";
				return false;
			}

			Eigen::Index accepted = 0;
			for (Eigen::Index j = 0; j < n && accepted < n - 2; ++j) {
				Eigen::VectorXd candidate =
				    Eigen::VectorXd::Unit(n, j) - basis(j, 0) * basis.col(0) - basis(j, 1) * basis.col(1);
				if (!normalizeAbove(candidate, tolerance))
					continue;

				for (Eigen::Index k = 0; k < accepted; ++k) {
					candidate -= candidate.dot(basis.col(2 + k)) * basis.col(2 + k);
				}
				if (!normalizeAbove(candidate, tolerance))
					continue;

				basis.col(2 + accepted) = candidate;
				++accepted;
			}

			if (accepted != n - 2 || !basis.array().isFinite().all()) {
				failureStage = "adapted-basis-completion";
				return false;
			}
			return true;
		}

		[[nodiscard]] bool semilinearCorrection(const Eigen::MatrixXd& dCoordsDParams,
		                                        const FitParameterCovariance& parameterCovariance,
		                                        const Eigen::Vector2d& targetDelta, Eigen::VectorXd& correction,
		                                        const char*& failureStage) noexcept
		{
			const Eigen::Index n = dCoordsDParams.cols();
			Eigen::MatrixXd basis;
			if (!buildTargetPlaneAdaptedBasis(dCoordsDParams, basis, failureStage)) {
				return false;
			}

			const FitParameterCovariance covariance = 0.5 * (parameterCovariance + parameterCovariance.transpose());
			const double covarianceScale = covariance.cwiseAbs().maxCoeff();
			if (!(covarianceScale > 0.0) || !std::isfinite(covarianceScale)) {
				failureStage = "covariance-scale";
				return false;
			}

			const FitParameterCovariance scaledCovariance = covariance / covarianceScale;
			const Eigen::MatrixXd covarianceInBasis = basis.transpose() * scaledCovariance * basis;
			const Eigen::MatrixXd c4 = covarianceInBasis.bottomRightCorner(n - 2, n - 2);
			const Eigen::MatrixXd c42 = covarianceInBasis.bottomLeftCorner(n - 2, 2);
			Eigen::LLT<Eigen::MatrixXd> c4Factor(0.5 * (c4 + c4.transpose()));
			if (c4Factor.info() != Eigen::Success) {
				failureStage = "c4-factor";
				return false;
			}

			const Eigen::MatrixXd ceicel = c4Factor.solve(c42);
			if (c4Factor.info() != Eigen::Success || !ceicel.array().isFinite().all()) {
				failureStage = "c4-solve";
				return false;
			}

			const Eigen::Matrix2d a = dCoordsDParams * basis.leftCols(2);
			const double determinant = a.determinant();
			if (determinant == 0.0 || !std::isfinite(determinant)) {
				failureStage = "target-map";
				return false;
			}

			Eigen::Matrix2d b;
			b << a(1, 1), -a(0, 1), -a(1, 0), a(0, 0);
			b /= determinant;

			const Eigen::Vector2d del2 = b * targetDelta;
			Eigen::VectorXd adaptedCorrection(n);
			adaptedCorrection.head(2) = del2;
			adaptedCorrection.tail(n - 2) = -ceicel * del2;
			correction = basis * adaptedCorrection;
			if (correction.size() != n || !correction.array().isFinite().all()) {
				failureStage = "correction";
				return false;
			}
			return true;
		}
	}  // namespace

	Eigen::Vector2d impactBoundaryTarget(const Eigen::Vector2d& coordsKm, double impactRadiusKm) noexcept
	{
		const double radius = std::max(0.0, impactRadiusKm);
		if (radius == 0.0)
			return Eigen::Vector2d::Zero();

		const double norm = coordsKm.norm();
		if (norm > 1.0e-12)
			return coordsKm * (radius / norm);

		return Eigen::Vector2d(radius, 0.0);
	}

	VirtualImpactorRefinementResult refineVirtualImpactorWithObservations(const VirtualImpactorRefinementInput& input)
	{
		VirtualImpactorRefinementResult out;
		out.parameters = input.initialParameters;

		const Eigen::Index n = input.initialParameters.size();
		if (!input.cfg || !input.cfg->rmsScaling.linearization || n < 2 || input.dCoordsDParams.rows() != 2 ||
		    input.dCoordsDParams.cols() != n || input.parameterCovariance.rows() != n ||
		    input.parameterCovariance.cols() != n || !input.initialParameters.array().isFinite().all() ||
		    !input.currentCoordsKm.array().isFinite().all() || !input.targetCoordsKm.array().isFinite().all() ||
		    !input.dCoordsDParams.array().isFinite().all() || !input.parameterCovariance.array().isFinite().all()) {
			return out;
		}

		out.attempted = true;

		constexpr double kmPerAu = Constants::AU_IN_METERS / 1000.0;
		const Eigen::MatrixXd dCoordsDParamsAu = input.dCoordsDParams / kmPerAu;

		const Eigen::Vector2d targetDelta = (input.targetCoordsKm - input.currentCoordsKm) / kmPerAu;
		Eigen::VectorXd correction;
		if (!semilinearCorrection(dCoordsDParamsAu, input.parameterCovariance, targetDelta, correction,
		                          out.failureStage)) {
			out.badConditioning = true;
			return out;
		}

		out.parameters = input.initialParameters + correction;
		if (!out.parameters.array().isFinite().all()) {
			out.failureStage = "parameters";
			return out;
		}

		const auto lin = input.cfg->rmsScaling.linearization(out.parameters);
		const double objectiveScale =
		    (std::isfinite(input.cfg->rmsScaling.objectiveScale) && input.cfg->rmsScaling.objectiveScale > 0.0)
		        ? input.cfg->rmsScaling.objectiveScale
		        : 1.0;
		if (!lin.valid || lin.jacobian.cols() != n ||
		    lin.jacobian.rows() != static_cast<Eigen::Index>(lin.residuals.size()) || !std::isfinite(lin.sumSquares)) {
			out.failureStage = "observation-linearization";
			return out;
		}

		const ObservationNormalMetric metric = makeObservationNormalMetric(lin.jacobian, lin.sumSquares);
		if (!metric.valid) {
			out.badConditioning = true;
			out.failureStage = "observation-normal";
			return out;
		}

		out.observationChiSquare = lin.sumSquares / objectiveScale;
		out.parameterCovariance = metric.covariance;
		out.valid = std::isfinite(out.observationChiSquare);
		return out;
	}

}  // namespace orbital_core::hazard
