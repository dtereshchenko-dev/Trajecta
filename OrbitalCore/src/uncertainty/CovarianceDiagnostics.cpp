#include "orbital_core/uncertainty/CovarianceDiagnostics.h"
#ifdef ORBITAL_CORE_HAZARD_AVAILABLE
#include "hazard/LovCoordinateChart.h"
#endif
#include <Eigen/Eigenvalues>
#include <cmath>

namespace orbital_core::uncertainty {

	std::optional<CovarianceSpectrum> covarianceSpectrum(const FitParameterCovariance& covariance)
	{
		if (covariance.rows() == 0 || covariance.rows() != covariance.cols() || !covariance.allFinite() ||
		    covariance.isZero(0.0))
			return std::nullopt;
		const Eigen::MatrixXd symmetric = covariance.selfadjointView<Eigen::Upper>();
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(symmetric);
		if (solver.info() != Eigen::Success)
			return std::nullopt;
		CovarianceSpectrum result;
		result.signedSigmas = solver.eigenvalues().unaryExpr(
		    [](double value) { return std::copysign(std::sqrt(std::abs(value)), value); });
		result.weakDirection = solver.eigenvectors().col(solver.eigenvalues().size() - 1);
		if (result.weakDirection(0) < 0.0)
			result.weakDirection = -result.weakDirection;
		return result;
	}

	std::optional<CovarianceDiagnostics> describeCovariance(const FitResult& result, CovarianceBasis basis)
	{
		if (!result.fitEpochOrbit)
			return std::nullopt;
		const int n = 6 + result.finalNonGrav.numFitParams();
		const auto& uncertainty = result.uncertainty;
		CovarianceDiagnostics out;
		const auto& full = uncertainty.fitParameterCovariance;
		if (full.rows() == n && full.cols() == n && full.allFinite() && !full.isZero(0.0))
			out.covariance = full;
		else if (n == 6 && uncertainty.hasStateCovariance())
			out.covariance = uncertainty.cartesianStateCovariance;
		else
			return std::nullopt;
		if (basis == CovarianceBasis::Equinoctial) {
#ifdef ORBITAL_CORE_HAZARD_AVAILABLE
			auto chart = hazard::makeEquinoctialLovChart(*result.fitEpochOrbit, out.covariance, result.finalNonGrav);
			if (!chart.valid)
				return std::nullopt;
			out.covariance = std::move(chart.coordinateCovariance);
#else
			return std::nullopt;
#endif
		}
		out.epochMjdTdb = result.fitEpochOrbit->epochMjd();
		out.objectiveScale = uncertainty.covarianceObjectiveScale;
		out.spectrum = covarianceSpectrum(out.covariance);
		return out;
	}

}  // namespace orbital_core::uncertainty
