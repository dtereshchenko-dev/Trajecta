#pragma once

#include "orbital_core/FitResult.h"

namespace orbital_core::uncertainty {

	enum class CovarianceBasis { Cartesian, Equinoctial };

	struct CovarianceSpectrum {
		Eigen::VectorXd signedSigmas;
		Eigen::VectorXd weakDirection;
	};

	struct CovarianceDiagnostics {
		FitParameterCovariance covariance;
		double epochMjdTdb{};
		double objectiveScale{1.0};
		std::optional<CovarianceSpectrum> spectrum;
	};

	[[nodiscard]] ORBITAL_CORE_API std::optional<CovarianceSpectrum> covarianceSpectrum(
	    const FitParameterCovariance& covariance);
	[[nodiscard]] ORBITAL_CORE_API std::optional<CovarianceDiagnostics> describeCovariance(const FitResult& result,
	                                                                                       CovarianceBasis basis);

}  // namespace orbital_core::uncertainty
