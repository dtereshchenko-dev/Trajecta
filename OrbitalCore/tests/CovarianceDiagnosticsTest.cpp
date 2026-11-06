#include <gtest/gtest.h>
#include "orbital_core/uncertainty/CovarianceDiagnostics.h"
#include "orbital_core/StateVectors.h"
#ifdef ORBITAL_CORE_HAZARD_AVAILABLE
#include "hazard/LovCoordinateChart.h"
#endif
#include <limits>

using namespace orbital_core;
using namespace orbital_core::uncertainty;

TEST(CovarianceDiagnostics, SpectrumRetainsExportConvention)
{
	FitParameterCovariance covariance(2, 2);
	covariance << 2, -1, -1, 2;
	const auto spectrum = covarianceSpectrum(covariance);
	ASSERT_TRUE(spectrum);
	EXPECT_NEAR(spectrum->signedSigmas(0), 1.0, 1e-14);
	EXPECT_NEAR(spectrum->signedSigmas(1), std::sqrt(3.0), 1e-14);
	EXPECT_NEAR(spectrum->weakDirection(0), std::sqrt(0.5), 1e-14);
	EXPECT_NEAR(spectrum->weakDirection(1), -std::sqrt(0.5), 1e-14);
	covariance << -4, 0, 0, 9;
	const auto indefinite = covarianceSpectrum(covariance);
	ASSERT_TRUE(indefinite);
	EXPECT_DOUBLE_EQ(indefinite->signedSigmas(0), -2.0);
	EXPECT_DOUBLE_EQ(indefinite->signedSigmas(1), 3.0);
	EXPECT_FALSE(covarianceSpectrum(FitParameterCovariance{}));
	EXPECT_FALSE(covarianceSpectrum(FitParameterCovariance::Zero(6, 6)));
	EXPECT_FALSE(covarianceSpectrum(FitParameterCovariance::Ones(2, 3)));
	covariance(0, 0) = std::numeric_limits<double>::quiet_NaN();
	EXPECT_FALSE(covarianceSpectrum(covariance));
}

TEST(CovarianceDiagnostics, FitEpochScaleAndNonGravCrossTermsStayTogether)
{
	FitResult fit;
	fit.fitEpochOrbit =
	    make_orbit_from_states_vectors(StateVectors(Vec3D(1, 0.2, 0.1), Vec3D(-0.1, 0.8, 0.1)), 2460000.5);
	ASSERT_TRUE(fit.fitEpochOrbit);
	fit.orbit = make_orbit_from_states_vectors(StateVectors(Vec3D(1, 0.2, 0.1), Vec3D(-0.1, 0.8, 0.1)), 2460030.5);
	fit.finalNonGrav.setEnabled(true);
	fit.finalNonGrav.setModel(NgModel::Yarkovsky);
	auto& covariance = fit.uncertainty.fitParameterCovariance;
	covariance = FitParameterCovariance::Identity(7, 7) * 1e-8;
	covariance(0, 6) = covariance(6, 0) = 2e-10;
	fit.uncertainty.covarianceObjectiveScale = 1.04;
	const auto car = describeCovariance(fit, CovarianceBasis::Cartesian);
	ASSERT_TRUE(car);
	EXPECT_DOUBLE_EQ(car->epochMjdTdb, 60000.0);
	EXPECT_DOUBLE_EQ(car->objectiveScale, 1.04);
	EXPECT_TRUE(car->covariance.isApprox(covariance, 0.0));
#ifdef ORBITAL_CORE_HAZARD_AVAILABLE
	const auto equ = describeCovariance(fit, CovarianceBasis::Equinoctial);
	ASSERT_TRUE(equ);
	const auto chart = hazard::makeEquinoctialLovChart(*fit.fitEpochOrbit, covariance, fit.finalNonGrav);
	ASSERT_TRUE(chart.valid);
	EXPECT_TRUE(equ->covariance.isApprox(chart.coordinateCovariance, 0.0));
	EXPECT_DOUBLE_EQ(equ->epochMjdTdb, car->epochMjdTdb);
	EXPECT_EQ(equ->covariance.rows(), 7);
	EXPECT_GT((equ->covariance.block<6, 1>(0, 6).norm()), 0.0);
	EXPECT_DOUBLE_EQ(equ->covariance(6, 6), covariance(6, 6));
#endif
	fit.uncertainty.fitParameterCovariance.resize(6, 6);
	EXPECT_FALSE(describeCovariance(fit, CovarianceBasis::Cartesian));
	EXPECT_FALSE(describeCovariance({}, CovarianceBasis::Cartesian));
}

TEST(CovarianceDiagnostics, CartesianFallbackDoesNotInventEquinoctialOrNonGravData)
{
	FitResult fit;
	fit.fitEpochOrbit = make_orbit_from_states_vectors(StateVectors(Vec3D(1, 0, 0), Vec3D(0, 2, 0)), 2460000.5);
	ASSERT_TRUE(fit.fitEpochOrbit);
	fit.uncertainty.cartesianStateCovariance = StateCov6x6::Identity() * 1e-8;
	fit.uncertainty.cartesianStateCovarianceFormal = StateCov6x6::Identity() * 1e-8;
	ASSERT_TRUE(describeCovariance(fit, CovarianceBasis::Cartesian));
	EXPECT_FALSE(describeCovariance(fit, CovarianceBasis::Equinoctial));
	fit.finalNonGrav.setModel(NgModel::Yarkovsky);
	fit.finalNonGrav.setEnabled(true);
	EXPECT_FALSE(describeCovariance(fit, CovarianceBasis::Cartesian));
}
