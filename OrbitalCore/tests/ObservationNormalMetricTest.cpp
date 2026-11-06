#include "hazard/ObservationNormalMetric.h"

#include <gtest/gtest.h>
#include <cmath>

namespace orbital_core::hazard {
	namespace {

		TEST(ObservationNormalMetricTest, PreservesResolvedShortArcWeakDirection)
		{
			Eigen::MatrixXd rotation = Eigen::MatrixXd::Identity(6, 6);
			const double c = std::sqrt(0.5);
			rotation(0, 0) = rotation(5, 5) = c;
			rotation(0, 5) = -c;
			rotation(5, 0) = c;
			Eigen::VectorXd singular(6);
			singular << 1.7e9, 4.1e8, 4.7e5, 1.47e5, 355.0, 24.0;
			Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(14, 6);
			jacobian.topRows(6) = singular.asDiagonal() * rotation.transpose();

			const auto metric = makeObservationNormalMetric(jacobian, 0.1);

			ASSERT_TRUE(metric.valid);
			EXPECT_NEAR(metric.eigenvalues(0), 576.0, 1.0e-4);
			EXPECT_LT(metric.eigenvalueFloor, metric.eigenvalues(0));
			EXPECT_NEAR(std::abs(metric.eigenvectors.col(0).dot(rotation.col(5))), 1.0, 1.0e-10);
			const double variance = rotation.col(5).dot(metric.covariance * rotation.col(5));
			EXPECT_NEAR(variance, 1.75 / 576.0, 1.0e-10);
		}

		TEST(ObservationNormalMetricTest, RejectsUnresolvedRankInsteadOfZeroVariance)
		{
			Eigen::MatrixXd jacobian = Eigen::MatrixXd::Identity(14, 6);
			jacobian.col(5).setZero();
			EXPECT_FALSE(makeObservationNormalMetric(jacobian, 0.1).valid);
		}

		TEST(ObservationNormalMetricTest, UsesAPosterioriScaleForLowChiSquare)
		{
			Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(156, 6);
			for (int i = 0; i < 6; ++i)
				jacobian(i, i) = static_cast<double>(i + 1);

			const auto metric = makeObservationNormalMetric(jacobian, 50.0);

			ASSERT_TRUE(metric.valid);
			EXPECT_DOUBLE_EQ(metric.covarianceScale, 1.04);
			EXPECT_TRUE((jacobian.transpose() * jacobian * metric.covariance)
			                .isApprox(1.04 * Eigen::MatrixXd::Identity(6, 6), 1.0e-12));
		}

		TEST(ObservationNormalMetricTest, UsesResidualRmsForHighChiSquare)
		{
			const Eigen::MatrixXd jacobian = Eigen::MatrixXd::Identity(100, 6);

			const auto metric = makeObservationNormalMetric(jacobian, 150.0);

			ASSERT_TRUE(metric.valid);
			EXPECT_DOUBLE_EQ(metric.covarianceScale, 150.0 / 94.0);
			EXPECT_TRUE(metric.covariance.isApprox((150.0 / 94.0) * Eigen::MatrixXd::Identity(6, 6), 1.0e-12));
		}

		TEST(ObservationNormalMetricTest, UsesResidualCountAndDegreesOfFreedom)
		{
			Eigen::MatrixXd shortJacobian = Eigen::MatrixXd::Zero(8, 6);
			shortJacobian.topRows(6) = Eigen::MatrixXd::Identity(6, 6);
			Eigen::MatrixXd longJacobian = Eigen::MatrixXd::Zero(156, 6);
			longJacobian.topRows(8) = shortJacobian;

			const auto shortMetric = makeObservationNormalMetric(shortJacobian, 10.0);
			const auto longMetric = makeObservationNormalMetric(longJacobian, 10.0);

			ASSERT_TRUE(shortMetric.valid);
			ASSERT_TRUE(longMetric.valid);
			EXPECT_DOUBLE_EQ(shortMetric.covarianceScale, 5.0);
			EXPECT_DOUBLE_EQ(longMetric.covarianceScale, 1.04);
			EXPECT_TRUE((shortMetric.covariance / shortMetric.covarianceScale)
			                .isApprox(longMetric.covariance / longMetric.covarianceScale, 1.0e-12));
		}

	}  // namespace
}  // namespace orbital_core::hazard
