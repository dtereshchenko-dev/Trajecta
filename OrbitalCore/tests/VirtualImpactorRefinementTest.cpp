#include "hazard/VirtualImpactorRefinement.h"

#include <gtest/gtest.h>

namespace orbital_core::hazard {
	namespace {

		ImpactRiskConfig refinementConfig()
		{
			ImpactRiskConfig cfg;
			cfg.rmsScaling.linearization = [](const Eigen::VectorXd& x) {
				ImpactRiskRmsScaling::ObservationLinearization out;
				out.residuals.assign(8, 0.0);
				out.jacobian = Eigen::MatrixXd::Zero(8, x.size());
				out.jacobian.topRows(x.size()) = Eigen::MatrixXd::Identity(x.size(), x.size());
				out.sumSquares = 42.0;
				out.valid = true;
				return out;
			};
			return cfg;
		}

		TEST(VirtualImpactorRefinementTest, SemilinearStepSatisfiesTargetPlaneConstraint)
		{
			ImpactRiskConfig cfg = refinementConfig();

			VirtualImpactorRefinementInput input;
			input.cfg = &cfg;
			input.initialParameters = Eigen::VectorXd::Zero(6);
			input.dCoordsDParams = Eigen::MatrixXd::Zero(2, 6);
			input.dCoordsDParams(0, 0) = 1.0;
			input.dCoordsDParams(1, 1) = 1.0;
			input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6);
			input.parameterCovariance(0, 0) = 4.0;
			input.parameterCovariance(1, 1) = 9.0;
			input.currentCoordsKm = Eigen::Vector2d(10.0, 20.0);
			input.targetCoordsKm = Eigen::Vector2d(2.0, 3.0);

			const VirtualImpactorRefinementResult result = refineVirtualImpactorWithObservations(input);

			ASSERT_TRUE(result.attempted);
			ASSERT_TRUE(result.valid);
			EXPECT_FALSE(result.badConditioning);
			EXPECT_NEAR(result.parameters(0), -8.0, 1.0e-12);
			EXPECT_NEAR(result.parameters(1), -17.0, 1.0e-12);
			const Eigen::Vector2d reached =
			    input.currentCoordsKm + input.dCoordsDParams * (result.parameters - input.initialParameters);
			EXPECT_TRUE(reached.isApprox(input.targetCoordsKm, 1.0e-12));
			EXPECT_DOUBLE_EQ(result.observationChiSquare, 42.0);
		}

		TEST(VirtualImpactorRefinementTest, RejectsRankOneTargetPlaneMap)
		{
			ImpactRiskConfig cfg = refinementConfig();

			VirtualImpactorRefinementInput input;
			input.cfg = &cfg;
			input.initialParameters = Eigen::VectorXd::Zero(6);
			input.dCoordsDParams = Eigen::MatrixXd::Zero(2, 6);
			input.dCoordsDParams(0, 0) = 1.0;
			input.dCoordsDParams(1, 0) = 2.0;
			input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6);
			input.currentCoordsKm = Eigen::Vector2d(10.0, 20.0);
			input.targetCoordsKm = Eigen::Vector2d(2.0, 3.0);

			const VirtualImpactorRefinementResult result = refineVirtualImpactorWithObservations(input);

			EXPECT_TRUE(result.attempted);
			EXPECT_FALSE(result.valid);
			EXPECT_TRUE(result.badConditioning);
		}

		TEST(VirtualImpactorRefinementTest, BasisAcceptsHighlyStretchedFullRankMap)
		{
			ImpactRiskConfig cfg = refinementConfig();

			VirtualImpactorRefinementInput input;
			input.cfg = &cfg;
			input.initialParameters = Eigen::VectorXd::Zero(6);
			input.dCoordsDParams = Eigen::MatrixXd::Zero(2, 6);
			input.dCoordsDParams(0, 0) = 1.0e8;
			input.dCoordsDParams(1, 1) = 1.0;
			input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6);
			input.currentCoordsKm = Eigen::Vector2d(10.0, 20.0);
			input.targetCoordsKm = Eigen::Vector2d(2.0, 3.0);

			const VirtualImpactorRefinementResult result = refineVirtualImpactorWithObservations(input);

			ASSERT_TRUE(result.attempted);
			ASSERT_TRUE(result.valid);
			EXPECT_FALSE(result.badConditioning);
			const Eigen::Vector2d reached =
			    input.currentCoordsKm + input.dCoordsDParams * (result.parameters - input.initialParameters);
			EXPECT_TRUE(reached.isApprox(input.targetCoordsKm, 1.0e-12));
		}

		TEST(VirtualImpactorRefinementTest, UsesAuUnitsForAdaptedBasisTolerance)
		{
			ImpactRiskConfig cfg = refinementConfig();

			VirtualImpactorRefinementInput input;
			input.cfg = &cfg;
			input.initialParameters = Eigen::VectorXd::Zero(6);
			input.dCoordsDParams = Eigen::MatrixXd::Zero(2, 6);
			input.dCoordsDParams(0, 0) = 5.0e15;
			input.dCoordsDParams(1, 0) = 1.0e15;
			input.dCoordsDParams(1, 1) = 1.0e10;
			input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6);
			input.currentCoordsKm = Eigen::Vector2d(10.0, 20.0);
			input.targetCoordsKm = Eigen::Vector2d(2.0, 3.0);

			const VirtualImpactorRefinementResult result = refineVirtualImpactorWithObservations(input);

			ASSERT_TRUE(result.valid);
			EXPECT_FALSE(result.badConditioning);
			const Eigen::Vector2d reached =
			    input.currentCoordsKm + input.dCoordsDParams * (result.parameters - input.initialParameters);
			EXPECT_TRUE(reached.isApprox(input.targetCoordsKm, 1.0e-9));
		}

		TEST(VirtualImpactorRefinementTest, SemilinearStepUsesCovarianceRegression)
		{
			ImpactRiskConfig cfg = refinementConfig();

			VirtualImpactorRefinementInput input;
			input.cfg = &cfg;
			input.initialParameters = Eigen::VectorXd::Zero(6);
			input.dCoordsDParams = Eigen::MatrixXd::Zero(2, 6);
			input.dCoordsDParams(0, 0) = 1.0;
			input.dCoordsDParams(1, 1) = 1.0;
			input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6);
			input.parameterCovariance(0, 0) = 2.0;
			input.parameterCovariance(1, 1) = 3.0;
			input.parameterCovariance(0, 2) = 0.5;
			input.parameterCovariance(2, 0) = 0.5;
			input.parameterCovariance(1, 3) = 0.75;
			input.parameterCovariance(3, 1) = 0.75;
			input.currentCoordsKm = Eigen::Vector2d(10.0, 20.0);
			input.targetCoordsKm = Eigen::Vector2d(2.0, 3.0);

			const VirtualImpactorRefinementResult result = refineVirtualImpactorWithObservations(input);

			ASSERT_TRUE(result.valid);
			EXPECT_NEAR(result.parameters(0), -8.0, 1.0e-12);
			EXPECT_NEAR(result.parameters(1), -17.0, 1.0e-12);
			EXPECT_NEAR(result.parameters(2), 4.0, 1.0e-12);
			EXPECT_NEAR(result.parameters(3), 12.75, 1.0e-12);
			EXPECT_NEAR(result.parameters(4), 0.0, 1.0e-12);
			EXPECT_NEAR(result.parameters(5), 0.0, 1.0e-12);
		}

	}  // namespace
}  // namespace orbital_core::hazard
