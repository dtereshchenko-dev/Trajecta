#include <gtest/gtest.h>

#include "orbital_core/optimizers/EigenLevenbergMarquardt.h"

#include <cmath>
#include <limits>
#include <utility>

namespace {
	namespace opt = orbital_core::optimizers;

	class ScalarProblem : public opt::ILSQProblem {
	   public:
		int numParams() const noexcept override { return 1; }

		int numResiduals() const noexcept override { return 1; }

		void residuals(const std::vector<double>& x, std::vector<double>& r) const override
		{
			++evaluations;
			r = {x[0] - target};
		}

		void jacobian(const std::vector<double>&, Eigen::MatrixXd& j) const override
		{
			++jacobians;
			j = Eigen::MatrixXd::Ones(1, 1);
		}

		double displayCost(const std::vector<double>& x) const override
		{
			++displayCalls;
			return std::abs(x[0] - target);
		}

		double target = 2.0;
		mutable int evaluations = 0;
		mutable int jacobians = 0;
		mutable int displayCalls = 0;
	};

	class RosenbrockProblem : public opt::ILSQProblem {
	   public:
		int numParams() const noexcept override { return 2; }

		int numResiduals() const noexcept override { return 2; }

		void residuals(const std::vector<double>& x, std::vector<double>& r) const override
		{
			r = {10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]};
		}

		void jacobian(const std::vector<double>& x, Eigen::MatrixXd& j) const override
		{
			j.resize(2, 2);
			j << -20.0 * x[0], 10.0, -1.0, 0.0;
		}

		double displayCost(const std::vector<double>& x) const override
		{
			std::vector<double> r;
			residuals(x, r);
			return std::hypot(r[0], r[1]);
		}
	};

}  // namespace

TEST(EigenLevenbergMarquardtTest, SolvesLinearResiduals)
{
	ScalarProblem problem;
	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 2.0, 1e-12);
	EXPECT_LT(result.iterations, 10);
	EXPECT_EQ(problem.displayCalls, 1);
	EXPECT_DOUBLE_EQ(result.finalCost, 0.0);
}

TEST(EigenLevenbergMarquardtTest, RestartsAfterAcceptedObjectiveChange)
{
	ScalarProblem problem;
	int changes = 0;
	const auto result =
	    opt::EigenLevenbergMarquardt().optimize({0.0}, problem, {}, [&](const std::vector<double>& accepted) {
		    if (changes == 0) {
			    EXPECT_NEAR(accepted[0], 2.0, 1e-12);
			    problem.target = 7.0;
			    ++changes;
			    return true;
		    }
		    return false;
	    });
	EXPECT_EQ(changes, 1);
	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 7.0, 1e-12);

	ScalarProblem invalid;
	const auto failed = opt::EigenLevenbergMarquardt().optimize({0.0}, invalid, {}, [&](const std::vector<double>&) {
		invalid.target = std::numeric_limits<double>::infinity();
		return true;
	});
	EXPECT_FALSE(failed.converged);
}

TEST(EigenLevenbergMarquardtTest, IncorrectJacobianDimensionsDoNotReachEigen)
{
	class InvalidProblem final : public ScalarProblem {
	   public:
		void jacobian(const std::vector<double>&, Eigen::MatrixXd& j) const override
		{
			j = Eigen::MatrixXd::Ones(rows, columns);
		}

		int rows = 1;
		int columns = 1;
	} problem;

	for (const auto dimensions : {std::pair{1, 2}, std::pair{2, 1}}) {
		problem.rows = dimensions.first;
		problem.columns = dimensions.second;
		const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
		EXPECT_FALSE(result.converged);
		EXPECT_FALSE(result.cancelled);
		EXPECT_DOUBLE_EQ(result.params[0], 0.0);
	}
}

TEST(EigenLevenbergMarquardtTest, ChangingTrialResidualCountKeepsAcceptedState)
{
	class InvalidProblem final : public ScalarProblem {
	   public:
		void residuals(const std::vector<double>& x, std::vector<double>& r) const override
		{
			ScalarProblem::residuals(x, r);
			if (x[0] != 0.0)
				r.push_back(x[0]);
		}
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_FALSE(result.converged);
	EXPECT_FALSE(result.cancelled);
	EXPECT_DOUBLE_EQ(result.params[0], 0.0);
}

TEST(EigenLevenbergMarquardtTest, InitialGuessDoesNotConstrainTheSolution)
{
	for (double initial : {-100.0, 0.0, 1.999, 2.0, 100.0}) {
		ScalarProblem problem;
		const auto result = opt::EigenLevenbergMarquardt().optimize({initial}, problem);
		EXPECT_TRUE(result.converged);
		EXPECT_NEAR(result.params[0], problem.target, 1e-12);
		EXPECT_NEAR(result.finalCost, 0.0, 1e-12);
	}
}

TEST(EigenLevenbergMarquardtTest, SolvesNonlinearResiduals)
{
	RosenbrockProblem problem;
	const auto result = opt::EigenLevenbergMarquardt().optimize({-1.2, 1.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 1.0, 1e-9);
	EXPECT_NEAR(result.params[1], 1.0, 1e-9);
}

TEST(EigenLevenbergMarquardtTest, HandlesDifferentParameterUnits)
{
	class ScaledProblem final : public opt::ILSQProblem {
	   public:
		int numParams() const noexcept override { return 2; }

		int numResiduals() const noexcept override { return 2; }

		void residuals(const std::vector<double>& x, std::vector<double>& r) const override
		{
			r = {x[0] - 2.0, x[1] * 1e14 + 3.0};
		}

		void jacobian(const std::vector<double>&, Eigen::MatrixXd& j) const override
		{
			j = Eigen::MatrixXd::Identity(2, 2);
			j(1, 1) = 1e14;
		}

		double displayCost(const std::vector<double>& x) const override
		{
			return std::hypot(x[0] - 2.0, x[1] * 1e14 + 3.0);
		}
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0, 0.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 2.0, 1e-12);
	EXPECT_NEAR(result.params[1], -3e-14, 1e-26);
}

TEST(EigenLevenbergMarquardtTest, InvalidTrialShrinksStepAndRetries)
{
	class GuardedProblem final : public ScalarProblem {
	   public:
		bool acceptParams(const std::vector<double>& x) const override
		{
			if (evaluations == 1 && x[0] > 1.5) {
				++rejections;
				return false;
			}
			return true;
		}

		mutable int rejections = 0;
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_GT(problem.rejections, 0);
	EXPECT_NEAR(result.params[0], 2.0, 1e-12);
}

TEST(EigenLevenbergMarquardtTest, AllInvalidTrialsKeepLastAcceptedState)
{
	class GuardedProblem final : public ScalarProblem {
	   public:
		bool acceptParams(const std::vector<double>& x) const override { return x[0] == 1.0; }
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({1.0}, problem);
	EXPECT_FALSE(result.converged);
	EXPECT_FALSE(result.cancelled);
	EXPECT_DOUBLE_EQ(result.params[0], 1.0);
}

TEST(EigenLevenbergMarquardtTest, CancelsBeforeEvaluation)
{
	ScalarProblem problem;
	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem, [] { return true; });
	EXPECT_TRUE(result.cancelled);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(problem.evaluations, 0);
	EXPECT_EQ(problem.displayCalls, 0);
}

TEST(EigenLevenbergMarquardtTest, CancelsAfterInitialResidualsWithoutPostprocessing)
{
	ScalarProblem problem;
	const auto result =
	    opt::EigenLevenbergMarquardt().optimize({0.0}, problem, [&] { return problem.evaluations > 0; });
	EXPECT_TRUE(result.cancelled);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(problem.jacobians, 0);
	EXPECT_EQ(problem.displayCalls, 0);
}

TEST(EigenLevenbergMarquardtTest, CancelsDuringJacobianWithoutTrial)
{
	ScalarProblem problem;
	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem, [&] { return problem.jacobians > 0; });
	EXPECT_TRUE(result.cancelled);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(problem.evaluations, 1);
	EXPECT_DOUBLE_EQ(result.params[0], 0.0);
}

TEST(EigenLevenbergMarquardtTest, CancelsDuringTrialWithoutAcceptingIt)
{
	ScalarProblem problem;
	const auto result =
	    opt::EigenLevenbergMarquardt().optimize({0.0}, problem, [&] { return problem.evaluations > 1; });
	EXPECT_TRUE(result.cancelled);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(problem.displayCalls, 0);
	EXPECT_DOUBLE_EQ(result.params[0], 0.0);
}

TEST(EigenLevenbergMarquardtTest, InvalidJacobianDoesNotClaimConvergence)
{
	class InvalidProblem final : public ScalarProblem {
	   public:
		void jacobian(const std::vector<double>&, Eigen::MatrixXd& j) const override
		{
			j = Eigen::MatrixXd::Constant(1, 1, std::numeric_limits<double>::quiet_NaN());
		}
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_FALSE(result.converged);
	EXPECT_FALSE(result.cancelled);
	EXPECT_DOUBLE_EQ(result.params[0], 0.0);
}

TEST(EigenLevenbergMarquardtTest, IncorrectResidualCountDoesNotReachEigen)
{
	class InvalidProblem final : public ScalarProblem {
	   public:
		void residuals(const std::vector<double>&, std::vector<double>& r) const override { r.clear(); }
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(problem.jacobians, 0);
}

TEST(EigenLevenbergMarquardtTest, FixedIterationLimitDoesNotClaimConvergence)
{
	class ExponentialProblem final : public ScalarProblem {
	   public:
		void residuals(const std::vector<double>& x, std::vector<double>& r) const override
		{
			++evaluations;
			r = {std::exp(x[0])};
		}

		void jacobian(const std::vector<double>& x, Eigen::MatrixXd& j) const override
		{
			j = Eigen::MatrixXd::Constant(1, 1, std::exp(x[0]));
		}

		double displayCost(const std::vector<double>& x) const override { return std::exp(x[0]); }
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(result.iterations, 100);
	EXPECT_LE(problem.evaluations, 1000);
}

TEST(EigenLevenbergMarquardtTest, NonfiniteInitialParametersAreRejected)
{
	ScalarProblem problem;
	EXPECT_THROW((void)opt::EigenLevenbergMarquardt().optimize({std::numeric_limits<double>::infinity()}, problem),
	             std::invalid_argument);
}

TEST(EigenLevenbergMarquardtTest, StationaryInitialPointConvergesWithoutMoving)
{
	ScalarProblem problem;
	const auto result = opt::EigenLevenbergMarquardt().optimize({2.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_DOUBLE_EQ(result.params[0], 2.0);
	EXPECT_EQ(problem.jacobians, 1);
}

TEST(EigenLevenbergMarquardtTest, CancellationDuringDisplayCostIsNotLost)
{
	ScalarProblem problem;
	const auto result =
	    opt::EigenLevenbergMarquardt().optimize({2.0}, problem, [&] { return problem.displayCalls > 0; });
	EXPECT_TRUE(result.cancelled);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(problem.displayCalls, 1);
}

TEST(EigenLevenbergMarquardtTest, DisplayCostRunsOnlyAfterOptimization)
{
	class FinalCostProblem final : public ScalarProblem {
	   public:
		double displayCost(const std::vector<double>& x) const override
		{
			EXPECT_DOUBLE_EQ(x[0], target);
			return ScalarProblem::displayCost(x);
		}
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_FALSE(result.cancelled);
	EXPECT_GT(result.iterations, 0);
	EXPECT_EQ(problem.displayCalls, 1);
}

TEST(EigenLevenbergMarquardtTest, NonfiniteTrialCanRetry)
{
	class NonfiniteTrial final : public ScalarProblem {
	   public:
		void residuals(const std::vector<double>& x, std::vector<double>& r) const override
		{
			ScalarProblem::residuals(x, r);
			if (evaluations == 2)
				r[0] = std::numeric_limits<double>::infinity();
		}
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_NEAR(result.params[0], 2.0, 1e-12);
}

TEST(EigenLevenbergMarquardtTest, SmallInitialCorrectionStillRunsToConvergence)
{
	ScalarProblem problem;
	const auto final = opt::EigenLevenbergMarquardt().optimize({1.95}, problem);
	EXPECT_TRUE(final.converged);
	EXPECT_GT(problem.evaluations, 1);
	EXPECT_DOUBLE_EQ(final.params[0], 2.0);
}

TEST(EigenLevenbergMarquardtTest, UsesEigenRelativeToleranceWithoutAdditionalCorrectionGate)
{
	class FlatObjective final : public ScalarProblem {
	   public:
		int numResiduals() const noexcept override { return 2; }

		void residuals(const std::vector<double>& x, std::vector<double>& r) const override { r = {1e12, x[0] - 2.0}; }

		void jacobian(const std::vector<double>&, Eigen::MatrixXd& j) const override
		{
			j.resize(2, 1);
			j << 0.0, 1.0;
		}
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({0.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_TRUE(std::isfinite(result.finalCost));
}

TEST(EigenLevenbergMarquardtTest, StationaryRankDeficientProblemUsesEigenStoppingStatus)
{
	class Degenerate final : public RosenbrockProblem {
	   public:
		void residuals(const std::vector<double>& x, std::vector<double>& r) const override
		{
			r = {x[0] + x[1] - 2.0, 0.0};
		}

		void jacobian(const std::vector<double>&, Eigen::MatrixXd& j) const override
		{
			j.resize(2, 2);
			j << 1.0, 1.0, 0.0, 0.0;
		}

		double displayCost(const std::vector<double>& x) const override { return std::abs(x[0] + x[1] - 2.0); }
	} problem;

	const auto result = opt::EigenLevenbergMarquardt().optimize({1.0, 1.0}, problem);
	EXPECT_TRUE(result.converged);
	EXPECT_DOUBLE_EQ(result.finalCost, 0.0);
}
