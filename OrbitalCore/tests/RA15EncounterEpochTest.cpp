#include <gtest/gtest.h>
#include "integrators/RA15Integrator.h"
#include "orbital_core/nongrav/NonGravModel.h"

namespace {
	struct OffsetObserved {};

	struct LegacyObserved {};

	class UniformAcceleration final : public IIntegrationDynamics {
	   public:
		Vec3D accelerationAtMjd(const Vec3D&, const Vec3D&, double, const NonGravModel&) const override
		{
			return {0.0001, -0.00002, 0.000003};
		}

		void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D& v, double t, const NonGravModel& ng, Vec3D& a,
		                                Eigen::Matrix3d& dr, Eigen::Matrix3d& dv) const override
		{
			a = accelerationAtMjd(r, v, t, ng);
			dr.setZero();
			dv.setZero();
		}
	};

	class RA15EncounterEpochTest : public ::testing::TestWithParam<int> {};
}  // namespace

TEST_P(RA15EncounterEpochTest, ControllerPreservesEpochAndSmallOffset)
{
	constexpr double epochMjd = 59169.771532767;
	constexpr double firstStep = 1.708207647081121;
	const std::array<double, 1> targets{epochMjd + 30.0};
	UniformAcceleration dynamics;
	NonGravModel ng;
	RA15Integrator integrator;
	IntegrationOptions options;
	options.initialStepDays = firstStep;
	bool called = false;
	options.nextEncounterStep = [](double, double, const IntegratedState&) -> std::optional<double> {
		throw LegacyObserved{};
	};
	options.nextEncounterStepOffset = [&](double step, double base, double offset,
	                                      const IntegratedState&) -> std::optional<double> {
		called = true;
		EXPECT_DOUBLE_EQ(base, epochMjd);
		EXPECT_DOUBLE_EQ(step, firstStep);
		EXPECT_DOUBLE_EQ(offset, firstStep);
		EXPECT_DOUBLE_EQ(base + offset, epochMjd + firstStep);
		volatile double asJd = epochMjd + 2400000.5;
		EXPECT_NE(static_cast<double>(asJd) - 2400000.5, epochMjd);
		throw OffsetObserved{};
	};
	const auto run = [&] {
		switch (GetParam()) {
			case 0: {
				const IntegratedState initial{{1.0, 0.1, 0.01}, {0.0, 0.01, 0.0}};
				std::array<IntegratedState, 1> output;
				integrator.integrateState(dynamics, ng, initial, epochMjd, targets, output, options);
				break;
			}
			case 1: {
				IntegratedStateWithStm initial;
				initial.R = {1.0, 0.1, 0.01};
				initial.V = {0.0, 0.01, 0.0};
				std::array<IntegratedStateWithStm, 1> output;
				integrator.integrateStateWithStm(dynamics, ng, initial, epochMjd, targets, output, options);
				break;
			}
			default: {
				IntegratedStateFull initial;
				initial.R = {1.0, 0.1, 0.01};
				initial.V = {0.0, 0.01, 0.0};
				std::array<IntegratedStateFull, 1> output;
				integrator.integrateStateWithStmAndSens(dynamics, ng, initial, epochMjd, targets, output, options);
				break;
			}
		}
	};
	EXPECT_THROW(run(), OffsetObserved);
	EXPECT_TRUE(called);
}

INSTANTIATE_TEST_SUITE_P(StateStmSens, RA15EncounterEpochTest, ::testing::Values(0, 1, 2));
