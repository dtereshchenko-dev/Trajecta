#include <gtest/gtest.h>

#include "NBodyPropagator.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Types.h"
#include "integrators/IAS15Integrator.h"
#include "integrators/RA15Integrator.h"
#include "orbital_core/forces/SunGravity.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <vector>

namespace {

	[[nodiscard]] double specificEnergy(const Vec3D& r, const Vec3D& v) noexcept
	{
		const double mu = Constants::GAUSS_K * Constants::GAUSS_K;
		return 0.5 * v.dot(v) - mu / r.magnitude();
	}

	[[nodiscard]] Vec3D angularMomentum(const Vec3D& r, const Vec3D& v) noexcept { return r.cross(v); }

	class KeplerDynamics final : public IIntegrationDynamics {
	   public:
		Vec3D accelerationAtMjd(const Vec3D& r, const Vec3D&, double, const NonGravModel&) const override
		{
			return r * (-Constants::GAUSS_K * Constants::GAUSS_K / std::pow(r.magnitude(), 3));
		}

		void accelerationAndPartialsMjd(const Vec3D& r, const Vec3D&, double, const NonGravModel&, Vec3D& a,
		                                Eigen::Matrix3d& dr, Eigen::Matrix3d& dv) const override
		{
			orbital_core::forces::ForceContext c{};
			c.R = r;
			c.centralNewtonianGm = Constants::GAUSS_K * Constants::GAUSS_K;
			a = {};
			dr.setZero();
			dv.setZero();
			orbital_core::forces::SunGravity{}.evaluate(c, a, dr, dv);
		}
	};

	// With x(0)=0 and vx=1, every stage position is an independent elapsed clock.
	class ClockDynamics final : public IIntegrationDynamics {
	   public:
		mutable double maximumStageDrift = 0.0;

		Vec3D accelerationAtMjd(const Vec3D&, const Vec3D&, double, const NonGravModel&) const override
		{
			ADD_FAILURE() << "IAS15 collapsed the force epoch";
			return {};
		}

		Vec3D accelerationAtOffsetMjd(const Vec3D& r, const Vec3D&, double, double offset,
		                              const NonGravModel&) const override
		{
			maximumStageDrift = std::max(maximumStageDrift, std::abs(r.x - offset));
			return {};
		}

		void accelerationAndPartialsMjd(const Vec3D&, const Vec3D&, double, const NonGravModel&, Vec3D& a,
		                                Eigen::Matrix3d& dr, Eigen::Matrix3d& dv) const override
		{
			ADD_FAILURE() << "IAS15 collapsed the variational force epoch";
			a = {};
			dr.setZero();
			dv.setZero();
		}

		void accelerationAndPartialsOffsetMjd(const Vec3D& r, const Vec3D& v, double epoch, double offset,
		                                      const NonGravModel& ng, Vec3D& a, Eigen::Matrix3d& dr,
		                                      Eigen::Matrix3d& dv) const override
		{
			a = accelerationAtOffsetMjd(r, v, epoch, offset, ng);
			dr.setZero();
			dv.setZero();
		}

		Eigen::Matrix<double, 3, 4> nonGravPartialsOffsetMjd(const Vec3D&, const Vec3D&, double, double,
		                                                     const NonGravModel&) const override
		{
			return Eigen::Matrix<double, 3, 4>::Zero();
		}
	};

	enum class ClockStateKind { Nominal, Stm, Sensitivities };

	void checkCompensatedClock(ClockStateKind kind, double direction, bool restart)
	{
		constexpr double epoch = 60000.0;
		constexpr double span = 16385.0;
		constexpr double smallStep = 0.1;
		constexpr int smallStepCount = 4000;
		constexpr double tolerance = 8.0 * std::numeric_limits<double>::epsilon() * span;
		const std::vector<double> targets{epoch + direction * (8192.0 + 0.03), epoch + direction * (8192.0 + 100.03),
		                                  epoch + direction * (8192.0 + 399.99), epoch + direction * span};
		ClockDynamics dynamics;
		IAS15Integrator engine;
		NonGravModel ng;
		ng.setEnabled(false);
		IntegrationOptions options;
		options.initialStepDays = 8192.0;
		options.restartOnEncounterExit = restart;
		int accepted = 0;
		double maximumAcceptedDrift = 0.0;
		options.nextEncounterStepOffset = [&](double, double base, double offset, const IntegratedState& state) {
			EXPECT_DOUBLE_EQ(base, epoch);
			maximumAcceptedDrift = std::max(maximumAcceptedDrift, std::abs(state.R.x - offset));
			return ++accepted <= smallStepCount ? std::optional<double>(direction * smallStep) : std::nullopt;
		};
		const auto verify = [&](const auto& states) {
			for (std::size_t i = 0; i < states.size(); ++i) {
				EXPECT_NEAR(states[i].R.x, targets[i] - epoch, tolerance) << i;
				EXPECT_DOUBLE_EQ(states[i].V.x, 1.0);
			}
		};
		IntegrationReport report;
		if (kind == ClockStateKind::Nominal) {
			std::vector<IntegratedState> out(targets.size());
			report = engine.integrateState(dynamics, ng, {{}, {1, 0, 0}}, epoch, targets, out, options);
			verify(out);
		} else if (kind == ClockStateKind::Stm) {
			IntegratedStateWithStm initial;
			initial.V = {1, 0, 0};
			std::vector<IntegratedStateWithStm> out(targets.size());
			report = engine.integrateStateWithStm(dynamics, ng, initial, epoch, targets, out, options);
			verify(out);
			for (std::size_t i = 0; i < out.size(); ++i)
				EXPECT_NEAR(out[i].Phi(0, 3), targets[i] - epoch, tolerance);
		} else {
			IntegratedStateFull initial;
			initial.V = {1, 0, 0};
			std::vector<IntegratedStateFull> out(targets.size());
			report = engine.integrateStateWithStmAndSens(dynamics, ng, initial, epoch, targets, out, options);
			verify(out);
			for (std::size_t i = 0; i < out.size(); ++i) {
				EXPECT_NEAR(out[i].Phi(0, 3), targets[i] - epoch, tolerance);
				EXPECT_DOUBLE_EQ(out[i].Sens.norm(), 0.0);
			}
		}
		EXPECT_EQ(report.completedOutputs, targets.size());
		EXPECT_EQ(report.encounterRestarts, restart ? 1u : 0u);
		EXPECT_NEAR(report.finalState.R.x, direction * span, tolerance);
		EXPECT_LT(dynamics.maximumStageDrift, tolerance);
		EXPECT_LT(maximumAcceptedDrift, tolerance);
	}

	std::vector<IntegratedState> propagateKepler(IntegratorKind kind, const Vec3D& r, const Vec3D& v, double fromJd,
	                                             double toJd, double step = 10.0, IntegrationReport* report = nullptr)
	{
		std::unique_ptr<IIntegrator> engine;
		if (kind == IntegratorKind::IAS15)
			engine = std::make_unique<IAS15Integrator>();
		else
			engine = std::make_unique<RA15Integrator>();
		IntegrationOptions options;
		options.initialStepDays = step;
		const std::array<double, 1> dates{toJd - 2400000.5};
		std::vector<IntegratedState> out(1);
		const auto result =
		    engine->integrateState(KeplerDynamics{}, NonGravModel{}, {r, v}, fromJd - 2400000.5, dates, out, options);
		EXPECT_EQ(result.completedOutputs, 1u);
		if (report)
			*report = result;
		return out;
	}

	TEST(IAS15ClockTest, RetainsSmallStepsAfterLargeElapsedTime)
	{
		for (const auto kind : {ClockStateKind::Nominal, ClockStateKind::Stm, ClockStateKind::Sensitivities}) {
			for (const double direction : {1.0, -1.0}) {
				for (const bool restart : {false, true}) {
					SCOPED_TRACE(::testing::Message() << "state=" << static_cast<int>(kind)
					                                  << " direction=" << direction << " restart=" << restart);
					checkCompensatedClock(kind, direction, restart);
				}
			}
		}
	}

}  // namespace

class IAS15IntegratorTest : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (getExecutableDir().empty()) {
			auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
			std::filesystem::path candidates[] = {
			    srcDir.parent_path() / "data",
			    srcDir / "data",
			    srcDir / "release",
			    srcDir / "debug",
			};
			for (const auto& dir : candidates) {
				if (std::filesystem::exists(dir / "linux_p1550p2650.440")) {
					getExecutableDir() = dir;
					break;
				}
			}
		}
		de440 = std::make_unique<DE440Source>();
		if (de440->isAvailable())
			setDefaultEphemerisSource(de440.get());

		sb441 = std::make_unique<SpkReader>();
		if (sb441->isAvailable())
			setDefaultAsteroidSource(sb441.get());
	}

	void TearDown() override
	{
		setDefaultAsteroidSource(nullptr);
		setDefaultEphemerisSource(nullptr);
	}

	[[nodiscard]] static NBodyPropagator keplerOnly(IntegratorKind kind)
	{
		NBodyPropagator p(kind);
		p.setAllPlanetsEnabled(false);
		return p;
	}

	std::unique_ptr<DE440Source> de440;
	std::unique_ptr<SpkReader> sb441;
	NBodyPropagator propagator{IntegratorKind::IAS15};
};

TEST_F(IAS15IntegratorTest, IntegratorSelectionIsExplicit)
{
	EXPECT_EQ(propagator.integratorKind(), IntegratorKind::IAS15);
	EXPECT_STREQ(propagator.integratorName(), "IAS15");

	propagator.setIntegratorKind(IntegratorKind::RA15);
	EXPECT_EQ(propagator.integratorKind(), IntegratorKind::RA15);
	EXPECT_STREQ(propagator.integratorName(), "RA15");

	propagator.setIntegratorKind(IntegratorKind::IAS15);
	EXPECT_STREQ(propagator.integratorName(), "IAS15");
}

TEST_F(IAS15IntegratorTest, ZeroPropagationReturnsInitial)
{
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, 0.0172, 0.0);
	constexpr double epoch = 2460000.5;

	const Vec3D result = propagator.propagatePosition(R, V, epoch, epoch);

	EXPECT_NEAR(result.x, R.x, 1e-10);
	EXPECT_NEAR(result.y, R.y, 1e-10);
	EXPECT_NEAR(result.z, R.z, 1e-10);
}

TEST_F(IAS15IntegratorTest, CircularOrbitPeriod)
{
	setDefaultEphemerisSource(nullptr);
	const double vCirc = Constants::GAUSS_K;
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, vCirc, 0.0);
	constexpr double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(false);

	const double period = 2.0 * Constants::PI / Constants::GAUSS_K;
	const Vec3D result = propagator.propagatePosition(R, V, epoch, epoch + period);

	EXPECT_NEAR(result.x, R.x, 0.001) << "x after full period";
	EXPECT_NEAR(result.y, R.y, 0.001) << "y after full period";
	EXPECT_NEAR(result.z, R.z, 0.001) << "z after full period";
}

TEST_F(IAS15IntegratorTest, QuarterPeriodCircular)
{
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K, 0.0);
	constexpr double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(false);

	const double quarter = 0.5 * Constants::PI / Constants::GAUSS_K;
	const Vec3D result = propagator.propagatePosition(R, V, epoch, epoch + quarter);

	EXPECT_NEAR(result.x, 0.0, 0.01);
	EXPECT_NEAR(result.y, 1.0, 0.01);
}

TEST_F(IAS15IntegratorTest, BackwardPropagation)
{
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K, 0.0);
	constexpr double epoch = 2460000.5;

	propagator.setAllPlanetsEnabled(false);

	const Vec3D forward = propagator.propagatePosition(R, V, epoch, epoch + 40.0);
	NBodyPropagator back = keplerOnly(IntegratorKind::IAS15);
	const auto forwardFull = propagator.propagateAllWithVelocity(R, V, epoch, {epoch + 40.0});
	ASSERT_EQ(forwardFull.size(), 1u);
	const Vec3D returned = back.propagatePosition(forwardFull[0].R, forwardFull[0].V, epoch + 40.0, epoch);

	EXPECT_LT((returned - R).magnitude(), 1.0e-11) << "forward endpoint " << forward.magnitude();
}

TEST_F(IAS15IntegratorTest, BatchMatchesIndependentTargets)
{
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K, 0.0);
	constexpr double epoch = 2460000.5;
	const std::vector<double> targets{epoch + 25.0, epoch - 7.0, epoch, epoch + 2.0, epoch - 19.0, epoch + 25.0};

	propagator.setAllPlanetsEnabled(false);
	const auto batch = propagator.propagateAllWithVelocity(R, V, epoch, targets);

	ASSERT_EQ(batch.size(), targets.size());
	for (std::size_t index = 0; index < targets.size(); ++index) {
		NBodyPropagator independent = keplerOnly(IntegratorKind::IAS15);
		const auto single = independent.propagateAllWithVelocity(R, V, epoch, {targets[index]});
		ASSERT_EQ(single.size(), 1u);
		EXPECT_LT((batch[index].R - single[0].R).magnitude(), 2.0e-11) << "position at target " << index;
		EXPECT_LT((batch[index].V - single[0].V).magnitude(), 2.0e-12) << "velocity at target " << index;
	}
}

TEST_F(IAS15IntegratorTest, BatchStmMatchesIndependentTargets)
{
	const Vec3D R(1.2, 0.1, -0.02);
	const Vec3D V(-0.002, 0.015, 0.0003);
	constexpr double epoch = 2460000.5;
	const std::vector<double> targets{epoch + 3.0, epoch + 12.0};

	propagator.setAllPlanetsEnabled(false);
	const auto batch = propagator.propagateAllWithStm(R, V, epoch, targets);

	ASSERT_EQ(batch.size(), targets.size());
	for (std::size_t index = 0; index < targets.size(); ++index) {
		NBodyPropagator independent = keplerOnly(IntegratorKind::IAS15);
		const auto single = independent.propagateAllWithStm(R, V, epoch, {targets[index]});
		ASSERT_EQ(single.size(), 1u);
		EXPECT_LT((batch[index].R - single[0].R).magnitude(), 2.0e-11);
		EXPECT_LT((batch[index].V - single[0].V).magnitude(), 2.0e-12);
		EXPECT_LT((batch[index].Phi - single[0].Phi).norm(), 2.0e-9);
	}
}

TEST_F(IAS15IntegratorTest, ChainedStmMatchesContinuousStm)
{
	const Vec3D R(1.2, 0.1, -0.02);
	const Vec3D V(-0.002, 0.015, 0.0003);
	constexpr double epoch = 2460000.5;
	const std::vector<double> targets{epoch + 3.0, epoch + 12.0, epoch + 40.0};

	propagator.setAllPlanetsEnabled(false);

	IntegrationReport continuousReport;
	IntegrationReport chainedReport;
	const auto continuous =
	    propagator.propagateAllWithStm(R, V, epoch, targets, propagator.nonGrav(), &continuousReport);
	const auto chained =
	    propagator.propagateAllWithStmChained(R, V, epoch, targets, propagator.nonGrav(), &chainedReport);

	ASSERT_EQ(chained.size(), continuous.size());
	EXPECT_EQ(chainedReport.completedOutputs, targets.size());
	for (std::size_t index = 0; index < targets.size(); ++index) {
		EXPECT_LT((chained[index].R - continuous[index].R).magnitude(), 2.0e-11);
		EXPECT_LT((chained[index].V - continuous[index].V).magnitude(), 2.0e-12);
		EXPECT_LT((chained[index].Phi - continuous[index].Phi).norm(), 2.0e-9);
	}
}

TEST_F(IAS15IntegratorTest, AgreesWithRa15OnAPerturbedArc)
{
	if (!de440->isAvailable())
		GTEST_SKIP() << "DE440 not available";

	const Vec3D R(1.2, 0.1, -0.02);
	const Vec3D V(-0.002, 0.015, 0.0003);
	constexpr double epoch = 2460000.5;
	const std::vector<double> targets{epoch + 100.0, epoch + 365.0, epoch + 1000.0};

	NBodyPropagator ra15(IntegratorKind::RA15);
	NBodyPropagator ias15(IntegratorKind::IAS15);

	const auto a = ra15.propagateAllWithVelocity(R, V, epoch, targets);
	const auto b = ias15.propagateAllWithVelocity(R, V, epoch, targets);

	ASSERT_EQ(a.size(), b.size());
	for (std::size_t i = 0; i < a.size(); ++i) {
		const double positionAu = (a[i].R - b[i].R).magnitude();
		EXPECT_LT(positionAu, 1.0e-9) << "target " << i << " differs by "
		                              << positionAu * Constants::AU_IN_METERS / 1000.0 << " km";
		EXPECT_LT((a[i].V - b[i].V).magnitude(), 1.0e-11) << "target " << i;
	}
}

TEST_F(IAS15IntegratorTest, RejectsAndRedoesAnOversizedStep)
{
	setDefaultEphemerisSource(nullptr);
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K, 0.0);
	constexpr double epoch = 2460000.5;

	NBodyPropagator p = keplerOnly(IntegratorKind::IAS15);
	p.setStepSize(400.0);

	IntegrationReport report;
	const auto out = p.propagateAllWithVelocity(R, V, epoch, {epoch + 500.0}, p.nonGrav(), &report);

	ASSERT_EQ(out.size(), 1u);
	EXPECT_GT(report.rejectedSteps, 0u) << "an oversized first step was accepted";
	EXPECT_GT(report.acceptedSteps, 0u);

	NBodyPropagator reference = keplerOnly(IntegratorKind::IAS15);
	reference.setStepSize(5.0);
	const auto expected = reference.propagateAllWithVelocity(R, V, epoch, {epoch + 500.0});
	EXPECT_LT((out[0].R - expected[0].R).magnitude(), 1.0e-10);
}

TEST_F(IAS15IntegratorTest, ConservesCircularInvariantsOverManyRevolutions)
{
	setDefaultEphemerisSource(nullptr);
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K, 0.0);
	constexpr double epoch = 2460000.5;
	const double period = 2.0 * Constants::PI / Constants::GAUSS_K;

	const double span = 300.0 * period;

	const auto out = propagateKepler(IntegratorKind::IAS15, R, V, epoch, epoch + span);
	ASSERT_EQ(out.size(), 1u);

	const double e0 = specificEnergy(R, V);
	const double e1 = specificEnergy(out[0].R, out[0].V);
	const Vec3D h0 = angularMomentum(R, V);
	const Vec3D h1 = angularMomentum(out[0].R, out[0].V);

	EXPECT_LT(std::abs((e1 - e0) / e0), 1.0e-13) << "energy drift over 300 orbits";
	EXPECT_LT((h1 - h0).magnitude() / h0.magnitude(), 1.0e-13) << "angular momentum drift over 300 orbits";
}

TEST_F(IAS15IntegratorTest, RoundTripReturnsToTheStart)
{
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K * 1.303840, 0.0);
	constexpr double epoch = 2460000.5;
	const double period = 2.0 * Constants::PI / Constants::GAUSS_K;
	const double span = 100.0 * period;

	const auto forward = propagateKepler(IntegratorKind::IAS15, R, V, epoch, epoch + span);
	ASSERT_EQ(forward.size(), 1u);

	const auto returned = propagateKepler(IntegratorKind::IAS15, forward[0].R, forward[0].V, epoch + span, epoch);
	ASSERT_EQ(returned.size(), 1u);

	const double closure = (returned[0].R - R).magnitude();
	EXPECT_LT(closure, 1.0e-11) << "round trip over 100 eccentric orbits closed to "
	                            << closure * Constants::AU_IN_METERS / 1000.0 << " km";
}

TEST_F(IAS15IntegratorTest, FixedAccuracyWithDifferentInitialSteps)
{
	const Vec3D R(1.0, 0.0, 0.0);
	const Vec3D V(0.0, Constants::GAUSS_K * 1.140175, 0.0);
	constexpr double epoch = 2460000.5;
	const double period = 2.0 * Constants::PI / Constants::GAUSS_K;
	const double span = 20.0 * period;

	const auto closureAt = [&](double initialStepDays, IntegrationReport& report) {
		EXPECT_DOUBLE_EQ(IAS15Integrator::kEpsilon, 1.0e-9);
		const auto f = propagateKepler(IntegratorKind::IAS15, R, V, epoch, epoch + span, initialStepDays, &report);
		const auto b = propagateKepler(IntegratorKind::IAS15, f[0].R, f[0].V, epoch + span, epoch, initialStepDays);
		return (b[0].R - R).magnitude();
	};

	for (const double initialStepDays : {0.1, 10.0}) {
		IntegrationReport report;
		EXPECT_LT(closureAt(initialStepDays, report), 1.0e-11);
		EXPECT_GT(report.acceptedSteps, 0u);
	}
}

TEST_F(IAS15IntegratorTest, MatchesRa15OnDifferentKeplerOrbits)
{
	const double period = 2.0 * Constants::PI / Constants::GAUSS_K;
	constexpr double epoch = 2460000.5;
	const double vc = Constants::GAUSS_K;

	struct Case {
		const char* name;
		Vec3D V;
		double spanPeriods;
	};

	const std::vector<Case> cases{
	    {"circular  e=0.0   10 rev", {0, vc, 0}, 10.0},
	    {"eccentric e=0.3   10 rev", {0, vc * 1.140175, 0}, 10.0},
	    {"eccentric e=0.7   10 rev", {0, vc * 1.303840, 0}, 10.0},
	    {"circular  e=0.0  300 rev", {0, vc, 0}, 300.0},
	    {"eccentric e=0.7  300 rev", {0, vc * 1.303840, 0}, 300.0},
	};

	const Vec3D R(1.0, 0.0, 0.0);

	for (const auto& c : cases) {
		const double span = c.spanPeriods * period;
		IntegrationReport ra15Report, ias15Report;
		const auto ra15 = propagateKepler(IntegratorKind::RA15, R, c.V, epoch, epoch + span, 10.0, &ra15Report);
		const auto ias15 = propagateKepler(IntegratorKind::IAS15, R, c.V, epoch, epoch + span, 10.0, &ias15Report);
		ASSERT_EQ(ra15.size(), 1u) << c.name;
		ASSERT_EQ(ias15.size(), 1u) << c.name;
		const double separationKm = (ra15[0].R - ias15[0].R).magnitude() * Constants::AU_IN_METERS / 1000.0;
		EXPECT_GT(ra15Report.acceptedSteps, 0u);
		EXPECT_GT(ias15Report.acceptedSteps, 0u);
		EXPECT_LT(separationKm, 1.0) << c.name;
	}
}
