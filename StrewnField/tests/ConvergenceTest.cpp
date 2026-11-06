#include "Convergence.h"
#include "strewn_core/Simulator.h"
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <numbers>

namespace {
	using namespace strewn;

	ConvergenceConfig config()
	{
		ConvergenceConfig c;
		c.centerDeg = {49.8, 0.8};
		c.halfSizeKm = {1, 1};
		c.batchSize = 2;
		return c;
	}

	RunResult impact(double lat = 49.8, double lon = 0.8, double mass = 1, double count = 1)
	{
		RunResult r;
		Fragment f;
		f.state.massKg = mass;
		f.multiplicity = count;
		r.impacts.push_back({f, lat, lon, 0});
		return r;
	}

	TEST(Convergence, BaseTwoDivergenceNotDistance)
	{
		const std::array<double, 2> p{1, 0}, q{0, 1}, r{.5, .5};
		EXPECT_DOUBLE_EQ(detail::jsDivergence(p, p), 0);
		EXPECT_DOUBLE_EQ(detail::jsDivergence(p, q), 1);
		EXPECT_NEAR(detail::jsDivergence(p, r), .31127812445913283, 1e-15);
		EXPECT_DOUBLE_EQ(detail::jsDivergence(p, r), detail::jsDivergence(r, p));
	}

	TEST(Convergence, BothMapsAndThreeConsecutiveChecksRequired)
	{
		auto c = config();
		c.batchSize = 1;
		c.threshold = 1e-5;
		detail::ConvergenceMonitor monitor(c);
		RunSummary s;
		auto a = impact(49.795, .8, 1, 10);
		const auto b = impact(49.805, .8, 1, 10);
		a.impacts.push_back(b.impacts.front());
		s.completed = 1;
		monitor.update(a, s);
		EXPECT_FALSE(s.converged);
		EXPECT_FALSE(s.convergenceHistory.back().jsCount);
		a.impacts.back().fragment.state.massKg = 100;
		s.completed = 2;
		monitor.update(a, s);
		EXPECT_LT(*s.convergenceHistory.back().jsCount, 1e-14);
		EXPECT_GT(*s.convergenceHistory.back().jsMass, .1);
		EXPECT_EQ(s.convergenceHistory.back().stableChecks, 0U);
		a.impacts.back().fragment.state.massKg = 50.5;
		for (int n = 3; n <= 5; ++n) {
			s.completed = n;
			monitor.update(a, s);
			EXPECT_EQ(s.converged, n == 5);
		}
	}

	TEST(Convergence, EmptyPartialAndOutsideCannotConverge)
	{
		for (int mode = 0; mode < 6; ++mode) {
			detail::ConvergenceMonitor monitor(config());
			RunSummary summary;
			auto r = mode == 0 ? RunResult{} : impact();
			if (mode == 1)
				r.status = RunStatus::ParticleLimit;
			if (mode == 2)
				r.impacts.front().longitudeDeg = 10;
			if (mode == 3)
				r.status = RunStatus::Timeout;
			if (mode == 4)
				r.status = RunStatus::Failed;
			if (mode == 5)
				r.status = RunStatus::WallTimeout;
			for (int i = 1; i <= 10; ++i) {
				summary.completed = i;
				monitor.update(r, summary);
			}
			EXPECT_FALSE(summary.converged);
			EXPECT_EQ(summary.convergenceBlocked, mode != 0);
			EXPECT_EQ(summary.convergenceOutOfGridImpacts, mode == 2 ? 10U : 0U);
			EXPECT_FALSE(summary.convergenceHistory.back().jsCount);
		}
	}

	TEST(Convergence, MultiplicityMatchesExpandedIdenticalImpacts)
	{
		auto c = config();
		c.batchSize = 1;
		detail::ConvergenceMonitor grouped(c), expanded(c);
		RunSummary a, b;
		for (std::size_t i = 1; i <= 4; ++i) {
			auto g = impact(49.8, .8 + .001 * i, static_cast<double>(i), 20);
			auto e = g;
			e.impacts.front().fragment.multiplicity = 1;
			e.impacts.resize(20, e.impacts.front());
			a.completed = b.completed = i;
			grouped.update(g, a);
			expanded.update(e, b);
			if (i > 1) {
				EXPECT_DOUBLE_EQ(*a.convergenceHistory.back().jsCount, *b.convergenceHistory.back().jsCount);
				EXPECT_DOUBLE_EQ(*a.convergenceHistory.back().jsMass, *b.convergenceHistory.back().jsMass);
			}
		}
	}

	TEST(Convergence, InvalidConfigurationIsBounded)
	{
		for (int mode = 0; mode < 6; ++mode) {
			auto c = config();
			if (mode == 0)
				c.threshold = 0;
			if (mode == 1)
				c.cellM = 0;
			if (mode == 2)
				c.halfSizeKm[0] = 1e10;
			if (mode == 3)
				c.centerDeg[0] = std::numeric_limits<double>::quiet_NaN();
			if (mode == 4)
				c.batchSize = 0;
			if (mode == 5)
				c.kernelSigmaM = 1e10;
			EXPECT_THROW(detail::ConvergenceMonitor{c}, std::invalid_argument);
		}
	}

	TEST(Convergence, ThreadCountKeepsTheSameStoppedPrefix)
	{
		const auto env = Environment::idealized();
		PhysicsConfig physics;
		physics.phase = FlightPhase::DarkFlight;
		const GeodeticPosition p{49.8 * std::numbers::pi / 180, .8 * std::numbers::pi / 180, 0};
		const std::vector<EntryState> entries{
		    {ecefPosition(p), enuToEcef({0, 0, -30}, p.latitudeRad, p.longitudeRad), 59988}};
		MonteCarloConfig mc;
		mc.runs = 50;
		mc.convergence = config();
		for (unsigned workers : {1U, 15U}) {
			mc.workers = workers;
			std::vector<std::uint64_t> ids;
			const auto summary =
			    Simulator(env, physics).run(entries, mc, {}, [&](const auto& r) { ids.push_back(r.runId); });
			EXPECT_TRUE(summary.converged);
			EXPECT_FALSE(summary.cancelled);
			EXPECT_EQ(summary.completed, 8U);
			EXPECT_EQ(ids, (std::vector<std::uint64_t>{0, 1, 2, 3, 4, 5, 6, 7}));
			EXPECT_EQ(summary.convergenceHistory.size(), 4U);
		}
		mc.runs = 6;
		const auto limited = Simulator(env, physics).run(entries, mc);
		EXPECT_FALSE(limited.converged);
		EXPECT_FALSE(limited.cancelled);
		EXPECT_EQ(limited.completed, 6U);
		EXPECT_EQ(limited.convergenceHistory.back().stableChecks, 2U);
	}

}  // namespace
