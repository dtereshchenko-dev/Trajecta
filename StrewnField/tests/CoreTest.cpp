#include "strewn_core/FragmentationModel.h"
#include "strewn_core/Simulator.h"
#include "RunDeadline.h"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <thread>

namespace {
	using namespace strewn;
	constexpr double deg = std::numbers::pi / 180.0;

	EntryState entry(double height = 1000.0)
	{
		const GeodeticPosition position{49.8 * deg, 0.8 * deg, height};
		return {ecefPosition(position), enuToEcef({10, 0, -30}, position.latitudeRad, position.longitudeRad), 59988.125,
		        1.0, 42};
	}

	PhysicsConfig dark()
	{
		PhysicsConfig physics;
		physics.phase = FlightPhase::DarkFlight;
		return physics;
	}

	double distance(const Vector3& a, const Vector3& b)
	{
		double sum = 0;
		for (std::size_t i = 0; i < 3; ++i)
			sum += (a[i] - b[i]) * (a[i] - b[i]);
		return std::sqrt(sum);
	}

	void balanced(const MassBalance& b)
	{
		EXPECT_NEAR(b.initialKg, b.landedKg + b.ablatedKg + b.dustKg + b.residualKg + b.belowCutoffKg + b.unresolvedKg,
		            std::max(1e-10, b.initialKg * 2e-12));
	}

	Fragment fragment(double height, double mass = 0.1)
	{
		const auto initial = entry(height);
		Fragment f;
		f.state = {initial.positionEcefM, initial.velocityEcefMps, mass};
		f.strengthPa = 4e6;
		return f;
	}

	TEST(Geodesy, RoundTripsAndPoles)
	{
		for (double lat : {-90.0, -49.8, 0.0, 49.8, 90.0}) {
			for (double height : {-20.0, 100.0, 100000.0}) {
				const GeodeticPosition p{lat * deg, 0.8 * deg, height};
				const auto got = geodetic(ecefPosition(p));
				EXPECT_NEAR(got.latitudeRad, p.latitudeRad, 1e-10);
				EXPECT_NEAR(got.heightWgs84M, height, 1e-3);
			}
		}
		EXPECT_THROW(geodetic({0, 0, 0}), std::invalid_argument);
	}

	TEST(Environment, StandardAtmosphereAndNoVerticalWind)
	{
		EXPECT_NEAR(standardAtmosphere(0).densityKgM3, 1.224999, 1e-5);
		EXPECT_DOUBLE_EQ(standardAtmosphere(100000).temperatureK, 195.08);
		auto env = Environment::idealized();
		const auto a = env.air({49 * deg, 0.8 * deg, 10000});
		EXPECT_EQ(a.windEnuMps, (Vector3{0, 0, 0}));
		EXPECT_THROW(standardAtmosphere(150000), std::out_of_range);
	}

	TEST(Environment, TrilinearFieldsAndBilinearTerrain)
	{
		AtmosphereGrid a;
		a.latitudesDeg = {0, 1};
		a.longitudesDeg = {10, 11};
		a.heightCount = 2;
		a.heightStepM = 1000;
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < 2; ++k)
					a.samples.push_back({1.0 + i + 2 * j + 4 * k, 200.0 + i + j + k, {10.0 * i, 20.0 * j, 3.0 * k}});
		TerrainGrid t{{0, 1}, {10, 11}, {0, 20, 10, 30}};
		Environment env(std::move(a), std::move(t));
		const GeodeticPosition p{0.5 * deg, 10.5 * deg, 500};
		const auto value = env.air(p);
		EXPECT_DOUBLE_EQ(value.densityKgM3, 4.5);
		EXPECT_DOUBLE_EQ(value.temperatureK, 201.5);
		EXPECT_EQ(value.windEnuMps, (Vector3{5, 10, 1.5}));
		EXPECT_DOUBLE_EQ(env.ground(p), 15);
		EXPECT_THROW(env.air({3 * deg, 10.5 * deg, 500}), std::out_of_range);
	}

	TEST(Environment, RejectsMalformedGrid)
	{
		AtmosphereGrid a;
		a.heightCount = 2;
		EXPECT_THROW(Environment(a, {}), std::invalid_argument);
		const auto env = Environment::idealized();
		EXPECT_THROW(env.air({0, 0, std::numeric_limits<double>::quiet_NaN()}), std::invalid_argument);
		a.samples = {{1, 280, {}}, {0, 270, {}}};
		EXPECT_THROW(Environment(a, {}), std::invalid_argument);
	}

	TEST(Integrator, DarkFlightPreservesMassAndRefinesGround)
	{
		const auto env = Environment::idealized(123.0);
		const auto r = TrajectoryIntegrator(env, dark()).advance(fragment(1000), RunParameters{});
		EXPECT_EQ(r.event, FlightEvent::Ground);
		EXPECT_NEAR(geodetic(r.state.positionEcefM).heightWgs84M, 123, 1e-5);
		EXPECT_DOUBLE_EQ(r.state.massKg, 0.1);
		EXPECT_GT(r.timeS, 0.0);
		EXPECT_LT(r.stats.acceptedSteps, 1000U);
	}

	TEST(Integrator, CancellationAndTimeoutAreDistinct)
	{
		const auto env = Environment::idealized();
		std::stop_source stop;
		stop.request_stop();
		EXPECT_EQ(TrajectoryIntegrator(env, dark()).advance(fragment(1000), {}, stop.get_token()).event,
		          FlightEvent::Cancelled);
		IntegrationOptions options;
		options.maxFlightTimeS = 0.01;
		const auto result = TrajectoryIntegrator(env, dark(), options).advance(fragment(1000), {});
		EXPECT_EQ(result.event, FlightEvent::Timeout);
		EXPECT_NEAR(result.timeS, 0.01, 1e-14);
	}

	TEST(Integrator, InitialCutoffAndBreakup)
	{
		const auto env = Environment::idealized();
		EXPECT_EQ(TrajectoryIntegrator(env, {}).advance(fragment(1000, 0.0005), {}).event, FlightEvent::MassCutoff);
		auto f = fragment(10000, 100);
		f.state.velocityEcefMps = {0, 0, -15000};
		f.strengthPa = 10000;
		EXPECT_EQ(TrajectoryIntegrator(env, {}).advance(f, {}).event, FlightEvent::Breakup);
	}

	TEST(Integrator, StepRefinementIsStable)
	{
		const auto env = Environment::idealized();
		IntegrationOptions coarse, fine;
		fine.relativeTolerance = 1e-10;
		fine.positionToleranceM = 1e-4;
		fine.velocityToleranceMps = 1e-6;
		fine.maxStepS = 0.2;
		const auto a = TrajectoryIntegrator(env, dark(), coarse).advance(fragment(25000), {});
		const auto b = TrajectoryIntegrator(env, dark(), fine).advance(fragment(25000), {});
		EXPECT_EQ(a.event, FlightEvent::Ground);
		EXPECT_EQ(b.event, FlightEvent::Ground);
		EXPECT_LT(distance(a.state.positionEcefM, b.state.positionEcefM), 0.05);
		EXPECT_NEAR(a.timeS, b.timeS, 0.005);
	}

	TEST(Integrator, HighTerrainIsNotSkipped)
	{
		const auto env = Environment::idealized(5500);
		const auto result = TrajectoryIntegrator(env, dark()).advance(fragment(6000), {});
		EXPECT_EQ(result.event, FlightEvent::Ground);
		EXPECT_NEAR(geodetic(result.state.positionEcefM).heightWgs84M, 5500, 1e-5);
	}

	TEST(Integrator, NonFiniteDerivativeCannotBecomeAnAcceptedStep)
	{
		const auto env = Environment::idealized();
		RunParameters p;
		p.dragScale = 1e308;
		EXPECT_THROW(TrajectoryIntegrator(env, dark()).advance(fragment(1000), p), std::runtime_error);
	}

	TEST(Fragmentation, GroupingPreservesMassAndCount)
	{
		PhysicsConfig grouped, individual;
		individual.individualFragments = {100000, 100000};
		RunParameters p;
		p.massKg = 100;
		const auto parent = fragment(30000, 100);
		const auto a = FragmentationModel(grouped).groups(parent, p);
		const auto b = FragmentationModel(individual).groups(parent, p);
		const auto mass = [](const auto& groups) {
			double sum = 0;
			for (const auto& g : groups)
				sum += g.massKg * g.multiplicity;
			return sum;
		};
		double count = 0;
		for (const auto& g : a)
			count += g.multiplicity;
		EXPECT_NEAR(mass(a), mass(b), 1e-10);
		EXPECT_EQ(count, static_cast<double>(b.size()));
		EXPECT_LT(a.size(), b.size());
	}

	TEST(Fragmentation, MoscatiExamplesCloseAvailableMass)
	{
		struct Example {
			double dust, largest, beta;
			std::size_t count;
			std::array<double, 5> first;
			double remainder;
		};

		const Example examples[]{{.6,
		                          .15,
		                          .55,
		                          117,
		                          {60, 44.18036031004198, 34.001022218635185, 27.04541297357298, 22.07090694135014},
		                          .031322109851249},
		                         {.2,
		                          .45,
		                          .35,
		                          6,
		                          {360, 193.67621463332728, 116.40641388460833, 75.57526236401169, 51.93439001152369},
		                          2.407719106529015}};
		for (const auto& e : examples) {
			SCOPED_TRACE(e.beta);
			auto parent = fragment(30000, 1000);
			parent.multiplicity = 7;
			RunParameters p;
			p.massKg = 1000;
			p.dustFraction = e.dust;
			p.largestFraction = e.largest;
			p.beta = e.beta;
			PhysicsConfig physics;
			physics.individualFragments = {1000, 1000};
			const FragmentationModel model(physics);
			const auto groups = model.groups(parent, p);
			ASSERT_EQ(groups.size(), e.count);
			for (std::size_t i = 0; i < e.first.size(); ++i)
				EXPECT_NEAR(groups[i].massKg, e.first[i], 1e-10);
			EXPECT_NEAR(groups.back().massKg, e.remainder, 1e-10);
			double total = 0;
			for (const auto& g : groups) {
				EXPECT_EQ(g.multiplicity, 1);
				EXPECT_GE(g.massKg, physics.minMassKg);
				total += g.massKg;
			}
			EXPECT_NEAR(total, 1000 * (1 - e.dust), 1e-10);
			const auto batch =
			    model.split(parent, p, standardAtmosphere(30000), groups, std::vector<double>(groups.size(), .5));
			EXPECT_NEAR(batch.residualKg, 0, 1e-9);
			EXPECT_NEAR(batch.dustKg, 7000 * e.dust, 1e-10);
			for (std::size_t reps : {1U, 3U, 6U}) {
				PhysicsConfig grouped;
				grouped.representativesPerBin = reps;
				double mass = 0, count = 0;
				for (const auto& g : FragmentationModel(grouped).groups(parent, p)) {
					mass += g.massKg * g.multiplicity;
					count += g.multiplicity;
				}
				EXPECT_NEAR(mass, total, 1e-10);
				EXPECT_NEAR(count, static_cast<double>(e.count), 1e-10);
			}
			physics.retainMassRemainder = false;
			const FragmentationModel legacy(physics);
			const auto oldGroups = legacy.groups(parent, p);
			ASSERT_EQ(oldGroups.size(), e.count - 1);
			const auto oldBatch = legacy.split(parent, p, standardAtmosphere(30000), oldGroups,
			                                   std::vector<double>(oldGroups.size(), .5));
			EXPECT_NEAR(oldBatch.residualKg, 7 * e.remainder, 1e-9);
		}
	}

	TEST(Fragmentation, CutoffDoesNotInventResolvedFragmentsFromUnresolvedTail)
	{
		const auto check = [](double parentMass, double dust, double largest, double beta, double cutoff,
		                      std::size_t expectedCount) {
			PhysicsConfig physics;
			physics.minMassKg = cutoff;
			physics.individualFragments = {1000, 1000};
			RunParameters p;
			p.dustFraction = dust;
			p.largestFraction = largest;
			p.beta = beta;
			const auto parent = fragment(30000, parentMass);
			const FragmentationModel model(physics);
			const auto groups = model.groups(parent, p);
			ASSERT_EQ(groups.size(), expectedCount);
			double mass = 0;
			for (const auto& g : groups) {
				EXPECT_GE(g.massKg, cutoff);
				mass += g.massKg;
			}
			const auto batch =
			    model.split(parent, p, standardAtmosphere(30000), groups, std::vector<double>(groups.size(), .5));
			EXPECT_GT(batch.residualKg, 0);
			EXPECT_NEAR(mass + batch.residualKg + batch.dustKg, parentMass, 1e-10);
		};
		check(1000, .2, .45, .35, 3, 5);
		check(.01, 0, .15, .55, .001, 2);
		check(.0015, 0, .1, .55, .001, 0);
	}

	TEST(Fragmentation, RemainderInheritsStrengthMultiplicityAndGeneration)
	{
		PhysicsConfig physics;
		physics.individualFragments = {1000, 1000};
		physics.strengthCapPa = 12e6;
		RunParameters p;
		p.dustFraction = .2;
		p.largestFraction = .45;
		p.beta = .35;
		auto parent = fragment(30000, 1000);
		parent.id = 17;
		parent.generation = 3;
		parent.multiplicity = 7;
		const FragmentationModel model(physics);
		const auto groups = model.groups(parent, p);
		const auto batch =
		    model.split(parent, p, standardAtmosphere(30000), groups, std::vector<double>(groups.size(), .5));
		ASSERT_EQ(batch.children.size(), 6U);
		for (const auto& child : batch.children) {
			EXPECT_EQ(child.parentId, parent.id);
			EXPECT_EQ(child.generation, 4U);
			EXPECT_EQ(child.multiplicity, 7);
			EXPECT_DOUBLE_EQ(
			    child.strengthPa,
			    std::min(physics.strengthCapPa,
			             parent.strengthPa * std::pow(parent.state.massKg / child.state.massKg, physics.weibullAlpha)));
		}
		EXPECT_DOUBLE_EQ(batch.children.back().strengthPa, physics.strengthCapPa);
	}

	TEST(Integrator, DefaultRamPressureAndExplicitConventionConversion)
	{
		const auto env = Environment::idealized();
		auto f = fragment(100000, 1000);
		f.state.velocityEcefMps = enuToEcef({9000, 0, -9000}, 49.8 * deg, .8 * deg);
		f.strengthPa = 1e6;
		PhysicsConfig physics;
		ASSERT_EQ(physics.ramPressure, RamPressure::RhoV2);
		const auto full = TrajectoryIntegrator(env, physics).advance(f, {});
		ASSERT_EQ(full.event, FlightEvent::Breakup);
		const auto position = geodetic(full.state.positionEcefM);
		double v2 = 0;
		for (double v : full.state.velocityEcefMps)
			v2 += v * v;
		EXPECT_NEAR(env.air(position).densityKgM3 * v2, f.strengthPa, 1);
		physics.ramPressure = RamPressure::HalfRhoV2;
		f.strengthPa *= .5;
		const auto half = TrajectoryIntegrator(env, physics).advance(f, {});
		EXPECT_EQ(half.event, full.event);
		EXPECT_EQ(half.state.positionEcefM, full.state.positionEcefM);
		EXPECT_DOUBLE_EQ(half.state.massKg, full.state.massKg);
		EXPECT_DOUBLE_EQ(half.timeS, full.timeS);
	}

	TEST(Fragmentation, ConservesMomentumAndWorksAtPole)
	{
		const auto env = Environment::idealized();
		RunParameters p;
		p.massKg = 100;
		auto parent = fragment(30000, 100);
		parent.state.velocityEcefMps = {0, 0, -12000};
		parent.multiplicity = 7;
		FragmentationModel model({});
		const auto groups = model.groups(parent, p);
		std::vector<double> draws(groups.size());
		for (std::size_t i = 0; i < draws.size(); ++i)
			draws[i] = static_cast<double>(i) / static_cast<double>(draws.size());
		const auto batch = model.split(parent, p, env.air(geodetic(parent.state.positionEcefM)), groups, draws);
		Vector3 momentum{};
		double mass = 0;
		for (const auto& f : batch.children) {
			const double weight = f.state.massKg * f.multiplicity;
			mass += weight;
			for (std::size_t c = 0; c < 3; ++c)
				momentum[c] += weight * f.state.velocityEcefMps[c];
		}
		for (std::size_t c = 0; c < 3; ++c)
			EXPECT_NEAR(momentum[c] / mass, parent.state.velocityEcefMps[c], 1e-10);
		EXPECT_NEAR(mass + batch.dustKg + batch.residualKg, 700, 1e-10);
	}

	TEST(Fragmentation, RepresentativesPreserveBinMassMultiplicityAndMomentum)
	{
		RunParameters p;
		p.massKg = 100;
		const auto parent = fragment(30000, 100);
		const auto original = FragmentationModel({}).groups(parent, p);
		for (std::size_t reps : {1U, 3U, 6U}) {
			PhysicsConfig physics;
			physics.representativesPerBin = reps;
			FragmentationModel model(physics);
			const auto groups = model.groups(parent, p);
			std::size_t j = 0;
			for (std::size_t i = 0; i < original.size(); ++i) {
				const auto n = i < physics.individualFragments[0]
				                   ? 1U
				                   : std::min(reps, static_cast<std::size_t>(original[i].multiplicity));
				double multiplicity = 0;
				for (std::size_t k = 0; k < n; ++k) {
					ASSERT_LT(j, groups.size());
					EXPECT_DOUBLE_EQ(groups[j].massKg, original[i].massKg);
					EXPECT_GE(groups[j].multiplicity, 1);
					multiplicity += groups[j++].multiplicity;
				}
				EXPECT_NEAR(multiplicity, original[i].multiplicity, 1e-12);
			}
			EXPECT_EQ(j, groups.size());
			std::vector<double> draws(groups.size());
			for (std::size_t i = 0; i < draws.size(); ++i)
				draws[i] = static_cast<double>(i) / draws.size();
			const auto batch = model.split(parent, p, standardAtmosphere(30000), groups, draws);
			double mass = 0;
			Vector3 momentum{};
			for (const auto& f : batch.children) {
				const double weight = f.state.massKg * f.multiplicity;
				mass += weight;
				for (std::size_t c = 0; c < 3; ++c)
					momentum[c] += weight * f.state.velocityEcefMps[c];
			}
			EXPECT_NEAR(mass + batch.dustKg + batch.residualKg, parent.state.massKg, 1e-10);
			for (std::size_t c = 0; c < 3; ++c)
				EXPECT_NEAR(momentum[c] / mass, parent.state.velocityEcefMps[c], 1e-9);
		}
	}

	TEST(Fragmentation, RejectsExcessMassAndInvalidParent)
	{
		const auto env = Environment::idealized();
		const FragmentationModel model({});
		auto parent = fragment(30000, 100);
		const std::vector<FragmentGroup> groups{{101, 1}};
		const std::vector<double> draws{0.5};
		EXPECT_THROW(model.split(parent, {}, env.air(geodetic(parent.state.positionEcefM)), groups, draws),
		             std::invalid_argument);
		parent.multiplicity = -1;
		EXPECT_THROW(model.groups(parent, {}), std::invalid_argument);
	}

	TEST(Simulator, FullCascadeBalancesMassAndReportsLimits)
	{
		const auto env = Environment::idealized();
		auto initial = entry(30000);
		initial.velocityEcefMps = enuToEcef({9000, 0, -9000}, 49.8 * deg, 0.8 * deg);
		RunParameters p;
		p.massKg = 500;
		p.strengthPa = 2e6;
		PhysicsConfig physics;
		const auto full = Simulator(env, physics).runOne(initial, p, 1);
		EXPECT_EQ(full.status, RunStatus::Complete) << full.error;
		EXPECT_GT(full.breakups.size(), 0U);
		EXPECT_GT(full.impacts.size(), 0U);
		balanced(full.mass);
		physics.maxParticles = 1;
		const auto limited = Simulator(env, physics).runOne(initial, p, 1);
		EXPECT_EQ(limited.status, RunStatus::ParticleLimit);
		EXPECT_GT(limited.mass.unresolvedKg, 0);
		balanced(limited.mass);
	}

	TEST(Simulator, ResidualFragmentReachesGround)
	{
		const auto env = Environment::idealized();
		auto initial = entry(1000);
		initial.velocityEcefMps = enuToEcef({0, 0, -12000}, 49.8 * deg, .8 * deg);
		RunParameters p;
		p.massKg = 1000;
		p.strengthPa = 10000;
		p.ablationKgJ = 0;
		p.dustFraction = .2;
		p.largestFraction = .45;
		p.beta = .35;
		PhysicsConfig physics;
		physics.individualFragments = {1000, 1000};
		physics.weibullAlpha = 10;
		physics.strengthCapPa = 1e20;
		const auto result = Simulator(env, physics).runOne(initial, p, 42);
		ASSERT_EQ(result.status, RunStatus::Complete) << result.error;
		ASSERT_EQ(result.breakups.size(), 1U);
		ASSERT_EQ(result.impacts.size(), 6U);
		EXPECT_NEAR(result.mass.landedKg, 800, 1e-10);
		EXPECT_NEAR(result.mass.dustKg, 200, 1e-10);
		EXPECT_NEAR(result.mass.residualKg, 0, 1e-10);
		const auto smallest = std::min_element(
		    result.impacts.begin(), result.impacts.end(),
		    [](const auto& a, const auto& b) { return a.fragment.state.massKg < b.fragment.state.massKg; });
		EXPECT_NEAR(smallest->fragment.state.massKg, 2.407719106529015, 1e-10);
		EXPECT_GT(smallest->fragment.timeS, 0);
		balanced(result.mass);
	}

	TEST(Simulator, RemainderCascadesAreDeterministicAcrossWorkers)
	{
		const auto env = Environment::idealized();
		auto initial = entry(100000);
		initial.velocityEcefMps = enuToEcef({9000, 0, -9000}, 49.8 * deg, .8 * deg);
		const std::vector<EntryState> entries{initial};
		PhysicsConfig physics;
		physics.representativesPerBin = 3;
		MonteCarloConfig mc;
		mc.runs = 24;
		mc.workers = 1;
		mc.parameters.size = MassSize{{500, 500}};
		mc.parameters.strongStrengthPa = {.3e6, .3e6};
		std::vector<RunResult> serial;
		const Simulator simulator(env, physics);
		const auto first = simulator.run(entries, mc, {}, [&](const auto& run) { serial.push_back(run); });
		ASSERT_EQ(first.failed, 0U);
		ASSERT_EQ(first.truncated, 0U);
		ASSERT_EQ(serial.size(), mc.runs);
		bool secondary = false;
		for (const auto& run : serial) {
			balanced(run.mass);
			for (const auto& breakup : run.breakups)
				secondary |= breakup.parent.generation > 0;
		}
		EXPECT_TRUE(secondary);
		mc.workers = 15;
		std::size_t index = 0;
		const auto second = simulator.run(entries, mc, {}, [&](const auto& run) {
			ASSERT_LT(index, serial.size());
			const auto& a = serial[index++];
			ASSERT_EQ(run.runId, a.runId);
			ASSERT_EQ(run.impacts.size(), a.impacts.size());
			ASSERT_EQ(run.breakups.size(), a.breakups.size());
			EXPECT_EQ(run.integration.evaluations, a.integration.evaluations);
			EXPECT_EQ(run.mass.residualKg, a.mass.residualKg);
			EXPECT_EQ(run.mass.landedKg, a.mass.landedKg);
			for (std::size_t i = 0; i < run.impacts.size(); ++i) {
				const auto& x = run.impacts[i].fragment;
				const auto& y = a.impacts[i].fragment;
				EXPECT_EQ(x.id, y.id);
				EXPECT_EQ(x.state.positionEcefM, y.state.positionEcefM);
				EXPECT_EQ(x.state.massKg, y.state.massKg);
				EXPECT_EQ(x.multiplicity, y.multiplicity);
			}
		});
		EXPECT_EQ(second.failed, 0U);
		EXPECT_EQ(second.truncated, 0U);
		EXPECT_EQ(index, serial.size());
	}

	TEST(Simulator, UnlimitedParticleCountPreservesCompletedCascade)
	{
		const auto env = Environment::idealized();
		auto initial = entry(30000);
		initial.velocityEcefMps = enuToEcef({9000, 0, -9000}, 49.8 * deg, 0.8 * deg);
		RunParameters p;
		p.massKg = 500;
		p.strengthPa = 2e6;
		PhysicsConfig physics;
		ASSERT_EQ(physics.maxParticles, 0U);
		const auto unlimited = Simulator(env, physics).runOne(initial, p, 1);
		ASSERT_EQ(unlimited.status, RunStatus::Complete);
		ASSERT_GT(unlimited.particles, 1U);
		physics.maxParticles = unlimited.particles;
		const auto limited = Simulator(env, physics).runOne(initial, p, 1);
		ASSERT_EQ(limited.status, RunStatus::Complete);
		ASSERT_EQ(limited.impacts.size(), unlimited.impacts.size());
		EXPECT_EQ(limited.integration.evaluations, unlimited.integration.evaluations);
		for (std::size_t i = 0; i < unlimited.impacts.size(); ++i) {
			EXPECT_EQ(limited.impacts[i].fragment.state.positionEcefM,
			          unlimited.impacts[i].fragment.state.positionEcefM);
			EXPECT_EQ(limited.impacts[i].fragment.state.massKg, unlimited.impacts[i].fragment.state.massKg);
		}
		balanced(unlimited.mass);
	}

	TEST(Integrator, WallDeadlineInterruptsAcceptedStepsWithoutChangingFlightHorizon)
	{
		const auto env = Environment::idealized();
		bool delayed = false;
		const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(100);
		const auto result = TrajectoryIntegrator(env, dark())
		                        .advance(
		                            fragment(1000), {}, {},
		                            [&](double time, const auto&) {
			                            if (time > 0 && !delayed) {
				                            delayed = true;
				                            std::this_thread::sleep_until(deadline + std::chrono::milliseconds(1));
			                            }
		                            },
		                            deadline);
		EXPECT_TRUE(delayed);
		EXPECT_EQ(result.event, FlightEvent::WallTimeout);
		EXPECT_GT(result.timeS, 0);
		EXPECT_LT(result.timeS, IntegrationOptions{}.maxFlightTimeS);
		EXPECT_GT(result.stats.acceptedSteps, 0U);
	}

	TEST(Fragmentation, GroupEnumerationAndSplitHonorDeadlineAndCancellation)
	{
		const auto env = Environment::idealized();
		PhysicsConfig physics;
		physics.minMassKg = 1e-30;
		const FragmentationModel model(physics);
		auto parent = fragment(30000, 1e12);
		RunParameters parameters;
		parameters.largestFraction = 1e-8;
		const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(1);
		try {
			model.groups(parent, parameters, {}, deadline);
			FAIL() << "Expected a wall-clock interruption during enumeration";
		} catch (const detail::RunInterrupted& error) {
			EXPECT_EQ(error.status(), RunStatus::WallTimeout);
		}
		std::stop_source stop;
		stop.request_stop();
		try {
			model.groups(parent, {}, stop.get_token(), deadline);
			FAIL() << "Expected cancellation to take priority";
		} catch (const detail::RunInterrupted& error) {
			EXPECT_EQ(error.status(), RunStatus::Cancelled);
		}
		const std::vector<FragmentGroup> groups{{1, 1}};
		const std::vector<double> draws{0.5};
		EXPECT_THROW(model.split(parent, {}, env.air(geodetic(parent.state.positionEcefM)), groups, draws, {},
		                         std::chrono::steady_clock::now()),
		             detail::RunInterrupted);
	}

	TEST(Simulator, TimeoutDuringTracePreservesMassAndIsNotUserCancellation)
	{
		const auto env = Environment::idealized();
		const auto initial = entry(1000);
		bool delayed = false;
		const auto result = Simulator(env, dark())
		                        .runOne(
		                            initial, {}, 42, 0, {},
		                            [&](const auto&, double time, const auto&) {
			                            if (time > 0 && !delayed) {
				                            delayed = true;
				                            std::this_thread::sleep_for(std::chrono::milliseconds(110));
			                            }
		                            },
		                            .1);
		ASSERT_TRUE(delayed);
		EXPECT_EQ(result.status, RunStatus::WallTimeout);
		EXPECT_GT(result.mass.unresolvedKg, 0);
		balanced(result.mass);
		EXPECT_FALSE(result.error.empty());
	}

	TEST(Simulator, RealizationDeadlinesRemainFailuresInTheOrderedEnsemble)
	{
		const auto env = Environment::idealized();
		const std::vector<EntryState> entries{entry()};
		MonteCarloConfig config;
		config.runs = 5;
		config.workers = 3;
		config.runTimeoutS = 1e-12;
		std::vector<std::uint64_t> ids;
		const auto summary = Simulator(env, dark()).run(entries, config, {}, [&](const auto& run) {
			EXPECT_EQ(run.status, RunStatus::WallTimeout);
			balanced(run.mass);
			ids.push_back(run.runId);
		});
		EXPECT_EQ(ids, (std::vector<std::uint64_t>{0, 1, 2, 3, 4}));
		EXPECT_EQ(summary.completed, 5U);
		EXPECT_EQ(summary.failed, 5U);
		EXPECT_EQ(summary.wallTimeouts, 5U);
		EXPECT_EQ(summary.truncated, 0U);
		EXPECT_FALSE(summary.cancelled);
		EXPECT_FALSE(summary.converged);
	}

	TEST(Simulator, ThreadCountDoesNotChangeResults)
	{
		const auto env = Environment::idealized();
		const Simulator simulator(env, dark());
		const std::vector<EntryState> entries{entry(1000)};
		MonteCarloConfig options;
		options.runs = 24;
		options.workers = 1;
		std::vector<RunResult> serial, parallel;
		simulator.run(entries, options, {}, [&](const auto& r) { serial.push_back(r); });
		options.workers = 15;
		const auto callingThread = std::this_thread::get_id();
		const auto summary = simulator.run(entries, options, {}, [&](const auto& r) {
			EXPECT_EQ(std::this_thread::get_id(), callingThread);
			parallel.push_back(r);
		});
		ASSERT_EQ(serial.size(), parallel.size());
		EXPECT_EQ(summary.workers, 15U);
		EXPECT_EQ(summary.completed, options.runs);
		for (std::size_t i = 0; i < serial.size(); ++i) {
			ASSERT_EQ(serial[i].impacts.size(), 1U);
			ASSERT_EQ(parallel[i].impacts.size(), 1U);
			EXPECT_EQ(serial[i].runId, i);
			EXPECT_EQ(serial[i].impacts[0].fragment.state.positionEcefM,
			          parallel[i].impacts[0].fragment.state.positionEcefM);
			EXPECT_EQ(serial[i].parameters.densityKgM3, parallel[i].parameters.densityKgM3);
		}
	}

	TEST(Simulator, DiagnosticObserverPreservesCascadeExactly)
	{
		const auto env = Environment::idealized();
		auto initial = entry(30000);
		initial.velocityEcefMps = enuToEcef({9000, 0, -9000}, 49.8 * deg, .8 * deg);
		RunParameters p;
		p.massKg = 500;
		p.strengthPa = 2e6;
		const Simulator simulator(env, {});
		const auto plain = simulator.runOne(initial, p, 17, 3);
		std::size_t samples = 0;
		std::map<std::uint64_t, std::pair<double, FlightState>> last;
		const auto observed = simulator.runOne(initial, p, 17, 3, {},
		                                       [&](const Fragment& fragment, double time, const FlightState& state) {
			                                       ++samples;
			                                       EXPECT_GE(time, fragment.timeS);
			                                       EXPECT_LE(state.massKg, fragment.state.massKg);
			                                       last[fragment.id] = {time, state};
		                                       });
		ASSERT_EQ(plain.status, RunStatus::Complete);
		ASSERT_EQ(observed.status, plain.status);
		ASSERT_EQ(observed.impacts.size(), plain.impacts.size());
		ASSERT_EQ(observed.breakups.size(), plain.breakups.size());
		EXPECT_GT(samples, observed.impacts.size());
		EXPECT_EQ(observed.integration.acceptedSteps, plain.integration.acceptedSteps);
		EXPECT_EQ(observed.integration.rejectedSteps, plain.integration.rejectedSteps);
		EXPECT_EQ(observed.integration.evaluations, plain.integration.evaluations);
		for (std::size_t i = 0; i < plain.impacts.size(); ++i) {
			const auto& a = plain.impacts[i].fragment;
			const auto& b = observed.impacts[i].fragment;
			EXPECT_EQ(a.id, b.id);
			EXPECT_EQ(a.parentId, b.parentId);
			EXPECT_EQ(a.state.positionEcefM, b.state.positionEcefM);
			EXPECT_EQ(a.state.velocityEcefMps, b.state.velocityEcefMps);
			EXPECT_EQ(a.state.massKg, b.state.massKg);
			EXPECT_EQ(a.timeS, b.timeS);
			EXPECT_EQ(a.multiplicity, b.multiplicity);
		}
		EXPECT_EQ(observed.mass.ablatedKg, plain.mass.ablatedKg);
		EXPECT_EQ(observed.mass.landedKg, plain.mass.landedKg);
		for (const auto& event : observed.breakups) {
			ASSERT_TRUE(last.contains(event.parent.id));
			const auto& [time, state] = last.at(event.parent.id);
			EXPECT_EQ(time, event.parent.timeS);
			EXPECT_EQ(state.massKg, event.parent.state.massKg);
			EXPECT_EQ(state.positionEcefM, event.parent.state.positionEcefM);
		}
	}

	TEST(Simulator, CallbackCancellationPreservesOrderedPrefix)
	{
		const auto env = Environment::idealized();
		std::stop_source stop;
		RunControl control;
		control.cancel = stop.get_token();
		control.progress = [&](const auto& p) {
			if (p.completed == 3)
				stop.request_stop();
		};
		MonteCarloConfig mc;
		mc.runs = 50;
		std::vector<std::uint64_t> ids;
		const std::vector<EntryState> entries{entry(100)};
		const auto result =
		    Simulator(env, dark()).run(entries, mc, control, [&](const auto& r) { ids.push_back(r.runId); });
		EXPECT_TRUE(result.cancelled);
		EXPECT_EQ(result.completed, 3U);
		EXPECT_EQ(ids, (std::vector<std::uint64_t>{0, 1, 2}));
	}

	TEST(Simulator, CallbackExceptionStopsWorkers)
	{
		const auto env = Environment::idealized();
		const std::vector<EntryState> entries{entry(100)};
		MonteCarloConfig mc;
		mc.runs = 50;
		EXPECT_THROW(
		    Simulator(env, dark()).run(entries, mc, {}, [](const auto&) { throw std::runtime_error("writer failed"); }),
		    std::runtime_error);
	}

	TEST(Simulator, ZeroWeightEntryIsNeverChosen)
	{
		const auto env = Environment::idealized();
		auto zero = entry(100);
		zero.weight = 0;
		zero.id = 1;
		auto valid = entry(100);
		valid.id = 2;
		const std::vector<EntryState> entries{zero, valid};
		MonteCarloConfig mc;
		mc.runs = 30;
		Simulator(env, dark()).run(entries, mc, {}, [](const auto& r) { EXPECT_EQ(r.entry.id, 2U); });
	}

	TEST(Simulator, InvalidConfigurationFailsBeforeWorkers)
	{
		const auto env = Environment::idealized();
		const std::vector<EntryState> entries{entry()};
		MonteCarloConfig mc;
		mc.parameters.beta = {0.8, 0.4};
		EXPECT_THROW(Simulator(env).run(entries, mc), std::invalid_argument);
		mc = {};
		mc.parameters.densityKgM3 = {-1, -1};
		EXPECT_THROW(Simulator(env).run(entries, mc), std::invalid_argument);
		mc = {};
		mc.workers = 0;
		EXPECT_THROW(Simulator(env).run(entries, mc), std::invalid_argument);
		for (double timeout :
		     {0.0, -1.0, 86401.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
			mc = {};
			mc.runTimeoutS = timeout;
			EXPECT_THROW(Simulator(env).run(entries, mc), std::invalid_argument);
			EXPECT_THROW(Simulator(env).runOne(entries.front(), {}, 42, 0, {}, timeout), std::invalid_argument);
		}
	}

	TEST(Regression, PythonSingleFragmentFixtures)
	{
		std::ifstream stream(std::string(STREWN_TEST_DATA) + "/python_flight.json");
		ASSERT_TRUE(stream.good());
		const auto fixture = nlohmann::json::parse(stream);
		AtmosphereGrid a;
		a.heightCount = fixture["atmosphere"].size();
		for (const auto& row : fixture["atmosphere"])
			a.samples.push_back(
			    {row[0].get<double>(), row[1].get<double>(), {row[2].get<double>(), row[3].get<double>(), 0}});
		const Environment env(std::move(a), {});
		for (const auto& test : fixture["cases"]) {
			SCOPED_TRACE(test["name"].get<std::string>());
			Fragment f;
			const auto state = test["initial"].get<std::array<double, 7>>();
			f.state = {{state[0], state[1], state[2]}, {state[3], state[4], state[5]}, state[6]};
			f.strengthPa = test["strength_pa"].get<double>();
			PhysicsConfig physics;
			physics.phase = test["dark"].get<bool>() ? FlightPhase::DarkFlight : FlightPhase::Entry;
			const auto result = TrajectoryIntegrator(env, physics).advance(f, RunParameters{});
			const auto expected = test["final"].get<std::array<double, 7>>();
			IntegrationOptions tight;
			tight.relativeTolerance = 1e-11;
			tight.positionToleranceM = 1e-4;
			tight.velocityToleranceMps = 1e-6;
			tight.massToleranceKg = 1e-12;
			tight.maxStepS = 0.05;
			const auto refined = TrajectoryIntegrator(env, physics, tight).advance(f, RunParameters{});
			const auto reference = test["refined_final"].get<std::array<double, 7>>();
			std::cout << test["name"] << std::setprecision(9) << " position_difference_m="
			          << distance(result.state.positionEcefM, {expected[0], expected[1], expected[2]})
			          << " mass_difference_kg=" << result.state.massKg - expected[6]
			          << " refined_mass_difference_kg=" << refined.state.massKg - reference[6] << '\n';
			EXPECT_EQ(static_cast<int>(result.event), test["event"].get<int>());
			EXPECT_LT(distance(result.state.positionEcefM, {expected[0], expected[1], expected[2]}), 0.1);
			EXPECT_LT(distance(result.state.velocityEcefMps, {expected[3], expected[4], expected[5]}), 0.005);
			EXPECT_NEAR(result.state.massKg, reference[6], std::max(1e-8, reference[6] * 4e-6));
			EXPECT_NEAR(expected[6], reference[6], std::max(1e-8, reference[6] * 4e-6));
			EXPECT_EQ(refined.event, result.event);
			EXPECT_NEAR(refined.state.massKg, reference[6], std::max(1e-10, reference[6] * 5e-8));
			EXPECT_LT(distance(refined.state.positionEcefM, {reference[0], reference[1], reference[2]}), 0.001);
			EXPECT_NEAR(refined.timeS, test["refined_time_s"].get<double>(), 0.0001);
			EXPECT_NEAR(result.timeS, test["time_s"].get<double>(), 0.005);
		}
	}

	TEST(Regression, PythonFragmentGroupingAndKicks)
	{
		std::ifstream stream(std::string(STREWN_TEST_DATA) + "/python_flight.json");
		const auto fixture = nlohmann::json::parse(stream);
		const auto& test = fixture["fragmentation"];
		AtmosphereGrid a;
		a.heightCount = fixture["atmosphere"].size();
		for (const auto& row : fixture["atmosphere"])
			a.samples.push_back(
			    {row[0].get<double>(), row[1].get<double>(), {row[2].get<double>(), row[3].get<double>(), 0}});
		const Environment env(std::move(a), {});
		const auto state = test["parent"].get<std::array<double, 7>>();
		Fragment parent;
		parent.state = {{state[0], state[1], state[2]}, {state[3], state[4], state[5]}, state[6]};
		parent.strengthPa = 4e6;
		parent.multiplicity = 7;
		PhysicsConfig physics;
		physics.retainMassRemainder = false;
		FragmentationModel model(physics);
		const auto groups = model.groups(parent, {});
		ASSERT_EQ(groups.size(), test["groups"].size());
		for (std::size_t i = 0; i < groups.size(); ++i) {
			EXPECT_NEAR(groups[i].massKg, test["groups"][i][0].get<double>(), 1e-12);
			EXPECT_EQ(groups[i].multiplicity, test["groups"][i][1].get<double>());
		}
		const auto draws = test["draws"].get<std::vector<double>>();
		const auto batch = model.split(parent, {}, env.air(geodetic(parent.state.positionEcefM)), groups, draws);
		ASSERT_EQ(batch.children.size(), test["children"].size());
		for (std::size_t i = 0; i < batch.children.size(); ++i) {
			const auto& expected = test["children"][i];
			const auto v = expected["state"].get<std::array<double, 7>>();
			EXPECT_LT(distance(batch.children[i].state.velocityEcefMps, {v[3], v[4], v[5]}), 1e-8);
			EXPECT_NEAR(batch.children[i].strengthPa, expected["strength"].get<double>(), 1e-5);
			EXPECT_EQ(batch.children[i].multiplicity, expected["mult"].get<double>());
		}
	}

}  // namespace
