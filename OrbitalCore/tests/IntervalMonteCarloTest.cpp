#include <gtest/gtest.h>
#include "orbital_core/hazard/IntervalMonteCarlo.h"
#include "hazard/McBoundaryScan.h"
#include "hazard/BoundaryScan.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/hazard/SurfaceGeometry.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include <nlohmann/json.hpp>
#include <fstream>

namespace {
	namespace hz = orbital_core::hazard;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;

	class FixedEarthMc : public IEphemerisSource {
	   public:
		std::optional<LunarLibration> libration;

		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, {1, 0, 0}, {}}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "interval-mc-test"; }

		std::optional<LunarLibration> lunarLibrationMjdTdb(double) override { return libration; }
	};

	class IntervalMonteCarloTest : public ::testing::Test {
	   protected:
		FixedEarthMc eph;
		IEphemerisSource* previous{};
		orbital_core::FitPipelineResult source;
		hz::IntervalMcConfig config;

		void SetUp() override
		{
			previous = getDefaultEphemerisSource();
			setDefaultEphemerisSource(&eph);
			source.fit.converged = true;
			source.fit.fitEpochOrbit = std::make_shared<Orbit>();
			source.fit.fitEpochOrbit->setEpoch(2460000.5);
			source.fit.fitEpochOrbit->setR({1.0 + 7000.0 / auKm, 0, 0});
			source.fit.fitEpochOrbit->setV({-10.0 * 86400.0 / auKm / Constants::GAUSS_K, 0, 0});
			source.fit.uncertainty.fitParameterCovariance = Eigen::MatrixXd::Identity(6, 6) * 1e-22;
			source.effectiveConfig.nbody.enabled = true;
			source.effectiveConfig.nbody.planetMask = 0;
			source.effectiveConfig.nbody.useSb441 = false;
			source.effectiveConfig.nbody.stepSize = .001;
			config.startMjdTdb = 60000;
			config.endMjdTdb = 60000.002;
			config.samples = 8;
			config.workers = 1;
		}

		void TearDown() override { setDefaultEphemerisSource(previous); }
	};

	TEST_F(IntervalMonteCarloTest, BoundaryAndParallelReproducibility)
	{
		const auto first = hz::runIntervalMonteCarlo(source, config);
		ASSERT_EQ(first.completed, 8);
		ASSERT_EQ(first.events.size(), 1);
		EXPECT_EQ(first.entryClones, 8);
		EXPECT_EQ(first.events[0].entryClones, 8);
		EXPECT_EQ(first.nominal.status, hz::McCloneStatus::Complete);
		EXPECT_DOUBLE_EQ(first.nominal.state[3], -10.0 * 86400.0 / auKm);
		config.workers = 3;
		const auto second = hz::runIntervalMonteCarlo(source, config);
		ASSERT_EQ(second.completed, 8);
		for (int i = 0; i < 8; ++i) {
			EXPECT_EQ(first.clones[i].state, second.clones[i].state);
			ASSERT_TRUE(first.clones[i].entry.has_value());
			const auto& entry = *first.clones[i].entry;
			EXPECT_DOUBLE_EQ(entry.mjdTdb, second.clones[i].entry->mjdTdb);
			EXPECT_NEAR(entry.relativePositionKm.magnitude(), Constants::EARTH_EQ_RADIUS_KM + 100, .05);
			EXPECT_LT(entry.mjdTdb, config.endMjdTdb);
		}
	}

	TEST_F(IntervalMonteCarloTest, RejectsIntervalWithoutEarthOrientation)
	{
		config.startMjdTdb = 234166.0;
		config.endMjdTdb = 234166.002;
		EXPECT_THROW((void)hz::runIntervalMonteCarlo(source, config), std::runtime_error);
	}

	TEST_F(IntervalMonteCarloTest, SharedBoundarySearchMatchesNominal)
	{
		const auto result = hz::runIntervalMonteCarlo(source, config);
		ASSERT_TRUE(result.nominal.entry);
		NBodyPropagator propagator(config.integrator);
		propagator.setStepSize(source.effectiveConfig.nbody.stepSize);
		propagator.setPlanetMask(source.effectiveConfig.nbody.planetMask);
		propagator = propagator.forEncounterMonitoring(eph);
		const auto crossing = hz::boundary::search(result.nominal.state, result.nominal.epochMjdTdb,
		                                           result.nominal.nonGrav, config.target, config.startMjdTdb,
		                                           config.endMjdTdb, config.altitudeKm, propagator, eph, {});
		ASSERT_TRUE(crossing.entry);
		EXPECT_DOUBLE_EQ(crossing.entry->mjdTdb, result.nominal.entry->mjdTdb);
		EXPECT_DOUBLE_EQ(crossing.entry->latitudeDeg, result.nominal.entry->latitudeDeg);
		EXPECT_DOUBLE_EQ(crossing.entry->longitudeDeg, result.nominal.entry->longitudeDeg);
	}

	TEST_F(IntervalMonteCarloTest, NearbyNominalEntryMatchesMonteCarloAtSameBoundary)
	{
		source.fit.orbit = source.fit.fitEpochOrbit;
		config.altitudeKm = 50;
		const auto entry = hz::boundary::nearbyNominalEarthEntry(source.fit, source.effectiveConfig.nbody, {});
		ASSERT_TRUE(entry);
		const auto mc = hz::runIntervalMonteCarlo(source, config);
		ASSERT_TRUE(mc.nominal.entry);
		EXPECT_DOUBLE_EQ(entry->altitudeKm, 50);
		EXPECT_DOUBLE_EQ(entry->crossing.mjdTdb, mc.nominal.entry->mjdTdb);
		EXPECT_DOUBLE_EQ(entry->crossing.latitudeDeg, mc.nominal.entry->latitudeDeg);
		EXPECT_DOUBLE_EQ(entry->crossing.longitudeDeg, mc.nominal.entry->longitudeDeg);
	}

	TEST_F(IntervalMonteCarloTest, NearbyNominalEntryDoesNotReportOutboundOrUnconvergedOrbit)
	{
		source.fit.orbit = source.fit.fitEpochOrbit;
		source.fit.orbit->setV({10.0 * 86400.0 / auKm / Constants::GAUSS_K, 0, 0});
		EXPECT_FALSE(hz::boundary::nearbyNominalEarthEntry(source.fit, source.effectiveConfig.nbody, {}));
		source.fit.converged = false;
		EXPECT_FALSE(hz::boundary::nearbyNominalEarthEntry(source.fit, source.effectiveConfig.nbody, {}));
	}

	TEST_F(IntervalMonteCarloTest, LunarOrientationChangesCoordinatesNotPhysicalEvents)
	{
		config.target = hz::Planet::Moon;
		source.fit.fitEpochOrbit->setR({1.0 + 2000.0 / auKm, 0, 0});
		EXPECT_THROW((void)hz::runIntervalMonteCarlo(source, config), std::runtime_error);
		eph.libration = LunarLibration{{.1, .4, .7}, 440};
		const auto first = hz::runIntervalMonteCarlo(source, config);
		ASSERT_EQ(first.completed, 8);
		ASSERT_EQ(first.entryClones, 8);
		EXPECT_EQ(first.bodyFixedFrame, "MOON_ME_DE440_ME421");
		eph.libration->angles.z += .2;
		config.workers = 3;
		const auto second = hz::runIntervalMonteCarlo(source, config);
		ASSERT_EQ(second.completed, 8);
		ASSERT_EQ(second.entryClones, 8);
		for (int i = 0; i < 8; ++i) {
			ASSERT_TRUE(first.clones[i].entry);
			ASSERT_TRUE(second.clones[i].entry);
			const auto& a = *first.clones[i].entry;
			const auto& b = *second.clones[i].entry;
			EXPECT_DOUBLE_EQ(a.mjdTdb, b.mjdTdb);
			EXPECT_DOUBLE_EQ(a.relativePositionKm.x, b.relativePositionKm.x);
			EXPECT_DOUBLE_EQ(a.relativePositionKm.y, b.relativePositionKm.y);
			EXPECT_DOUBLE_EQ(a.relativePositionKm.z, b.relativePositionKm.z);
			EXPECT_GT(std::abs(a.longitudeDeg - b.longitudeDeg), 10.0);
		}
	}

	TEST_F(IntervalMonteCarloTest, CancellationPreservesCompletedTrials)
	{
		std::stop_source stop;
		const auto result = hz::runIntervalMonteCarlo(source, config, stop.get_token(), [&](int done, int) {
			if (done == 2)
				stop.request_stop();
		});
		EXPECT_TRUE(result.cancelled);
		EXPECT_EQ(result.completed, 2);
		EXPECT_EQ(result.unstarted, 6);
		EXPECT_EQ(result.failedPropagation, 0);
		const auto before = hz::runIntervalMonteCarlo(source, config, stop.get_token());
		EXPECT_EQ(before.completed, 0);
		EXPECT_EQ(before.unstarted, 8);
	}

	TEST_F(IntervalMonteCarloTest, EarlierEntryDoesNotBecomeAnImpactInsideWindow)
	{
		config.startMjdTdb = 60000.001;
		const auto result = hz::runIntervalMonteCarlo(source, config);
		EXPECT_EQ(result.completed, 8);
		EXPECT_EQ(result.entriesBeforeInterval, 8);
		EXPECT_EQ(result.entryClones, 0);
		EXPECT_TRUE(result.events.empty());
	}

	TEST_F(IntervalMonteCarloTest, BadSourceAndCovarianceFailExplicitly)
	{
		source.fit.converged = false;
		EXPECT_THROW((void)hz::runIntervalMonteCarlo(source, config), std::invalid_argument);
		source.fit.converged = true;
		source.fit.uncertainty.fitParameterCovariance(0, 0) = -1;
		EXPECT_THROW((void)hz::runIntervalMonteCarlo(source, config), std::invalid_argument);
	}

	TEST_F(IntervalMonteCarloTest, ImportedSolutionUsesSameEngineAndRejectsInvalidCovariance)
	{
		orbital_core::OrbitSolution imported;
		imported.nominal = *source.fit.fitEpochOrbit;
		imported.covariance = source.fit.uncertainty.fitParameterCovariance;
		imported.nonGrav = source.fit.finalNonGrav;
		const auto fromFit = hz::runIntervalMonteCarlo(source, config);
		const auto fromFile = hz::runIntervalMonteCarlo(imported, source.effectiveConfig.nbody, config);
		ASSERT_EQ(fromFit.completed, fromFile.completed);
		ASSERT_EQ(fromFit.entryClones, fromFile.entryClones);
		for (int i = 0; i < config.samples; ++i) {
			EXPECT_EQ(fromFit.clones[i].state, fromFile.clones[i].state);
			ASSERT_EQ(fromFit.clones[i].entry.has_value(), fromFile.clones[i].entry.has_value());
			if (fromFit.clones[i].entry)
				EXPECT_DOUBLE_EQ(fromFit.clones[i].entry->mjdTdb, fromFile.clones[i].entry->mjdTdb);
		}
		imported.covariance(0, 0) = -1;
		EXPECT_THROW((void)hz::runIntervalMonteCarlo(imported, source.effectiveConfig.nbody, config),
		             std::invalid_argument);
		config.method = hz::McMethod::Observations;
		EXPECT_THROW((void)hz::runIntervalMonteCarlo(imported, source.effectiveConfig.nbody, config),
		             std::invalid_argument);
	}

	TEST_F(IntervalMonteCarloTest, CloseFlybyDoesNotCreateBoundaryEvent)
	{
		source.fit.fitEpochOrbit->setR({1, -.001, .0001});
		source.fit.fitEpochOrbit->setV({0, .02 / Constants::GAUSS_K, 0});
		config.endMjdTdb = 60000.1;
		const auto result = hz::runIntervalMonteCarlo(source, config);
		EXPECT_EQ(result.completed, 8);
		EXPECT_EQ(result.entryClones, 0);
		EXPECT_TRUE(result.events.empty());
		for (const auto& clone : result.clones)
			EXPECT_FALSE(clone.entry.has_value());
	}

	TEST_F(IntervalMonteCarloTest, InvalidTrajectoriesAreNotMisses)
	{
		source.fit.fitEpochOrbit->setR({1, 0, 0});
		const auto result = hz::runIntervalMonteCarlo(source, config);
		EXPECT_EQ(result.completed, 0);
		EXPECT_EQ(result.failedPropagation, 8);
		EXPECT_EQ(result.nominal.status, hz::McCloneStatus::PropagationFailed);
		EXPECT_TRUE(result.events.empty());
	}

	TEST_F(IntervalMonteCarloTest, SamplesFullNonGravCovarianceWithoutMutatingSource)
	{
		source.fit.finalNonGrav.setModel(NgModel::Yarkovsky);
		source.fit.finalNonGrav.setEnabled(true);
		source.fit.finalNonGrav.setParams(0, -3e-13, 0);
		source.fit.uncertainty.fitParameterCovariance = Eigen::MatrixXd::Identity(7, 7) * 1e-22;
		source.fit.uncertainty.fitParameterCovariance(6, 6) = 1e-32;
		const auto result = hz::runIntervalMonteCarlo(source, config);
		ASSERT_EQ(result.completed, 8);
		EXPECT_DOUBLE_EQ(source.fit.finalNonGrav.A2(), -3e-13);
		for (const auto& clone : result.clones) {
			EXPECT_EQ(clone.nonGrav.fitParameterCount(), 1);
			EXPECT_NE(clone.nonGrav.A2(), source.fit.finalNonGrav.A2());
			EXPECT_DOUBLE_EQ(clone.nonGrav.A1(), 0);
			EXPECT_DOUBLE_EQ(clone.nonGrav.A3(), 0);
		}
	}

	TEST_F(IntervalMonteCarloTest, ProgressCallbackFailureStopsWithoutTerminatingWorker)
	{
		EXPECT_THROW(
		    (void)hz::runIntervalMonteCarlo(source, config, {}, [](int, int) { throw std::runtime_error("callback"); }),
		    std::runtime_error);
	}

	TEST(IntervalMcEvents, CrossingsOnlyAndCompleteTrialDenominator)
	{
		hz::IntervalMcResult result;
		result.clones.resize(4);
		for (int i = 0; i < 4; ++i)
			result.clones[i].index = i;
		result.clones[0].status = hz::McCloneStatus::Complete;
		result.clones[1].status = hz::McCloneStatus::Complete;
		result.clones[2].status = hz::McCloneStatus::PropagationFailed;
		result.clones[3].status = hz::McCloneStatus::Cancelled;
		result.clones[0].entry = hz::BoundaryCrossing{.mjdTdb = 60000.5};
		result.clones[2].entry = hz::BoundaryCrossing{.mjdTdb = 60001.5};
		result.nominal.status = hz::McCloneStatus::Complete;
		result.nominal.entry = hz::BoundaryCrossing{.mjdTdb = 60002.5};
		hz::mc::summarize(result);
		ASSERT_EQ(result.events.size(), 1);
		EXPECT_EQ(result.completed, 2);
		EXPECT_EQ(result.failedPropagation, 1);
		EXPECT_DOUBLE_EQ(result.events[0].probability, .5);
		EXPECT_EQ(result.events[0].entryClones, 1);
		EXPECT_GT(result.events[0].probabilityHigh, .5);
	}

	TEST(IntervalMcEvents, MidnightDoesNotSplitOneEpisode)
	{
		hz::IntervalMcResult result;
		result.clones.resize(2);
		for (int i = 0; i < 2; ++i) {
			auto& clone = result.clones[i];
			clone.index = i;
			clone.status = hz::McCloneStatus::Complete;
			clone.entry = hz::BoundaryCrossing{
			    .mjdTdb = orbital_core::time::toMjdTdb(60001 + (i ? .0001 : -.0001), orbital_core::time::Scale::UTC)};
		}
		hz::mc::summarize(result);
		ASSERT_EQ(result.events.size(), 1);
		EXPECT_EQ(result.events[0].entryClones, 2);
		EXPECT_DOUBLE_EQ(result.events[0].probability, 1);
		EXPECT_DOUBLE_EQ(result.events[0].startMjdTdb, result.clones[0].entry->mjdTdb);
		EXPECT_DOUBLE_EQ(result.events[0].endMjdTdb, result.clones[1].entry->mjdTdb);
	}

	TEST(IntervalMcEvents, SortedConsecutiveGapsSeparateEpisodes)
	{
		hz::IntervalMcResult result;
		for (const double epoch : {60003.0, 60000.0, 60001.0, 60000.5}) {
			hz::McClone clone;
			clone.status = hz::McCloneStatus::Complete;
			clone.entry = hz::BoundaryCrossing{.mjdTdb = epoch};
			result.clones.push_back(clone);
		}
		hz::mc::summarize(result);
		ASSERT_EQ(result.events.size(), 2);
		EXPECT_EQ(result.events[0].entryClones, 3);
		EXPECT_DOUBLE_EQ(result.events[0].startMjdTdb, 60000);
		EXPECT_DOUBLE_EQ(result.events[0].endMjdTdb, 60001);
		EXPECT_DOUBLE_EQ(result.events[0].probability, .75);
		EXPECT_EQ(result.events[1].entryClones, 1);
		EXPECT_DOUBLE_EQ(result.events[1].probability, .25);
	}

	TEST(IntervalMcReport, ExistingArchiveRendersWithoutPropagationOrNumericalChanges)
	{
		using Json = nlohmann::json;
		const auto dir =
		    std::filesystem::temp_directory_path() /
		    ("interval-report-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::create_directory(dir);

		struct Cleanup {
			std::filesystem::path path;

			~Cleanup() { std::filesystem::remove_all(path); }
		} cleanup{dir};

		const double tt = orbital_core::time::toMjdTt(54746.115, orbital_core::time::Scale::UTC);
		Json visit = {{"start", tt - .1},
		              {"end", tt},
		              {"closest", tt},
		              {"entry",
		               {{"epoch", tt},
		                {"latitude_deg", 21.09},
		                {"longitude_deg", 30.53},
		                {"relative_velocity_km_s", {1, 2, 3}}}}};
		Json data = {{"schema", "orbital-interval-mc-1"},
		             {"time_scale", "MJD TT"},
		             {"target", 2},
		             {"start", tt - .1},
		             {"end", tt + .1},
		             {"events", {{{"start", tt - .1}, {"end", tt}, {"closest", tt}}}},
		             {"nominal", {{"index", -1}, {"status", 1}, {"encounters", {visit}}}},
		             {"clones",
		              {{{"index", 0}, {"status", 1}, {"encounters", {visit}}},
		               {{"index", 1}, {"status", 3}, {"encounters", {visit}}}}},
		             {"untrusted_text", "</script><script>alert(1)</script>"}};
		const auto original = data.dump();
		{
			std::ofstream file(dir / "result.json");
			file << original;
		}
		ASSERT_NO_THROW(hz::renderIntervalMonteCarloReport(dir));
		std::ifstream file(dir / "report.html");
		const std::string html{std::istreambuf_iterator<char>(file), {}};
		const std::string marker = "<script id=\"mc-data\" type=\"application/json\">";
		const auto start = html.find(marker);
		ASSERT_NE(start, std::string::npos);
		const auto end = html.find("</script>", start + marker.size());
		const auto payload = Json::parse(html.substr(start + marker.size(), end - start - marker.size()));
		EXPECT_EQ(payload["target_name"], "Earth");
		ASSERT_EQ(payload["entries"].size(), 2);
		EXPECT_EQ(payload["entries"][0]["clone"], -1);
		EXPECT_EQ(payload["entries"][1]["clone"], 0);
		EXPECT_EQ(payload["entries"][1]["event"], 0);
		EXPECT_NEAR(payload["entries"][1]["epoch_utc"].get<double>(), 54746.115, 1e-10);
		EXPECT_DOUBLE_EQ(payload["entries"][1]["latitude_deg"], 21.09);
		EXPECT_EQ(html.find("</script><script>alert(1)"), std::string::npos);
		EXPECT_EQ(html.find("__LEAFLET"), std::string::npos);
		EXPECT_NE(html.find("Satellite imagery"), std::string::npos);
		std::ifstream numerical(dir / "result.json");
		EXPECT_EQ(std::string(std::istreambuf_iterator<char>(numerical), {}), original);
		data["nominal"]["status"] = 3;
		data["clones"] = Json::array();
		{
			std::ofstream output(dir / "result.json");
			output << data;
		}
		EXPECT_NO_THROW(hz::renderIntervalMonteCarloReport(dir));
		data["schema"] = "orbital-interval-mc-2";
		data["nominal"]["entry"] = nullptr;
		data["nominal"].erase("encounters");
		data["clones"] = {{{"index", 0}, {"status", 1}, {"entry", visit.at("entry")}}};
		{
			std::ofstream output(dir / "result.json");
			output << data;
		}
		ASSERT_NO_THROW(hz::renderIntervalMonteCarloReport(dir));
		std::ifstream revised(dir / "report.html");
		const std::string newHtml{std::istreambuf_iterator<char>(revised), {}};
		EXPECT_NE(newHtml.find("\"entries\":1"), std::string::npos);
		EXPECT_EQ(newHtml.find("reporting_radius_au"), std::string::npos);
		EXPECT_NE(newHtml.find("#ff5b57"), std::string::npos);
		data["schema"] = "orbital-interval-mc-4";
		data["time_scale"] = "MJD TDB";
		const auto toTdb = [](double value) {
			return orbital_core::time::toMjdTdb(value, orbital_core::time::Scale::TT);
		};
		data["start"] = toTdb(data["start"].get<double>());
		data["end"] = toTdb(data["end"].get<double>());
		data["clones"][0]["entry"]["epoch"] = toTdb(tt);
		{
			std::ofstream output(dir / "result.json");
			output << data;
		}
		ASSERT_NO_THROW(hz::renderIntervalMonteCarloReport(dir));
		std::ifstream tdbFile(dir / "report.html");
		const std::string tdbHtml{std::istreambuf_iterator<char>(tdbFile), {}};
		const auto tdbStart = tdbHtml.find(marker) + marker.size();
		const auto tdbEnd = tdbHtml.find("</script>", tdbStart);
		const auto tdbPayload = Json::parse(tdbHtml.substr(tdbStart, tdbEnd - tdbStart));
		EXPECT_EQ(tdbPayload["time_scale"], "MJD TDB");
		EXPECT_DOUBLE_EQ(tdbPayload["entries"][0]["epoch"], toTdb(tt));
		EXPECT_DOUBLE_EQ(tdbPayload["entries"][0]["epoch_utc"], payload["entries"][1]["epoch_utc"]);
		data["time_scale"] = "MJD TT";
		{
			std::ofstream output(dir / "result.json");
			output << data;
		}
		EXPECT_THROW(hz::renderIntervalMonteCarloReport(dir), std::runtime_error);
	}
}  // namespace
