#include "strewn_core/StrewnJob.h"
#include "strewn_core/AtmosphericInput.h"
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <gdal_priv.h>
#include <numeric>

namespace {
	using namespace strewn;
	using Json = nlohmann::json;

	Json jsonScript(const std::string& html, const std::string& id)
	{
		const std::string marker = "<script id=\"" + id + "\" type=\"application/json\">";
		const auto position = html.find(marker);
		if (position == std::string::npos)
			throw std::runtime_error("Missing report payload: " + id);
		const auto first = position + marker.size();
		return Json::parse(html.substr(first, html.find("</script>", first) - first));
	}

	Json reportData(const std::string& html)
	{
		auto data = jsonScript(html, "field-data");
		for (auto& run : data["runs"])
			run["flight"] = jsonScript(html, run["flight"]["element"]);
		return data;
	}

	class JobTest : public testing::Test {
	   protected:
		std::filesystem::path m_directory;

		void SetUp() override
		{
			static std::atomic<unsigned> counter{0};
			m_directory =
			    std::filesystem::temp_directory_path() /
			    ("strewn-test-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "-" +
			     std::to_string(counter++));
			std::filesystem::create_directory(m_directory);
		}

		void TearDown() override
		{
			const auto parent = std::filesystem::canonical(std::filesystem::temp_directory_path());
			if (std::filesystem::canonical(m_directory).parent_path() == parent &&
			    m_directory.filename().string().starts_with("strewn-test-"))
				std::filesystem::remove_all(m_directory);
		}

		Json config() const
		{
			return {{"schema_version", 2},
			        {"input",
			         {{"kind", "fixed_state"},
			          {"state",
			           {{"epoch_mjd_utc", 59988.125},
			            {"latitude_deg", 49.8},
			            {"longitude_deg", 0.8},
			            {"height_wgs84_m", 100.0},
			            {"velocity_enu_mps", {0, 0, -20}}}}}},
			        {"environment", {{"mode", "idealized"}}},
			        {"physics", {{"phase", "dark_flight"}}},
			        {"monte_carlo", {{"runs", 3}, {"workers", 2}}},
			        {"output_directory", "result"}};
		}

		JobSpec load(const Json& json)
		{
			const auto path = m_directory / "input.json";
			{
				std::ofstream out(path);
				out << json.dump(2);
			}
			return loadJob(path);
		}
	};

	TEST_F(JobTest, WritesReplayableResults)
	{
		const auto spec = load(config());
		const auto summary = StrewnJob{}.run(spec).simulation;
		EXPECT_EQ(summary.completed, 3U);
		EXPECT_EQ(summary.failed, 0U);
		EXPECT_TRUE(std::filesystem::exists(spec.outputDirectory / "impacts.csv"));
		const auto replay = loadJob(spec.outputDirectory / "job.json");
		const auto repeated = StrewnJob{}.run(replay).simulation;
		EXPECT_EQ(repeated.trackedImpacts, summary.trackedImpacts);
		const auto read = [](const auto& path) {
			std::ifstream in(path);
			return std::string(std::istreambuf_iterator<char>(in), {});
		};
		EXPECT_EQ(read(spec.outputDirectory / "impacts.csv"), read(replay.outputDirectory / "impacts.csv"));
		EXPECT_EQ(read(spec.outputDirectory / "runs.jsonl"), read(replay.outputDirectory / "runs.jsonl"));
		EXPECT_THROW(StrewnJob{}.run(spec), std::invalid_argument);
	}

	TEST_F(JobTest, FragmentRemainderDefaultsAreVersionedWithoutChangingPressure)
	{
		for (int schema : {1, 2}) {
			auto c = config();
			c["schema_version"] = schema;
			for (bool omitPhysics : {false, true}) {
				if (omitPhysics)
					c.erase("physics");
				auto job = load(c);
				EXPECT_EQ(job.physics.ramPressure, RamPressure::RhoV2);
				EXPECT_EQ(job.physics.retainMassRemainder, schema == 2);
				const auto file = m_directory / "saved.json";
				saveJobConfiguration(file, job);
				const auto replay = loadJob(file);
				EXPECT_EQ(replay.physics.ramPressure, job.physics.ramPressure);
				EXPECT_EQ(replay.physics.retainMassRemainder, job.physics.retainMassRemainder);
				std::ifstream in(file);
				const auto saved = Json::parse(in);
				EXPECT_EQ(saved.at("schema_version"), 2);
				EXPECT_EQ(saved.at("physics").at("retain_mass_remainder"), schema == 2);
			}
			c["physics"] = {{"ram_pressure", "half_rho_v2"}, {"retain_mass_remainder", true}};
			const auto explicitJob = load(c);
			EXPECT_EQ(explicitJob.physics.ramPressure, RamPressure::HalfRhoV2);
			EXPECT_TRUE(explicitJob.physics.retainMassRemainder);
		}
		auto invalid = config();
		invalid["schema_version"] = 3;
		EXPECT_THROW(load(invalid), std::invalid_argument);
	}

	TEST_F(JobTest, BothFragmentRemainderPoliciesReplayFullCascades)
	{
		const auto read = [](const auto& path) {
			std::ifstream in(path);
			return std::string(std::istreambuf_iterator<char>(in), {});
		};
		for (int schema : {1, 2}) {
			auto c = config();
			c["schema_version"] = schema;
			c["physics"] = {{"phase", "entry"}};
			c["input"]["state"]["height_wgs84_m"] = 30000;
			c["input"]["state"]["velocity_enu_mps"] = {9000, 0, -9000};
			c["output_directory"] = "cascade-" + std::to_string(schema);
			c["monte_carlo"]["parameters"] = {{"size", {{"mass_kg", 500}}}, {"strong_strength_pa", 2e6}};
			const auto job = load(c);
			const auto result = StrewnJob{}.run(job).simulation;
			ASSERT_EQ(result.failed, 0U);
			ASSERT_EQ(result.truncated, 0U);
			ASSERT_GT(result.trackedImpacts, 0U);
			const auto replay = loadJob(job.outputDirectory / "job.json");
			EXPECT_EQ(replay.physics.retainMassRemainder, schema == 2);
			ASSERT_EQ(StrewnJob{}.run(replay).simulation.failed, 0U);
			for (const auto* file : {"impacts.csv", "breakups.jsonl", "runs.jsonl"}) {
				ASSERT_TRUE(std::filesystem::exists(job.outputDirectory / file));
				ASSERT_TRUE(std::filesystem::exists(replay.outputDirectory / file));
				EXPECT_EQ(read(job.outputDirectory / file), read(replay.outputDirectory / file)) << file;
			}
		}
	}

	TEST_F(JobTest, GeneratedConfigurationRoundTripsAllEffectiveSettings)
	{
		auto job = load(config());
		job.environment = AutomaticEnvironment{};
		job.environmentDirectory = m_directory / "shared-cache";
		job.physics.representativesPerBin = 3;
		job.monteCarlo.runs = 1000;
		job.monteCarlo.workers = 15;
		job.monteCarlo.runTimeoutS = 123.5;
		job.monteCarlo.parameters.densityKgM3 = {2700, 3400};
		job.monteCarlo.parameters.size = PhotometricSize{30.311, .15, .15, .15};
		job.monteCarlo.convergence = ConvergenceConfig{};
		job.monteCarlo.convergence->centerDeg = {49.8, .8};
		job.maps = MapOptions{};
		const auto file = m_directory / "generated.json";
		saveJobConfiguration(file, job);
		const auto loaded = loadJob(file);
		EXPECT_EQ(loaded.environmentDirectory, job.environmentDirectory);
		EXPECT_EQ(loaded.outputDirectory, job.outputDirectory);
		EXPECT_EQ(loaded.physics.representativesPerBin, 3);
		EXPECT_EQ(loaded.monteCarlo.runs, 1000);
		EXPECT_EQ(loaded.monteCarlo.workers, 15);
		EXPECT_EQ(loaded.monteCarlo.runTimeoutS, 123.5);
		EXPECT_EQ(loaded.physics.maxParticles, 0U);
		EXPECT_EQ(loaded.monteCarlo.parameters.densityKgM3.min, 2700);
		EXPECT_EQ(std::get<PhotometricSize>(loaded.monteCarlo.parameters.size).absoluteMagnitude, 30.311);
		EXPECT_EQ(loaded.monteCarlo.convergence->threshold, .02);
		EXPECT_TRUE(loaded.maps->autoCellSize);
		job.monteCarlo.parameters.densityKgM3 = {4000, 2500};
		EXPECT_THROW(saveJobConfiguration(file, job), std::invalid_argument);
		EXPECT_EQ(loadJob(file).monteCarlo.parameters.densityKgM3.min, 2700);
	}

	TEST_F(JobTest, MaterialPresetsMatchTableAndPreserveRunControls)
	{
		const std::array<Json, 3> expected{{{{"density_kg_m3", {2500, 3500}},
		                                     {"strong_strength_pa", {.9e6, 5e6}},
		                                     {"dust_fraction", {0, .9}},
		                                     {"largest_fraction", {.05, .4}},
		                                     {"beta", {.4, .8}}},
		                                    {{"density_kg_m3", {2100, 3100}},
		                                     {"strong_strength_pa", {.25e6, .7e6}},
		                                     {"dust_fraction", {.3, .95}},
		                                     {"largest_fraction", {.01, .25}},
		                                     {"beta", {.5, .95}}},
		                                    {{"density_kg_m3", {7000, 7800}},
		                                     {"strong_strength_pa", {5e6, 50e6}},
		                                     {"dust_fraction", {0, .25}},
		                                     {"largest_fraction", {.4, .9}},
		                                     {"beta", {.1, .4}}}}};
		const std::array<Json, 3> albedo{{{.075, .208, .33}, {.02, .06, .18}, {.07, .18, .54}}};
		auto job = load(config());
		job.monteCarlo.parameters.size = PhotometricSize{30.311, .15, .15, .15};
		job.monteCarlo.seed = 12345;
		job.monteCarlo.workers = 15;
		job.physics.weibullAlpha = .2;
		job.maps = MapOptions{};
		const auto file = m_directory / "preset.json";
		saveJobConfiguration(file, job);
		const auto read = [&] {
			std::ifstream in(file);
			return Json::parse(in);
		};
		const auto original = read();
		for (int i = 0; i < 3; ++i) {
			const auto preset = static_cast<MaterialPreset>(i);
			applyMaterialPreset(job, preset);
			EXPECT_EQ(materialPreset(job), preset);
			saveJobConfiguration(file, job);
			const auto saved = read();
			const auto& parameters = saved.at("monte_carlo").at("parameters");
			for (const auto& [key, value] : expected[i].items())
				EXPECT_EQ(parameters.at(key), value) << key;
			EXPECT_EQ(parameters.at("size").at("albedo"), albedo[i]);
			EXPECT_EQ(parameters.at("size").at("absolute_magnitude"), 30.311);
			EXPECT_EQ(parameters.at("weak_probability"), 0);
			EXPECT_EQ(parameters.at("weak_strength_pa"), Json({.04e6, .12e6}));
			EXPECT_EQ(parameters.at("wind_scale"), Json({.85, 1.15}));
			EXPECT_EQ(parameters.at("drag_scale"), Json({.85, 1.15}));
			EXPECT_EQ(parameters.at("ablation_kg_j"), Json({7e-8, 9e-8}));
			for (const auto* key : {"input", "environment", "physics", "integration", "maps", "output_directory"})
				EXPECT_EQ(saved.at(key), original.at(key)) << key;
			for (const auto* key : {"seed", "runs", "workers"})
				EXPECT_EQ(saved.at("monte_carlo").at(key), original.at("monte_carlo").at(key));
			EXPECT_EQ(materialPreset(loadJob(file)), preset);
		}
		job.monteCarlo.parameters.densityKgM3.min = 7100;
		EXPECT_EQ(materialPreset(job), MaterialPreset::Custom);
		EXPECT_THROW(applyMaterialPreset(job, MaterialPreset::Custom), std::invalid_argument);
	}

	TEST_F(JobTest, MaterialPresetRequiresPhotometricSizeBeforeChangingParameters)
	{
		auto job = load(config());
		job.monteCarlo.parameters.size = MassSize{{5, 5}};
		const auto density = job.monteCarlo.parameters.densityKgM3;
		EXPECT_EQ(materialPreset(job), MaterialPreset::Custom);
		EXPECT_THROW(applyMaterialPreset(job, MaterialPreset::Iron), std::invalid_argument);
		EXPECT_EQ(std::get<MassSize>(job.monteCarlo.parameters.size).kg.min, 5);
		EXPECT_EQ(job.monteCarlo.parameters.densityKgM3.min, density.min);
	}

	TEST_F(JobTest, AutomaticOfflineFailureAndExplicitFallback)
	{
		auto c = config();
		c["environment"] = {{"mode", "automatic"}, {"offline", true}};
		auto spec = load(c);
		EXPECT_THROW(StrewnJob{}.run(spec), std::runtime_error);
		EXPECT_FALSE(std::filesystem::exists(spec.outputDirectory / "impacts.csv"));
		c["output_directory"] = "fallback";
		c["environment"]["on_unavailable"] = "idealized";
		spec = load(c);
		EXPECT_EQ(StrewnJob{}.run(spec).simulation.failed, 0U);
		std::ifstream in(spec.outputDirectory / "manifest.json");
		EXPECT_TRUE(Json::parse(in).at("idealized_environment").get<bool>());
		std::ifstream environment(spec.outputDirectory / "environment.json");
		EXPECT_TRUE(Json::parse(environment).at("idealized_fallback").get<bool>());
	}

	TEST_F(JobTest, ConvergenceAndRepresentativesRoundTrip)
	{
		auto c = config();
		c["physics"]["representatives_per_bin"] = 3;
		c["monte_carlo"]["runs"] = 50;
		c["monte_carlo"]["convergence"] = {{"center_deg", {49.8, .8}},
		                                   {"half_size_km", {1, 1}},
		                                   {"batch_size", 2},
		                                   {"consecutive_checks", 3},
		                                   {"threshold", .02}};
		const auto spec = load(c);
		const auto s = StrewnJob{}.run(spec).simulation;
		EXPECT_TRUE(s.converged);
		EXPECT_EQ(s.completed, 8U);
		const auto replay = loadJob(spec.outputDirectory / "job.json");
		EXPECT_EQ(replay.physics.representativesPerBin, 3U);
		ASSERT_TRUE(replay.monteCarlo.convergence);
		EXPECT_EQ(replay.monteCarlo.convergence->batchSize, 2U);
		EXPECT_DOUBLE_EQ(replay.monteCarlo.convergence->threshold, .02);
		const auto repeated = StrewnJob{}.run(replay).simulation;
		EXPECT_EQ(repeated.completed, s.completed);
		std::ifstream stream(spec.outputDirectory / "summary.json");
		const auto summary = Json::parse(stream);
		EXPECT_EQ(summary["stop_reason"], "js_converged");
		EXPECT_TRUE(summary["convergence_evaluated"].get<bool>());
		EXPECT_EQ(summary["convergence_history"].size(), 4U);
	}

	TEST_F(JobTest, RejectsInvalidRepresentativesAndConvergence)
	{
		for (Json bad : {Json(0), Json(-1), Json(1.5), Json(1025)}) {
			auto c = config();
			c["physics"]["representatives_per_bin"] = bad;
			EXPECT_THROW(load(c), std::invalid_argument);
		}
		for (const auto& [name, bad] : std::vector<std::pair<std::string, Json>>{
		         {"threshold", 0}, {"batch_size", 1.5}, {"cell_m", 0}, {"half_size_km", {1e10, 1e10}}}) {
			auto c = config();
			c["monte_carlo"]["convergence"] = {{"center_deg", {49.8, .8}}, {name, bad}};
			EXPECT_THROW(load(c), std::invalid_argument);
		}
	}

	TEST_F(JobTest, AutomaticRejectsBadBoundsAndAmbiguousSourcesBeforeNetwork)
	{
		auto c = config();
		c["environment"] = {{"mode", "automatic"}, {"terrain_step_deg", 0}};
		EXPECT_THROW(load(c), std::invalid_argument);
		c["environment"] = {{"mode", "automatic"}, {"bounds_deg", {49, 50, 170, -170}}};
		EXPECT_THROW(load(c), std::invalid_argument);
		c["environment"] = {{"mode", "automatic"}, {"cycle_mjd_utc", 59988}, {"gfs_file", "weather.grb2"}};
		EXPECT_THROW(load(c), std::invalid_argument);
		c["environment"] = {
		    {"mode", "automatic"}, {"offline", true}, {"on_unavailable", "idealized"}, {"bounds_deg", {40, 41, 0, 1}}};
		const auto spec = load(c);
		EXPECT_THROW(StrewnJob{}.run(spec), std::invalid_argument);
	}

	TEST_F(JobTest, RejectsUnknownOrAmbiguousInputs)
	{
		auto c = config();
		c["physics"]["strenght_pa"] = 1;
		EXPECT_THROW(load(c), std::invalid_argument);
		c = config();
		c["input"]["state"]["position_ecef_m"] = {1, 2, 3};
		EXPECT_THROW(load(c), std::invalid_argument);
		c = config();
		c["monte_carlo"]["parameters"]["size"] = {{"mass_kg", 1}, {"diameter_m", 1}};
		EXPECT_THROW(load(c), std::invalid_argument);
		c = config();
		c["input"]["kind"] = "orbit_solution";
		EXPECT_THROW(load(c), std::invalid_argument);
		c = config();
		c["input"]["state"]["epoch_jd_utc"] = 2459988.625;
		EXPECT_THROW(load(c), std::invalid_argument);
		c = config();
		c["input"]["state"].erase("epoch_mjd_utc");
		EXPECT_THROW(load(c), std::invalid_argument);
	}

	TEST_F(JobTest, MjdPrecisionSurvivesArchiveAndReplay)
	{
		const double epoch = std::nextafter(59988.125, 60000.0);
		ASSERT_NE(epoch, (epoch + 2400000.5) - 2400000.5);
		auto c = config();
		c["input"]["state"]["epoch_mjd_utc"] = epoch;
		const auto spec = load(c);
		EXPECT_DOUBLE_EQ(std::get<EntryState>(spec.input).epochMjdUtc, epoch);
		ASSERT_EQ(StrewnJob{}.run(spec).simulation.failed, 0U);
		std::ifstream in(spec.outputDirectory / "entries.json");
		const auto states = Json::parse(in);
		EXPECT_EQ(states.at("entries").at(0).at("epoch_mjd_utc").get<double>(), epoch);
		EXPECT_FALSE(states.at("entries").at(0).contains("epoch_jd_utc"));
		const auto replay = loadJob(spec.outputDirectory / "job.json");
		ASSERT_EQ(StrewnJob{}.run(replay).simulation.failed, 0U);
		std::ifstream repeated(replay.outputDirectory / "entries.json");
		EXPECT_EQ(Json::parse(repeated), states);
	}

	TEST_F(JobTest, LegacyJdConvertsOnlyOnImport)
	{
		auto c = config();
		c["input"]["state"].erase("epoch_mjd_utc");
		c["input"]["state"]["epoch_jd_utc"] = 2459988.625;
		EXPECT_DOUBLE_EQ(std::get<EntryState>(load(c).input).epochMjdUtc, 59988.125);
	}

	TEST_F(JobTest, PreparedEntriesUseTheSameArchiveAndReplay)
	{
		const auto original = load(config());
		const auto input = m_directory / "input.json";
		const auto spec = loadJob(input, EntrySamples{{std::get<EntryState>(original.input)}});
		ASSERT_EQ(StrewnJob{}.run(spec).simulation.failed, 0U);
		std::ifstream in(spec.outputDirectory / "entries.json");
		const auto entries = Json::parse(in);
		EXPECT_EQ(entries.at("entries").size(), 1);
		EXPECT_EQ(entries.at("entries").at(0).at("epoch_mjd_utc"), 59988.125);
		const auto replay = loadJob(spec.outputDirectory / "job.json");
		EXPECT_EQ(StrewnJob{}.run(replay).simulation.failed, 0U);
	}

	TEST_F(JobTest, RejectsMissingEnvironmentInsteadOfSilentFallback)
	{
		auto c = config();
		c["environment"] = {{"mode", "file"}, {"path", "missing.json"}};
		EXPECT_THROW(StrewnJob{}.run(load(c)), std::runtime_error);
		EXPECT_FALSE(std::filesystem::exists(m_directory / "result"));
	}

	TEST_F(JobTest, RemoteEnvironmentConfigurationIsExplicitAndBounded)
	{
		auto c = config();
		c["environment"] = {{"mode", "remote"},
		                    {"weather_url", "https://example.org/weather.json"},
		                    {"terrain_url", "https://example.org/terrain.json"},
		                    {"timeout_seconds", 20},
		                    {"max_download_bytes", 4096}};
		const auto spec = load(c);
		const auto& remote = std::get<EnvironmentUrls>(spec.environment);
		EXPECT_EQ(remote.timeoutSeconds, 20U);
		EXPECT_EQ(remote.maxDownloadBytes, 4096U);
		auto invalid = c;
		invalid["environment"]["path"] = "ambiguous.json";
		EXPECT_THROW(load(invalid), std::invalid_argument);
		for (const auto& [name, value] :
		     std::vector<std::pair<std::string, Json>>{{"timeout_seconds", 0},
		                                               {"timeout_seconds", 3601},
		                                               {"timeout_seconds", true},
		                                               {"max_download_bytes", -1},
		                                               {"max_download_bytes", 1073741825}}) {
			invalid = c;
			invalid["environment"][name] = value;
			EXPECT_THROW(load(invalid), std::invalid_argument);
		}
		std::stop_source stop;
		stop.request_stop();
		EXPECT_TRUE(StrewnJob{}.run(spec, {stop.get_token(), {}}).simulation.cancelled);
		EXPECT_FALSE(std::filesystem::exists(spec.outputDirectory));
	}

	TEST_F(JobTest, OrbitAdapterCanChooseOutputDirectoryButFixedInputCannot)
	{
		auto c = config();
		c.erase("output_directory");
		EXPECT_THROW(load(c), std::invalid_argument);
		c["input"] = {{"kind", "orbit_solution"}, {"path", "work/object"}};
		const auto path = m_directory / "orbit.json";
		{
			std::ofstream out(path);
			out << c.dump();
		}
		const auto spec = loadJob(path, EntrySamples{});
		EXPECT_TRUE(spec.outputDirectory.empty());
	}

	TEST_F(JobTest, RejectsInvalidGroupingCounts)
	{
		for (const Json& value : {Json::array({-1, 10}), Json::array({30, 10, 5}), Json::array({3.5, 10})}) {
			auto c = config();
			c["physics"]["individual_fragments"] = value;
			EXPECT_THROW(load(c), std::invalid_argument);
		}
		auto c = config();
		c["physics"]["bins_per_decade"] = {4294967297ULL, 2};
		EXPECT_THROW(load(c), std::invalid_argument);
	}

	TEST_F(JobTest, RejectsExtraVectorAndRangeComponents)
	{
		auto c = config();
		c["input"]["state"]["velocity_enu_mps"] = {1, 2, 3, 4};
		EXPECT_THROW(load(c), std::invalid_argument);
		c = config();
		c["monte_carlo"]["parameters"]["density_kg_m3"] = {2500, 3000, 3500};
		EXPECT_THROW(load(c), std::invalid_argument);
	}

	TEST_F(JobTest, FullCascadeIsByteIdenticalWithOneAndFifteenWorkers)
	{
		auto c = config();
		c["input"]["state"]["height_wgs84_m"] = 30000;
		c["input"]["state"]["velocity_enu_mps"] = {9000, 0, -9000};
		c["physics"]["phase"] = "entry";
		c["monte_carlo"] = {{"runs", 24}, {"workers", 1}, {"parameters", {{"size", {{"mass_kg", 100}}}}}};
		const auto serial = load(c);
		const auto first = StrewnJob{}.run(serial).simulation;
		ASSERT_EQ(first.failed, 0U);
		ASSERT_EQ(first.truncated, 0U);
		ASSERT_GT(first.trackedImpacts, first.completed);
		const auto read = [](const auto& path) {
			std::ifstream in(path);
			return std::string(std::istreambuf_iterator<char>(in), {});
		};
		ASSERT_FALSE(read(serial.outputDirectory / "breakups.jsonl").empty());
		auto parallel = serial;
		parallel.outputDirectory = m_directory / "parallel";
		parallel.monteCarlo.workers = 15;
		const auto second = StrewnJob{}.run(parallel).simulation;
		EXPECT_EQ(second.workers, 15U);
		EXPECT_EQ(second.completed, first.completed);
		for (const auto* file : {"impacts.csv", "runs.jsonl", "breakups.jsonl"})
			EXPECT_EQ(read(serial.outputDirectory / file), read(parallel.outputDirectory / file)) << file;
	}

	TEST_F(JobTest, CancelledOutputIsLabelled)
	{
		std::stop_source cancel;
		RunControl control{cancel.get_token(), [&](const auto&) { cancel.request_stop(); }};
		const auto spec = load(config());
		const auto summary = StrewnJob{}.run(spec, control).simulation;
		EXPECT_TRUE(summary.cancelled);
		std::ifstream in(spec.outputDirectory / "manifest.json");
		EXPECT_EQ(Json::parse(in)["status"], "cancelled");
	}

	TEST_F(JobTest, MapsAreLibraryOwnedAndReplayable)
	{
		auto c = config();
		c["maps"] = {{"method", "fixed"},      {"title", "CX1 <test>"}, {"center_deg", {49.8, .8}},
		             {"half_size_km", {1, 1}}, {"image_width", 800},    {"basemap", "esri_imagery"},
		             {"offline", true}};
		const auto spec = load(c);
		const auto result = StrewnJob{}.run(spec);
		ASSERT_EQ(result.artifacts.status, ArtifactStatus::Complete) << result.artifacts.error;
		EXPECT_FALSE(result.artifacts.warnings.empty());
		EXPECT_EQ(result.simulation.completed, 3U);
		EXPECT_EQ(result.directory, spec.outputDirectory);
		const auto maps = spec.outputDirectory / "maps";
		for (const auto* name : {"density.tif", "density.json", "contours.geojson", "field.html", "artifacts.json"})
			EXPECT_GT(std::filesystem::file_size(maps / name), 50U) << name;
		for (const auto& file : std::filesystem::directory_iterator(maps))
			EXPECT_NE(file.path().extension(), ".png");
		std::ifstream html(maps / "field.html", std::ios::binary);
		const std::string content(std::istreambuf_iterator<char>(html), {});
		EXPECT_NE(content.find("CX1 &lt;test&gt;"), std::string::npos);
		EXPECT_EQ(content.find("<img src=\"data:image/png;base64,"), std::string::npos);
		EXPECT_EQ(content.find("<script src="), std::string::npos);
		EXPECT_NE(content.find("L.map('map'"), std::string::npos);
		EXPECT_NE(content.find("L.geoJSON(data.contours"), std::string::npos);
		EXPECT_NE(content.find("L.circleMarker"), std::string::npos);
		EXPECT_EQ(content.find("Auto: 1x / 40x"), std::string::npos);
		EXPECT_EQ(content.find("time<model.firstBreakup+5"), std::string::npos);
		EXPECT_NE(content.find("<option value=\"1\">1x</option>"), std::string::npos);
		const auto view = reportData(content);
		EXPECT_EQ(view["entry_summary"]["positive_weight_entries"], 1);
		EXPECT_EQ(view["entry_summary"]["time_sigma_s"], nullptr);
		EXPECT_NEAR(view["entry_summary"]["latitude_deg"]["mean"].get<double>(), 49.8, 1e-9);
		EXPECT_NEAR(view["entry_summary"]["flight_path_deg"]["mean"].get<double>(), -90., 1e-9);
		EXPECT_EQ(view["entry_summary"]["heading_deg"]["mean"], nullptr);
		ASSERT_EQ(view["trace_columns"].size(), 8U);
		ASSERT_EQ(view["runs"].size(), 3U);
		for (const auto& run : view["runs"]) {
			const auto& flight = run["flight"];
			ASSERT_TRUE(flight["verified"].get<bool>()) << flight.dump();
			EXPECT_LE(flight["max_impact_delta_m"].get<double>(), .01);
			EXPECT_EQ(flight["traces"].size(), run["particles"]);
			EXPECT_FALSE(flight["samples"].get<std::string>().empty());
			for (const auto& trace : flight["traces"]) {
				EXPECT_GT(trace["last"][0].get<double>(), trace["first"][0].get<double>());
				EXPECT_LE(trace["count"].get<std::size_t>(), 96U);
				EXPECT_EQ(trace["end_event"], "ground");
				ASSERT_EQ(trace["first"].size(), 8U);
				ASSERT_EQ(trace["last"].size(), 8U);
				for (const auto& impact : view["impacts"])
					if (impact[5] == run["run_id"] && impact[6] == trace["id"]) {
						EXPECT_NEAR(trace["last"][5].get<double>(), impact[0].get<double>(), 1e-8);
						EXPECT_NEAR(trace["last"][6].get<double>(), impact[1].get<double>(), 1e-8);
						EXPECT_EQ(trace["last"][2], impact[2]);
					}
			}
		}
		GDALAllRegister();
		auto* raster = static_cast<GDALDataset*>(GDALOpen((maps / "density.tif").string().c_str(), GA_ReadOnly));
		ASSERT_NE(raster, nullptr);
		EXPECT_EQ(raster->GetRasterCount(), 2);
		double transform[6]{};
		EXPECT_EQ(raster->GetGeoTransform(transform), CE_None);
		EXPECT_GT(transform[1], 0);
		EXPECT_LT(transform[5], 0);
		std::vector<double> values(raster->GetRasterXSize() * raster->GetRasterYSize());
		for (int i = 1; i <= 2; ++i) {
			EXPECT_EQ(raster->GetRasterBand(i)->RasterIO(
			              GF_Read, 0, 0, raster->GetRasterXSize(), raster->GetRasterYSize(), values.data(),
			              raster->GetRasterXSize(), raster->GetRasterYSize(), GDT_Float64, 0, 0),
			          CE_None);
			EXPECT_NEAR(std::accumulate(values.begin(), values.end(), 0.), 1., 1e-12);
			const double cutoff = *std::max_element(values.begin(), values.end()) * 2e-4;
			for (int y = 0; y < raster->GetRasterYSize(); ++y)
				for (int x = 0; x < raster->GetRasterXSize(); ++x) {
					if (values[y * raster->GetRasterXSize() + x] < cutoff)
						continue;
					const double lon = transform[0] + (x + .5) * transform[1],
					             lat = transform[3] + (y + .5) * transform[5];
					EXPECT_GT(lon, view["bounds"][0][1].get<double>());
					EXPECT_LT(lon, view["bounds"][1][1].get<double>());
					EXPECT_GT(lat, view["bounds"][0][0].get<double>());
					EXPECT_LT(lat, view["bounds"][1][0].get<double>());
				}
		}
		std::ifstream contourStream(maps / "contours.geojson");
		const auto contours = Json::parse(contourStream);
		ASSERT_GE(contours["features"].size(), 6U);
		for (const auto& feature : contours["features"]) {
			const auto& p = feature["geometry"]["coordinates"][0];
			EXPECT_NEAR(p[0].get<double>(), .8, .02);
			EXPECT_NEAR(p[1].get<double>(), 49.8, .02);
			const int band = feature["properties"]["quantity"] == "count" ? 1 : 2;
			ASSERT_EQ(raster->GetRasterBand(band)->RasterIO(
			              GF_Read, 0, 0, raster->GetRasterXSize(), raster->GetRasterYSize(), values.data(),
			              raster->GetRasterXSize(), raster->GetRasterYSize(), GDT_Float64, 0, 0),
			          CE_None);
			const auto nx = raster->GetRasterXSize();
			const double peak = *std::max_element(values.begin(), values.end());
			for (const auto& point : feature["geometry"]["coordinates"]) {
				const double x = (point[0].get<double>() - transform[0]) / transform[1] - .5;
				const double y = (point[1].get<double>() - transform[3]) / transform[5] - .5;
				const int ix = static_cast<int>(std::floor(x)), iy = static_cast<int>(std::floor(y));
				ASSERT_GE(ix, 0);
				ASSERT_GE(iy, 0);
				ASSERT_LT(ix, nx - 1);
				ASSERT_LT(iy, raster->GetRasterYSize() - 1);
				const double u = x - ix, v = y - iy;
				const auto i = iy * nx + ix;
				const double value = (1 - v) * ((1 - u) * values[i] + u * values[i + 1]) +
				                     v * ((1 - u) * values[i + nx] + u * values[i + nx + 1]);
				EXPECT_NEAR(value, feature["properties"]["threshold"].get<double>(), peak * 2e-6);
			}
		}
		GDALClose(raster);
		const auto replay = loadJob(spec.outputDirectory / "job.json");
		ASSERT_TRUE(replay.maps);
		EXPECT_EQ(replay.maps->method, KdeMethod::Fixed);
		EXPECT_EQ(replay.maps->title, "CX1 <test>");
		const auto again = renderSavedRun(spec.outputDirectory, m_directory / "new-map", *spec.maps);
		EXPECT_EQ(again.status, ArtifactStatus::Complete) << again.error;
		EXPECT_THROW(renderSavedRun(spec.outputDirectory, m_directory / "new-map", *spec.maps), std::invalid_argument);
	}

	TEST_F(JobTest, FailedRendererPreservesSuccessfulNumericalOutput)
	{
		auto c = config();
		c["maps"] = {{"method", "fixed"}, {"center_deg", {48., .8}}, {"half_size_km", {1, 1}}};
		const auto spec = load(c);
		const auto result = StrewnJob{}.run(spec);
		EXPECT_EQ(result.simulation.completed, 3U);
		EXPECT_EQ(result.simulation.failed, 0U);
		EXPECT_EQ(result.artifacts.status, ArtifactStatus::Failed);
		EXPECT_FALSE(result.artifacts.error.empty());
		std::ifstream stream(spec.outputDirectory / "manifest.json");
		const auto manifest = Json::parse(stream);
		EXPECT_EQ(manifest["status"], "complete");
		EXPECT_EQ(manifest["artifacts"]["status"], "failed");
		EXPECT_TRUE(std::filesystem::exists(spec.outputDirectory / "impacts.csv"));
	}

	TEST_F(JobTest, ParticleLimitExplainsMissingMapsAndPreservesPartialArchive)
	{
		auto c = config();
		c["input"]["state"]["height_wgs84_m"] = 30000;
		c["input"]["state"]["velocity_enu_mps"] = {9000, 0, -9000};
		c["physics"] = {{"phase", "entry"}, {"max_particles", 1}};
		c["monte_carlo"]["parameters"] = {{"size", {{"mass_kg", 500}}}, {"strong_strength_pa", 2e6}};
		c["maps"] = {{"offline", true}};
		const auto spec = load(c);
		const auto result = StrewnJob{}.run(spec);
		ASSERT_EQ(result.simulation.truncated, 3U);
		EXPECT_EQ(result.simulation.failed, 0U);
		EXPECT_EQ(result.artifacts.status, ArtifactStatus::Skipped);
		EXPECT_NE(result.artifacts.error.find("3 particle-limited"), std::string::npos);
		EXPECT_NE(result.artifacts.error.find("physics.max_particles (currently 1)"), std::string::npos);
		EXPECT_TRUE(result.artifacts.files.empty());
		EXPECT_FALSE(std::filesystem::exists(spec.outputDirectory / "maps"));
		EXPECT_TRUE(std::filesystem::exists(spec.outputDirectory / "breakups.jsonl"));
		std::ifstream stream(spec.outputDirectory / "manifest.json");
		const auto manifest = Json::parse(stream);
		EXPECT_EQ(manifest["status"], "partial");
		EXPECT_EQ(manifest["artifacts"]["status"], "skipped");
		EXPECT_EQ(manifest["artifacts"]["error"], result.artifacts.error);
	}

	TEST_F(JobTest, FailedIntegrationExplainsMissingMapsWithoutParticleLimitAdvice)
	{
		auto c = config();
		c["integration"] = {{"max_flight_time_s", 0.001}};
		c["maps"] = {{"offline", true}};
		const auto spec = load(c);
		const auto result = StrewnJob{}.run(spec);
		ASSERT_EQ(result.simulation.failed, 3U);
		EXPECT_EQ(result.simulation.truncated, 0U);
		EXPECT_EQ(result.artifacts.status, ArtifactStatus::Skipped);
		EXPECT_NE(result.artifacts.error.find("3 failed"), std::string::npos);
		EXPECT_EQ(result.artifacts.error.find("physics.max_particles"), std::string::npos);
		EXPECT_FALSE(std::filesystem::exists(spec.outputDirectory / "maps"));
		std::ifstream stream(spec.outputDirectory / "manifest.json");
		EXPECT_EQ(Json::parse(stream)["artifacts"]["error"], result.artifacts.error);
	}

	TEST_F(JobTest, WallTimeoutIsArchivedAndBlocksMapExport)
	{
		auto c = config();
		c["monte_carlo"]["run_timeout_s"] = 1e-12;
		c["maps"] = {{"offline", true}};
		const auto spec = load(c);
		const auto result = StrewnJob{}.run(spec);
		EXPECT_EQ(result.simulation.failed, 3U);
		EXPECT_EQ(result.simulation.wallTimeouts, 3U);
		EXPECT_EQ(result.artifacts.status, ArtifactStatus::Skipped);
		EXPECT_NE(result.artifacts.error.find("monte_carlo.run_timeout_s"), std::string::npos);
		EXPECT_FALSE(std::filesystem::exists(spec.outputDirectory / "maps"));
		std::ifstream summary(spec.outputDirectory / "summary.json");
		EXPECT_EQ(Json::parse(summary)["wall_timeouts"], 3);
		std::ifstream runs(spec.outputDirectory / "runs.jsonl");
		std::string line;
		ASSERT_TRUE(static_cast<bool>(std::getline(runs, line)));
		EXPECT_EQ(Json::parse(line)["status"], "wall_timeout");
		const auto replay = loadJob(spec.outputDirectory / "job.json");
		EXPECT_EQ(replay.monteCarlo.runTimeoutS, 1e-12);
	}

	TEST_F(JobTest, RejectsInvalidWallTimeoutAndDefaultsToUnlimitedParticleCount)
	{
		const auto defaults = load(config());
		EXPECT_EQ(defaults.physics.maxParticles, 0U);
		EXPECT_EQ(defaults.monteCarlo.runTimeoutS, 300);
		for (const auto& value : std::vector<Json>{0, -1, 86401, "300", true, nullptr}) {
			auto c = config();
			c["monte_carlo"]["run_timeout_s"] = value;
			EXPECT_THROW(load(c), std::exception);
		}
		auto c = config();
		c["physics"]["max_particles"] = 5000;
		EXPECT_EQ(load(c).physics.maxParticles, 5000U);
		c["physics"]["max_particles"] = 0;
		EXPECT_EQ(load(c).physics.maxParticles, 0U);
	}

	TEST_F(JobTest, CancelledCalculationRecordsCancelledMapExport)
	{
		auto c = config();
		c["maps"] = {{"offline", true}};
		const auto spec = load(c);
		std::stop_source cancel;
		const auto result = StrewnJob{}.run(spec, {cancel.get_token(), [&](const auto&) { cancel.request_stop(); }});
		ASSERT_TRUE(result.simulation.cancelled);
		EXPECT_EQ(result.artifacts.status, ArtifactStatus::Cancelled);
		EXPECT_FALSE(result.artifacts.error.empty());
		EXPECT_FALSE(std::filesystem::exists(spec.outputDirectory / "maps"));
		std::ifstream stream(spec.outputDirectory / "manifest.json");
		const auto manifest = Json::parse(stream);
		EXPECT_EQ(manifest["status"], "cancelled");
		EXPECT_EQ(manifest["artifacts"]["status"], "cancelled");
	}

	TEST_F(JobTest, DiagnosticReplayMismatchDoesNotInvalidateProbabilityMaps)
	{
		const auto spec = load(config());
		ASSERT_EQ(StrewnJob{}.run(spec).simulation.failed, 0U);
		Json run;
		{
			std::ifstream in(spec.outputDirectory / "runs.jsonl");
			std::string line;
			std::getline(in, line);
			run = Json::parse(line);
		}
		run["parameters"]["mass_kg"] = 10.;
		{
			std::ofstream out(spec.outputDirectory / "runs.jsonl");
			out << run.dump() << '\n';
		}
		MapOptions options;
		options.offline = true;
		options.method = KdeMethod::Fixed;
		const auto maps = m_directory / "diagnostic-mismatch";
		const auto report = renderSavedRun(spec.outputDirectory, maps, options);
		ASSERT_EQ(report.status, ArtifactStatus::Complete) << report.error;
		EXPECT_TRUE(std::filesystem::exists(maps / "density.tif"));
		EXPECT_TRUE(std::filesystem::exists(maps / "field.html"));
		EXPECT_FALSE(std::filesystem::exists(maps / "fragment_time_series.png"));
		EXPECT_TRUE(std::any_of(report.warnings.begin(), report.warnings.end(), [](const auto& message) {
			return message.starts_with("Ensemble diagnostics unavailable:");
		}));
	}

	TEST_F(JobTest, RejectsScottAndInvalidMapOptions)
	{
		for (const auto& bad : std::vector<Json>{{{"method", "scott"}},
		                                         {{"method", "unknown"}},
		                                         {{"cell_m", 0}},
		                                         {{"adaptive_classes", 1.5}},
		                                         {{"image_width", -1}},
		                                         {{"mystery", true}}}) {
			auto c = config();
			c["maps"] = bad;
			EXPECT_THROW(load(c), std::invalid_argument);
		}
	}

	TEST_F(JobTest, RendererCancellationAndIncompleteArchiveAreExplicit)
	{
		const auto spec = load(config());
		ASSERT_EQ(StrewnJob{}.run(spec).simulation.failed, 0U);
		std::stop_source stop;
		stop.request_stop();
		const auto cancelled =
		    renderSavedRun(spec.outputDirectory, m_directory / "cancelled-map", {}, stop.get_token());
		EXPECT_EQ(cancelled.status, ArtifactStatus::Cancelled);
		EXPECT_FALSE(std::filesystem::exists(m_directory / "cancelled-map"));
		{
			std::ofstream out(spec.outputDirectory / "manifest.json");
			out << Json{{"status", "partial"}};
		}
		const auto failed = renderSavedRun(spec.outputDirectory, m_directory / "partial-map", {});
		EXPECT_EQ(failed.status, ArtifactStatus::Failed);
		EXPECT_FALSE(std::filesystem::exists(m_directory / "partial-map" / "impact_probability.png"));
	}

	TEST_F(JobTest, EmptyImpactArchiveProducesNoProbabilityMap)
	{
		const auto spec = load(config());
		ASSERT_EQ(StrewnJob{}.run(spec).simulation.failed, 0U);
		std::string header;
		{
			std::ifstream in(spec.outputDirectory / "impacts.csv");
			std::getline(in, header);
		}
		{
			std::ofstream out(spec.outputDirectory / "impacts.csv");
			out << header << '\n';
		}
		Json summary;
		{
			std::ifstream in(spec.outputDirectory / "summary.json");
			summary = Json::parse(in);
		}
		summary["tracked_impacts"] = 0;
		{
			std::ofstream out(spec.outputDirectory / "summary.json");
			out << summary;
		}
		const auto result = renderSavedRun(spec.outputDirectory, m_directory / "empty-map", {});
		EXPECT_EQ(result.status, ArtifactStatus::Empty) << result.error;
		EXPECT_FALSE(std::filesystem::exists(m_directory / "empty-map" / "impact_probability.png"));
		EXPECT_TRUE(std::filesystem::exists(m_directory / "empty-map" / "artifacts.json"));
		EXPECT_TRUE(std::filesystem::exists(m_directory / "empty-map" / "field.html"));
	}

	TEST_F(JobTest, HtmlEscapesInputAndNeedsNoFontOrExternalChartAssets)
	{
		auto c = config();
		c["maps"] = {
		    {"title", "</script><script>alert(1)</script>"}, {"font_file", "missing.ttf"}, {"method", "fixed"}};
		const auto spec = load(c);
		const auto result = StrewnJob{}.run(spec);
		ASSERT_EQ(result.artifacts.status, ArtifactStatus::Complete) << result.artifacts.error;
		std::ifstream in(spec.outputDirectory / "maps" / "field.html", std::ios::binary);
		const std::string html(std::istreambuf_iterator<char>(in), {});
		EXPECT_EQ(html.find("</script><script>alert(1)</script>"), std::string::npos);
		EXPECT_EQ(html.find("__PRINT_ICON__"), std::string::npos);
		EXPECT_EQ(html.find("__PRINT_JS__"), std::string::npos);
		EXPECT_EQ(html.find("__TRAJECTORY_JS__"), std::string::npos);
		EXPECT_EQ(html.find("__MAP_FLIGHT_JS__"), std::string::npos);
		EXPECT_EQ(html.find("__PLAY_ICON__"), std::string::npos);
		EXPECT_NE(html.find("StrewnMapFlight.create"), std::string::npos);
		EXPECT_NE(html.find("id=\"print-pages\""), std::string::npos);
		EXPECT_NE(html.find("StrewnPrint.prepare"), std::string::npos);
		EXPECT_EQ(html.find("__DATA__"), std::string::npos);
		const auto data = reportData(html);
		EXPECT_EQ(data["title"], c["maps"]["title"]);
		EXPECT_FALSE(data["options"].contains("font_file"));
		EXPECT_FALSE(data["options"].contains("image_width"));
	}

	TEST_F(JobTest, ReplayMismatchIsIsolatedToOneRealization)
	{
		const auto spec = load(config());
		ASSERT_EQ(StrewnJob{}.run(spec).simulation.failed, 0U);
		std::vector<Json> rows;
		{
			std::ifstream in(spec.outputDirectory / "runs.jsonl");
			std::string line;
			while (std::getline(in, line))
				rows.push_back(Json::parse(line));
		}
		rows[1]["parameters"]["mass_kg"] = 10.;
		{
			std::ofstream out(spec.outputDirectory / "runs.jsonl");
			for (const auto& row : rows)
				out << row.dump() << '\n';
		}
		const auto output = m_directory / "isolated-replay";
		MapOptions options;
		options.method = KdeMethod::Fixed;
		const auto result = renderSavedRun(spec.outputDirectory, output, options);
		ASSERT_EQ(result.status, ArtifactStatus::Complete) << result.error;
		std::ifstream in(output / "field.html", std::ios::binary);
		const std::string html(std::istreambuf_iterator<char>(in), {});
		const auto data = reportData(html);
		EXPECT_TRUE(data["runs"][0]["flight"]["verified"].get<bool>());
		EXPECT_FALSE(data["runs"][1]["flight"]["verified"].get<bool>());
		EXPECT_TRUE(data["runs"][2]["flight"]["verified"].get<bool>());
		EXPECT_FALSE(data["runs"][1]["flight"].contains("samples"));
	}

}  // namespace
