#include "orbital_core/hazard/MonteCarloArchive.h"
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "NBodyPropagator.h"
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <chrono>
#include <fstream>

namespace {
	namespace hz = orbital_core::hazard;
	using Json = nlohmann::json;

	class MonteCarloArchiveTest : public testing::Test {
	   protected:
		std::filesystem::path m_root;
		hz::MonteCarloArchiveInput m_input;
		hz::ObsMcResult m_result;
		hz::CloseApproachEvent m_event;
		NBodyPropagator m_prop;

		void SetUp() override
		{
			m_root = std::filesystem::temp_directory_path() /
			         ("mc-archive-test-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
			std::filesystem::create_directory(m_root);
			m_input.solution.objectName = "TC3-test";
			m_input.solution.nominal.setEpoch(2454746.5);
			m_input.solution.nominal.setR({1, 0, 0});
			m_input.solution.nominal.setV({0, 1, 0});
			m_input.solution.covariance = orbital_core::StateCov6x6::Identity() * 1e-18;
			m_input.solution.absoluteMagnitude = 30.311;
			m_input.albedo = 0.15;
			m_input.covarianceObjectiveScale = 2.5;
			m_input.sampling.numClones = 2;
			m_input.sampling.altitudeKm = 100;
			m_event.planet = hz::Planet::Earth;
			m_event.jd = 2454747;
			m_result.earthFrameModel = hz::earthBodyFixedFrameModel();
			m_result.altitudeKm = 100;
			m_result.numAttempted = m_result.numConverged = m_result.numAnyImpactClones = 2;
			m_result.aggregateImpactProb = 1;
			m_result.hitsPerEvent = {2};
			for (int i = 0; i < 2; ++i) {
				hz::ImpactCloneInfo point;
				point.cloneIndex = i;
				point.eventIndex = 0;
				point.planet = hz::Planet::Earth;
				point.epochMjdTdb = 54746.5;
				point.epochMjdUtc = 54746.4992;
				point.bodyFixedValid = true;
				point.rBodyFixedKm = {Constants::EARTH_EQ_RADIUS_KM + 100, 0, 0};
				point.vBodyFixedKmS = {-12, 0, 0};
				m_result.corridor.push_back(point);
			}
		}

		void TearDown() override
		{
			if (std::filesystem::canonical(m_root).parent_path() ==
			        std::filesystem::canonical(std::filesystem::temp_directory_path()) &&
			    m_root.filename().string().starts_with("mc-archive-test-"))
				std::filesystem::remove_all(m_root);
		}

		std::filesystem::path save()
		{
			return hz::saveMonteCarloArchive(m_root, m_input, m_prop, {&m_event, 1}, m_result).directory;
		}

		static Json read(const std::filesystem::path& path)
		{
			std::ifstream in(path);
			return Json::parse(in);
		}
	};

	TEST_F(MonteCarloArchiveTest, ExactScaledCovarianceAndMjdEntriesNoAutoLaunch)
	{
		const auto directory = save();
		const auto request = read(directory / "request.json");
		EXPECT_EQ(request["covariance"]["matrix"][0][0].get<double>(), m_input.solution.covariance(0, 0));
		EXPECT_EQ(request["covariance"]["objective_scale_already_applied"], 2.5);
		EXPECT_EQ(request["nominal"]["epoch_mjd_tdb"], m_input.solution.nominal.epochMjd());
		EXPECT_EQ(request["nominal"]["velocity_au_day"][1], Constants::GAUSS_K);
		const auto entries = read(directory / "entries.json");
		EXPECT_EQ(entries["entries"][0]["epoch_mjd_utc"], m_result.corridor[0].epochMjdUtc);
		EXPECT_EQ(entries["entries"][0]["velocity_ecef_mps"][0], -12000.0);
		EXPECT_EQ(entries["frame_model"], m_result.earthFrameModel);
		for (const auto& entry : entries["entries"])
			EXPECT_FALSE(entry.contains("weight"));
		const auto manifest = read(directory / "manifest.json");
		EXPECT_EQ(manifest["entry_states_ready"], true);
		EXPECT_EQ(manifest["atmospheric_input_ready"], true);
		EXPECT_EQ(manifest["automatic_launch"], false);
		EXPECT_FALSE(manifest.contains("probability_scope"));
		EXPECT_NE(directory, save());
	}

	TEST_F(MonteCarloArchiveTest, ExportsRecordedFrameWithoutSubstitutingCurrentModel)
	{
		m_result.earthFrameModel = "test frame: IAU2006/2000A; polar motion: applied";
		EXPECT_EQ(read(save() / "entries.json")["frame_model"], m_result.earthFrameModel);
		m_result.earthFrameModel.clear();
		EXPECT_EQ(read(save() / "entries.json")["frame_model"], "unknown (not recorded)");
	}

	TEST_F(MonteCarloArchiveTest, FrameDescriptionDistinguishesAppliedChainFromEopProfile)
	{
		const auto frame = hz::earthBodyFixedFrameModel();
		EXPECT_NE(frame.find("IAU1976/1980"), std::string::npos);
		EXPECT_NE(frame.find("polar motion: not applied"), std::string::npos);
		EXPECT_NE(frame.find("EOP profile:"), std::string::npos);
		EXPECT_NE(frame.find("UTC fallback"), std::string::npos);
	}

	TEST_F(MonteCarloArchiveTest, FullNonGravCrossCovarianceOrder)
	{
		auto& ng = m_input.solution.nonGrav;
		ng.setModel(NgModel::Yarkovsky);
		ng.setEnabled(true);
		ng.setFitParameterMask(NgParameterMask::A2);
		ASSERT_EQ(ng.numFitParams(), 1);
		m_input.kind = hz::CovarianceSampling::Parameters;
		m_input.solution.covariance = Eigen::MatrixXd::Identity(7, 7) * 1e-18;
		m_input.solution.covariance(0, 6) = m_input.solution.covariance(6, 0) = 1e-20;
		const auto directory = save();
		const auto json = read(directory / "request.json");
		EXPECT_EQ(json["covariance"]["parameter_order"][6], "a2_au_day2");
		EXPECT_EQ(json["covariance"]["matrix"][0][6], 1e-20);
		EXPECT_EQ(read(directory / "manifest.json")["entry_states_ready"], true);
	}

	TEST_F(MonteCarloArchiveTest, RoundedProbabilityDoesNotGrantReadiness)
	{
		m_input.sampling.numClones = 1000;
		EXPECT_EQ(read(save() / "manifest.json")["entry_states_ready"], false);
	}

	TEST_F(MonteCarloArchiveTest, PartialCancelledAndWrongBoundaryAreArchivedButBlocked)
	{
		m_result.cancelled = true;
		EXPECT_EQ(read(save() / "manifest.json")["status"], "cancelled");
		m_result.cancelled = false;
		m_result.numFailedFit = 1;
		EXPECT_EQ(read(save() / "manifest.json")["entry_states_ready"], false);
		m_result.numFailedFit = 0;
		m_input.sampling.altitudeKm = m_result.altitudeKm = 0;
		EXPECT_EQ(read(save() / "manifest.json")["entry_states_ready"], false);
	}

	TEST_F(MonteCarloArchiveTest, RejectsDuplicateOutboundAndWrongHeight)
	{
		m_result.corridor[1].cloneIndex = 0;
		EXPECT_EQ(read(save() / "manifest.json")["entry_states_ready"], false);
		m_result.corridor[1].cloneIndex = 1;
		m_result.corridor[0].vBodyFixedKmS.x = 12;
		EXPECT_EQ(read(save() / "manifest.json")["entry_states_ready"], false);
		m_result.corridor[0].vBodyFixedKmS.x = -12;
		m_result.corridor[0].rBodyFixedKm.x += 1;
		EXPECT_EQ(read(save() / "manifest.json")["entry_states_ready"], false);
	}

	TEST_F(MonteCarloArchiveTest, InvalidCovarianceAndMissingSizeAreExplicit)
	{
		m_input.solution.covariance(0, 0) = -1;
		EXPECT_EQ(read(save() / "manifest.json")["entry_states_ready"], false);
		m_input.solution.covariance(0, 0) = 1e-18;
		m_input.albedo.reset();
		const auto manifest = read(save() / "manifest.json");
		EXPECT_EQ(manifest["entry_states_ready"], true);
		EXPECT_EQ(manifest["atmospheric_input_ready"], false);
	}
}  // namespace
