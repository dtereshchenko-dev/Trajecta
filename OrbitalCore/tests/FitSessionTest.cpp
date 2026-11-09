#include <gtest/gtest.h>

#include "orbital_core/FitSession.h"
#include "orbital_core/FitPipeline.h"
#include "fit/FitObservations.h"
#include "fit/OrbitFitProblem.h"
#include "fit/OrbitFitUncertainty.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/optimizers/EigenLevenbergMarquardt.h"
#include "fit/FitPipelineDetail.h"
#include "fit/ObservationSelection.h"
#include "fit/HystereticRejector.h"
#include "fit/OppositionCoverage.h"
#include "orbital_core/SolutionRefinement.h"

#include "orbital_core/ObservationSet.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/SpkReader.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/StateVectors.h"
#include "AstrometryBias.h"
#include "InitialOrbit.h"
#include "IODExtender.h"
#include "OrbitOptimizer.h"
#include "logger/Logger.h"

#include <array>
#include <algorithm>
#include <filesystem>
#include <limits>
#include <memory>
#include <type_traits>
#include <Eigen/QR>
#include <Eigen/SVD>

namespace oc = orbital_core;

constexpr auto solverStart = "[EigenLM] start";
constexpr auto solverFinish = "[EigenLM] finish";

TEST(ObservationSetTest, RejectsEmptyAndMissingPaths)
{
	EXPECT_THROW((void)oc::loadObservations({}), std::invalid_argument);
	const auto missing =
	    std::filesystem::path(__FILE__).parent_path() / "data" / "nonexistent-observation-input-for-test.xml";
	ASSERT_FALSE(std::filesystem::exists(missing));
	EXPECT_THROW((void)oc::loadObservations(missing), std::runtime_error);
}

TEST(FitRequestTest, RejectsInvalidSelectionBeforePropagation)
{
	const oc::FitSession session{oc::FitConfig{}};
	const auto orbit = std::make_shared<Orbit>();
	const std::vector<ObsNode> observations(4);
	const std::array<std::uint8_t, 3> shortMask{1, 1, 1};
	const std::array<std::uint8_t, 4> sparseMask{1, 0, 1, 0};
	EXPECT_THROW(
	    (void)session.optimize({.initialOrbit = orbit, .observations = observations, .initialIncluded = shortMask}),
	    std::invalid_argument);
	EXPECT_THROW(
	    (void)session.optimize({.initialOrbit = orbit, .observations = observations, .initialIncluded = sparseMask}),
	    std::invalid_argument);
}

TEST(FitRequestTest, PreCancelledRequestDoesNoWork)
{
	std::stop_source cancel;
	cancel.request_stop();
	oc::FitSession session{oc::FitConfig{}};
	int progressCalls = 0;
	const auto result =
	    session.optimize({}, {.cancel = cancel.get_token(), .progress = [&](const auto&) { ++progressCalls; }});
	EXPECT_TRUE(result.cancelled);
	EXPECT_FALSE(result.converged);
	EXPECT_FALSE(result.orbit);
	EXPECT_FALSE(result.fitEpochOrbit);
	EXPECT_EQ(result.residuals.count(), 0u);
	EXPECT_FALSE(result.uncertainty.hasElementCovariance);
	EXPECT_EQ(progressCalls, 0);
}

TEST(FitRequestTest, RejectsInvalidOutputEpochBeforePropagation)
{
	const oc::FitSession session{oc::FitConfig{}};
	const auto orbit = std::make_shared<Orbit>();
	const std::vector<ObsNode> observations(3);
	for (double epoch :
	     {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
		EXPECT_THROW(
		    (void)session.optimize({.initialOrbit = orbit, .observations = observations, .outputEpochJdTdb = epoch}),
		    std::invalid_argument);
	}
}

TEST(FitPipelineTest, RejectsIncompleteObservationPools)
{
	const oc::ObservationSet shortAll{std::vector<ObsNode>(2), std::vector<ObsNode>(3)};
	const oc::ObservationSet shortIod{std::vector<ObsNode>(3), std::vector<ObsNode>(2)};
	const oc::FitPipeline pipeline{oc::FitConfig{}};
	EXPECT_THROW((void)pipeline.run({shortAll, {}}), std::invalid_argument);
	EXPECT_THROW((void)pipeline.run({shortIod, {}}), std::invalid_argument);
}

TEST(FitPipelineTest, RejectsInvalidOutputEpochBeforeIod)
{
	const oc::ObservationSet observations{std::vector<ObsNode>(3), std::vector<ObsNode>(3)};
	const oc::FitPipeline pipeline{oc::FitConfig{}};
	for (double epoch :
	     {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
		EXPECT_THROW((void)pipeline.run({observations, epoch}), std::invalid_argument);
	}
}

TEST(FitPipelineTest, CancellationPreventsIodAndCallbacks)
{
	const oc::ObservationSet observations{std::vector<ObsNode>(3), std::vector<ObsNode>(3)};
	std::stop_source cancel;
	cancel.request_stop();
	int calls = 0;
	const auto result = oc::FitPipeline{oc::FitConfig{}}.run(
	    {observations, {}}, {.cancel = cancel.get_token(), .progress = [&](const auto&) { ++calls; }},
	    [&](const auto&) { ++calls; });
	EXPECT_TRUE(result.fit.cancelled);
	EXPECT_FALSE(result.fit.orbit);
	EXPECT_FALSE(result.initial.orbit);
	EXPECT_FALSE(result.extension);
	EXPECT_EQ(result.fitInputObservations, 0u);
	EXPECT_TRUE(result.observations.all.empty());
	EXPECT_TRUE(result.observations.iod.empty());
	EXPECT_EQ(calls, 0);
}

TEST(FitPipelineTest, FittedObservationViewRetainsRejectedRows)
{
	oc::FitPipelineResult result;
	EXPECT_TRUE(result.fittedObservations().empty());
	result.observations.all.resize(4);
	for (std::size_t i = 0; i < result.observations.all.size(); ++i)
		result.observations.all[i].JD = 2450000.0 + i;
	const auto all = result.fittedObservations();
	static_assert(std::is_same_v<std::ranges::range_reference_t<decltype(all)>, const ObsNode&>);
	ASSERT_EQ(all.size(), 4u);
	EXPECT_EQ(&all[2], &result.observations.all[2]);

	result.fit.keptObsIndices = {1, 3};
	const auto selected = result.fittedObservations();
	ASSERT_EQ(selected.size(), 2u);
	EXPECT_EQ(&selected[0], &result.observations.all[1]);
	EXPECT_EQ(&selected[1], &result.observations.all[3]);
	EXPECT_EQ(result.observations.all.size(), 4u);
	EXPECT_DOUBLE_EQ(result.observations.all[0].JD, 2450000.0);
}

TEST(FitPipelineTest, FittedObservationViewUsesItsOwnerAndChecksIndices)
{
	oc::FitPipelineResult original;
	original.observations.all.resize(3);
	original.observations.all[2].JD = 2450002.0;
	original.fit.keptObsIndices = {0, 2};
	auto copy = original;
	EXPECT_NE(&copy.fittedObservations()[1], &original.fittedObservations()[1]);
	EXPECT_EQ(&copy.fittedObservations()[1], &copy.observations.all[2]);
	auto moved = std::move(copy);
	EXPECT_EQ(&moved.fittedObservations()[1], &moved.observations.all[2]);
	EXPECT_DOUBLE_EQ(moved.fittedObservations()[1].JD, 2450002.0);

	moved.fit.keptObsIndices = {3};
	EXPECT_THROW((void)moved.fittedObservations().front(), std::out_of_range);
}

TEST(FitPipelineTest, ProgressCanCancelBeforeIodOrPropagateItsException)
{
	const oc::ObservationSet observations{std::vector<ObsNode>(3), std::vector<ObsNode>(3)};
	const oc::FitPipeline pipeline{oc::FitConfig{}};
	std::stop_source cancel;
	int initialCalls = 0;
	const auto result = pipeline.run({observations, {}},
	                                 {.cancel = cancel.get_token(),
	                                  .progress =
	                                      [&](const std::string& message) {
		                                      EXPECT_EQ(message, "Starting IOD");
		                                      cancel.request_stop();
	                                      }},
	                                 [&](const auto&) { ++initialCalls; });
	EXPECT_TRUE(result.fit.cancelled);
	EXPECT_FALSE(result.initial.orbit);
	EXPECT_EQ(initialCalls, 0);

	EXPECT_THROW((void)pipeline.run({observations, {}},
	                                {.progress = [](const auto&) { throw std::logic_error("callback failure"); }}),
	             std::logic_error);
}

TEST(FitPipelineTest, MapsIodArcIncludingBoundaryTimestampTies)
{
	std::vector<ObsNode> all(6), iod(3);
	const std::array epochs{1.0, 2.0, 2.0, 3.0, 4.0, 4.0};
	for (std::size_t i = 0; i < all.size(); ++i)
		all[i].JD = epochs[i];
	iod[0].JD = 2.0;
	iod[1].JD = 3.0;
	iod[2].JD = 4.0;
	const auto arc = oc::fit::mapIodSubarc(iod, all, {0, 2});
	ASSERT_TRUE(arc);
	EXPECT_EQ(arc->first, 1u);
	EXPECT_EQ(arc->last, 5u);
	const auto singleEpoch = oc::fit::mapIodSubarc(iod, all, {0, 0});
	ASSERT_TRUE(singleEpoch);
	EXPECT_EQ(singleEpoch->first, 1u);
	EXPECT_EQ(singleEpoch->last, 2u);
}

TEST(FitPipelineTest, RejectsInvalidOrUnmappedIodArc)
{
	std::vector<ObsNode> all(2), iod(2);
	all[0].JD = 1.0;
	all[1].JD = 4.0;
	iod[0].JD = 2.0;
	iod[1].JD = 3.0;
	EXPECT_FALSE(oc::fit::mapIodSubarc(iod, all, {1, 0}));
	EXPECT_FALSE(oc::fit::mapIodSubarc(iod, all, {0, 2}));
	EXPECT_FALSE(oc::fit::mapIodSubarc(iod, all, {0, 1}));
	EXPECT_FALSE(oc::fit::mapIodSubarc(iod, {}, {0, 1}));
}

namespace {
	[[nodiscard]] bool resolveDataDir()
	{
		if (!getExecutableDir().empty())
			return true;

		auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
		const std::filesystem::path candidates[] = {
		    srcDir.parent_path() / "data",
		    srcDir / "data",
		    srcDir / "release",
		    srcDir / "debug",
		};
		for (const auto& dir : candidates) {
			if (std::filesystem::exists(dir / "linux_p1550p2650.440")) {
				getExecutableDir() = dir;
				return true;
			}
		}
		return false;
	}
}  // namespace

class FitSessionTest : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (!resolveDataDir())
			GTEST_SKIP() << "data/ directory not found";

		m_de440 = std::make_unique<DE440Source>();
		if (!m_de440->isAvailable())
			GTEST_SKIP() << "DE440 kernel not loadable";
		setDefaultEphemerisSource(m_de440.get());

		m_sb441 = std::make_unique<SpkReader>();
		if (m_sb441->isAvailable())
			setDefaultAsteroidSource(m_sb441.get());
	}

	void TearDown() override
	{
		setDefaultEphemerisSource(nullptr);
		setDefaultAsteroidSource(nullptr);
	}

	std::unique_ptr<DE440Source> m_de440;
	std::unique_ptr<SpkReader> m_sb441;
};

namespace {
	struct FitContractCase {
		oc::ObservationSet observations;
		Orbit_p initial;
		oc::FitConfig config;
	};

	FitContractCase fitContractCase()
	{
		FitContractCase data;
		data.observations = oc::loadObservations(std::filesystem::path(__FILE__).parent_path() /
		                                         "data/2020VV_ADES_2022.ESA11_table_replay.xml");
		data.config.nbody.enabled = true;
		data.config.useWeights = true;
		data.config.weightProfile = oc::AstrometricWeightProfile::Vfcc17;
		data.config.obsWeighting.useObserverSigmas = true;
		oc::prepareObservations(data.observations, {.debias = true, .weightProfile = data.config.weightProfile});
		data.initial = make_orbit_from_states_vectors(
		    StateVectors(Vec3D(.62434444840234238, .83392068676583608, .0034709928450530553),
		                 Vec3D(-.012853569243756456, .011750580852295998, .000093051552130828674) / Constants::GAUSS_K),
		    2459170.2715327661);
		return data;
	}

	struct FitLogCallback {
		explicit FitLogCallback(Logger::Callback callback) { Logger::setCallback(std::move(callback)); }

		~FitLogCallback() { Logger::setCallback({}); }
	};
}  // namespace

TEST_F(FitSessionTest, OpticalFitJacobianAndCovarianceMatchIndependentResidualPerturbations)
{
	for (int scenario = 0; scenario < 3; ++scenario) {
		SCOPED_TRACE(scenario);
		OrbitOptimizerConfig config;
		config.nbody.enabled = true;
		config.nbody.useSb441 = true;
		config.useWeights = true;
		config.obsWeighting.useObserverSigmas = true;
		config.obsWeighting.timingSigmaDefaultSeconds = 0;
		config.logDiagnostics = config.logCacheStats = false;
		auto orbit = std::make_shared<Orbit>();
		std::vector<ObsNode> observations;
		if (scenario == 0) {
			auto data = oc::loadObservations(std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() /
			                                 "data/examples/NEO/impacts/2014_AA.txt");
			oc::prepareObservations(data, {.debias = true});
			observations = std::move(data.all);
			orbit->setEpochMjd(56658.28746759193);
			orbit->setR({-.18023297527792959, .9692391793977728, -.0004468935734315787});
			orbit->setV(Vec3D(-.01757246516386107, -.006131033968496691, .00046305301202514553) / Constants::GAUSS_K);
		} else {
			oc::cometary::Vector6 cot;
			cot << 1.3, .6, .41, .82, .37, -.7;
			const auto state = oc::cometary::toCartesian(cot).value();
			orbit->setEpochMjd(60000);
			orbit->setR({state[0], state[1], state[2]});
			orbit->setV(Vec3D(state[3], state[4], state[5]) / Constants::GAUSS_K);
			for (int i = 0; i < 24; ++i) {
				ObsNode o;
				o.JD = 2460000.5 - 120 + i * 10;
				o.timeScale = oc::time::Scale::TDB;
				o.RA = 1;
				o.DEC = .3;
				o.obsCode = "500";
				o.hasObserverSigma = true;
				o.rmsRaArcsec = .5;
				o.rmsDecArcsec = .7;
				o.rmsCorr = .3;
				observations.push_back(o);
			}
			config.nonGrav.enabled = scenario == 2;
			config.nonGrav.model = oc::NonGravModelKind::InverseSquare;
			config.nonGrav.solveMask = static_cast<std::uint8_t>(NgParameterMask::A1 | NgParameterMask::A2);
		}
		oc::fit::OrbitFitModel model(orbit, oc::fit::prepareFitObservations(observations, config), config);
		auto params = oc::fit::OrbitFitModel::stateVectorFromOrbit(*orbit);
		if (config.nonGrav.enabled) {
			params.push_back(1e-10);
			params.push_back(-3e-12);
		}
		const auto analytic = model.computeResidualsAndStateJac(params);
		ASSERT_EQ(analytic.stateJ.cols(), params.size());
		const auto nominal = model.computeResiduals(params);
		for (std::size_t i = 0; i < nominal.size(); ++i)
			EXPECT_NEAR(nominal[i], analytic.residuals[i], 2e-12);
		oc::fit::OrbitFitProblem problem(model);
		problem.beginPass(static_cast<int>(params.size()));
		Eigen::MatrixXd weighted = analytic.stateJ;
		problem.weightJacobian(weighted, analytic.transforms);
		Eigen::MatrixXd finite(weighted.rows(), weighted.cols());
		for (int column = 0; column < weighted.cols(); ++column) {
			SCOPED_TRACE(column);
			const double h = column < 3 ? 1e-5 : column < 6 ? 1e-4 : 1e-9;
			std::array<Eigen::VectorXd, 3> difference;
			for (int k = 0; k < 3; ++k) {
				const double step = h / std::pow(2., k);
				auto plus = params, minus = params;
				plus[column] += step;
				minus[column] -= step;
				auto rp = model.computeResiduals(plus), rm = model.computeResiduals(minus);
				problem.weightResiduals(rp, analytic.transforms);
				problem.weightResiduals(rm, analytic.transforms);
				difference[k].resize(weighted.rows());
				for (int row = 0; row < weighted.rows(); ++row)
					difference[k](row) = (rp[row] - rm[row]) / (2 * step);
			}
			const Eigen::VectorXd first = (4 * difference[1] - difference[0]) / 3;
			finite.col(column) = (4 * difference[2] - difference[1]) / 3;
			ASSERT_LT((first - finite.col(column)).norm() / finite.col(column).norm(), 2e-7);
			EXPECT_LT((weighted.col(column) - finite.col(column)).norm() / finite.col(column).norm(), 2e-7);
		}
		problem.saveLinearizationFor(params);
		oc::fit::OrbitFitUncertainty uncertainty(model, problem, orbit->getEpoch());
		const auto covariance = uncertainty.computeFitCovariance(params, analytic.residuals);
		Eigen::VectorXd scale(finite.cols());
		for (int i = 0; i < finite.cols(); ++i)
			scale(i) = 1 / finite.col(i).norm();
		const Eigen::MatrixXd scaled = finite * scale.asDiagonal();
		const Eigen::JacobiSVD<Eigen::MatrixXd> svd(scaled, Eigen::ComputeThinV);
		const Eigen::MatrixXd reference = scale.asDiagonal() * svd.matrixV() *
		                                  svd.singularValues().array().square().inverse().matrix().asDiagonal() *
		                                  svd.matrixV().transpose() * scale.asDiagonal();
		for (int i = 0; i < reference.rows(); ++i)
			for (int j = 0; j < reference.cols(); ++j)
				EXPECT_LT(std::abs(covariance.formalFullCovariance(i, j) - reference(i, j)) /
				              std::sqrt(reference(i, i) * reference(j, j)),
				          scenario == 0 ? 2e-5 : 1e-6);
	}
}

TEST_F(FitSessionTest, ExtensionDoesNotRepeatConvergedFitsForUnchangedForceModel)
{
	for (const auto mode :
	     {Calculations::PerturberSelectionMode::ManualExact, Calculations::PerturberSelectionMode::AutoByArcLength}) {
		SCOPED_TRACE(static_cast<int>(mode));
		const auto data = fitContractCase();
		OrbitOptimizerConfig config;
		config.nbody = data.config.nbody;
		config.nbody.planetMask = 2303;
		config.useWeights = true;
		config.weightProfile = data.config.weightProfile;
		config.obsWeighting = data.config.obsWeighting;
		config.perturberSelectionMode = mode;
		const auto prepared = oc::fit::prepareFitObservations(data.observations.all, config);
		int passes = 0;
		bool lastFitConverged = false;
		bool refitRequired = false;
		bool autoReject = false;
		const FitLogCallback log([&](LogLevel, const std::string& message) {
			if (message.starts_with("[ATT-EXT] auto-reject:"))
				autoReject = true;
			if (message.starts_with(solverFinish))
				lastFitConverged = message.find("converged=1") != std::string::npos;
			if (!autoReject && message.starts_with("[ATT-EXT] iter ") &&
			    message.find(" refitRequired=") != std::string::npos) {
				EXPECT_TRUE(lastFitConverged);
				if (passes > 0)
					EXPECT_TRUE(refitRequired);
				refitRequired = message.find("refitRequired=1") != std::string::npos;
				++passes;
			}
		});
		const auto result = Calculations::attemptExtensions(data.initial, data.observations.all, 0,
		                                                    data.observations.all.size() - 1, config, prepared);
		EXPECT_TRUE(result.converged);
		EXPECT_EQ(passes, 1);
		EXPECT_FALSE(refitRequired);
		EXPECT_TRUE(Calculations::spansInputRange(result, data.observations.all.size()));
	}
}

TEST_F(FitSessionTest, PreparationCopyIsIndependentOfPreviousSettings)
{
	const auto raw = oc::loadObservations(std::filesystem::path(__FILE__).parent_path() /
	                                      "data/2020VV_ADES_2022.ESA11_table_replay.xml");
	auto previous = raw;
	for (const auto profile : {oc::AstrometricWeightProfile::Vfcc17, oc::AstrometricWeightProfile::Uniform,
	                           oc::AstrometricWeightProfile::Vfcc17}) {
		for (const bool debias : {true, false, true}) {
			const oc::ObservationPreparationOptions options{.debias = debias, .weightProfile = profile};
			auto fresh = raw;
			oc::prepareObservations(fresh, options);
			auto actual = oc::prepareObservationCopy(previous, options);
			for (const auto& [a, b] : {std::pair{&actual.all, &fresh.all}, std::pair{&actual.iod, &fresh.iod}}) {
				ASSERT_EQ(a->size(), b->size());
				for (std::size_t i = 0; i < a->size(); ++i) {
					EXPECT_EQ((*a)[i].id, (*b)[i].id);
					EXPECT_EQ((*a)[i].RA, (*b)[i].RA);
					EXPECT_EQ((*a)[i].DEC, (*b)[i].DEC);
					EXPECT_EQ((*a)[i].vector.x, (*b)[i].vector.x);
					EXPECT_EQ((*a)[i].vector.y, (*b)[i].vector.y);
					EXPECT_EQ((*a)[i].vector.z, (*b)[i].vector.z);
					EXPECT_EQ((*a)[i].biasRaArcsec, (*b)[i].biasRaArcsec);
					EXPECT_EQ((*a)[i].biasDecArcsec, (*b)[i].biasDecArcsec);
				}
			}
			previous = std::move(actual);
		}
	}
	EXPECT_FALSE(raw.preparation);
}

TEST_F(FitSessionTest, PipelineHonorsManualSelectionFromIodThroughFinalFit)
{
	auto data = fitContractCase();
	data.config.sigmaClip.enabled = true;
	const auto excluded = data.observations.all[1].id;
	const std::vector<oc::ObservationOverride> overrides{
	    {excluded, oc::ObservationDisposition::Exclude},
	    {data.observations.all[2].id, oc::ObservationDisposition::Include}};
	const auto result = oc::FitPipeline{data.config}.run({data.observations, {}, overrides});
	ASSERT_TRUE(result.fit.orbit);
	EXPECT_EQ(result.fitInputObservations, data.observations.all.size() - 1);
	EXPECT_EQ(result.observationOverrides, overrides);
	EXPECT_EQ(result.observations.all.size(), data.observations.all.size());
	EXPECT_EQ(result.observations.all[1].id, excluded);
	EXPECT_EQ(result.fit.allResiduals.count(), data.observations.all.size());
	EXPECT_EQ(std::ranges::count(result.fit.keptObsIndices, 1u), 0);
	EXPECT_EQ(std::ranges::count(result.fit.keptObsIndices, 2u), 1);
	auto filtered = data.observations;
	std::erase_if(filtered.all, [&](const auto& o) { return o.id == excluded; });
	std::erase_if(filtered.iod, [&](const auto& o) { return o.id == excluded; });
	const auto reference = oc::FitPipeline{data.config}.run({filtered, {}});
	ASSERT_TRUE(reference.initial.orbit);
	ASSERT_TRUE(result.initial.orbit);
	EXPECT_EQ(result.initial.method, reference.initial.method);
	EXPECT_EQ(result.initial.orbit->getEpoch(), reference.initial.orbit->getEpoch());
	EXPECT_EQ(result.initial.orbit->getR().x, reference.initial.orbit->getR().x);
	EXPECT_EQ(result.initial.orbit->getR().y, reference.initial.orbit->getR().y);
	EXPECT_EQ(result.initial.orbit->getR().z, reference.initial.orbit->getR().z);
	EXPECT_EQ(result.initial.orbit->getV().x, reference.initial.orbit->getV().x);
	EXPECT_EQ(result.initial.orbit->getV().y, reference.initial.orbit->getV().y);
	EXPECT_EQ(result.initial.orbit->getV().z, reference.initial.orbit->getV().z);

	const std::vector<oc::ObservationOverride> invalid{{"missing", oc::ObservationDisposition::Exclude}};
	bool startedIod = false;
	EXPECT_THROW((void)oc::FitPipeline{data.config}.run({data.observations, {}, invalid},
	                                                    {.progress = [&](const auto&) { startedIod = true; }}),
	             std::invalid_argument);
	EXPECT_FALSE(startedIod);
}

TEST_F(FitSessionTest, PipelineResultOwnsPreparedInputAfterCallerReleasesIt)
{
	oc::FitPipelineResult result;
	oc::ObservationSet reference;
	{
		auto data = fitContractCase();
		reference = data.observations;
		result = oc::FitPipeline{data.config}.run({data.observations, {}});
		ASSERT_TRUE(result.fit.orbit);
		EXPECT_NE(result.observations.all.data(), data.observations.all.data());
		EXPECT_NE(result.observations.iod.data(), data.observations.iod.data());
		EXPECT_EQ(data.observations.all.size(), reference.all.size());
		EXPECT_EQ(data.observations.iod.size(), reference.iod.size());
		data.observations.all.clear();
		data.observations.iod.clear();
	}
	const auto& actual = result.observations;
	EXPECT_EQ(actual.objectName, reference.objectName);
	EXPECT_EQ(actual.source, reference.source);
	EXPECT_EQ(actual.format, reference.format);
	ASSERT_TRUE(actual.preparation);
	ASSERT_TRUE(reference.preparation);
	EXPECT_EQ(actual.preparation->debiasRequested, reference.preparation->debiasRequested);
	EXPECT_EQ(actual.preparation->weightProfile, reference.preparation->weightProfile);
	EXPECT_EQ(actual.preparation->biasTable, reference.preparation->biasTable);
	EXPECT_EQ(actual.preparation->correctedAll, reference.preparation->correctedAll);
	EXPECT_EQ(actual.preparation->correctedIod, reference.preparation->correctedIod);
	const auto check = [](const auto& rows, const auto& expected) {
		ASSERT_EQ(rows.size(), expected.size());
		for (std::size_t i = 0; i < rows.size(); ++i) {
			EXPECT_EQ(rows[i].JD, expected[i].JD);
			EXPECT_EQ(rows[i].RA, expected[i].RA);
			EXPECT_EQ(rows[i].DEC, expected[i].DEC);
			EXPECT_EQ(rows[i].obsCode, expected[i].obsCode);
			EXPECT_EQ(rows[i].catalogCode, expected[i].catalogCode);
			EXPECT_EQ(rows[i].astCat, expected[i].astCat);
			EXPECT_EQ(rows[i].submittedRaRad(), expected[i].submittedRaRad());
			EXPECT_EQ(rows[i].submittedDecRad(), expected[i].submittedDecRad());
			EXPECT_EQ(rows[i].biasRaArcsec, expected[i].biasRaArcsec);
			EXPECT_EQ(rows[i].biasDecArcsec, expected[i].biasDecArcsec);
			EXPECT_EQ(rows[i].rmsRaArcsec, expected[i].rmsRaArcsec);
			EXPECT_EQ(rows[i].rmsDecArcsec, expected[i].rmsDecArcsec);
		}
	};
	check(actual.all, reference.all);
	check(actual.iod, reference.iod);
	const auto fitted = result.fittedObservations();
	ASSERT_EQ(fitted.size(), result.fit.residuals.count());
	for (std::size_t i = 0; i < fitted.size(); ++i)
		EXPECT_EQ(&fitted[i], &actual.all.at(result.fit.keptObsIndices.at(i)));
}

TEST_F(FitSessionTest, PipelinePreparesOneFullPoolWithManualExclusions)
{
	auto data = fitContractCase();
	data.config.sigmaClip.enabled = true;
	OrbitOptimizerConfig weights;
	weights.useWeights = data.config.useWeights;
	weights.weightProfile = data.config.weightProfile;
	weights.obsWeighting = data.config.obsWeighting;
	const auto prepared = oc::fit::prepareFitObservations(data.observations.all, weights);
	const std::array overrides{
	    oc::ObservationOverride{data.observations.all[1].id, oc::ObservationDisposition::Exclude}};
	int preparations = 0;
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		preparations += message.starts_with("ObservationPipeline: ready");
	});
	const auto result = oc::FitPipeline{data.config}.run({data.observations, {}, overrides});
	ASSERT_TRUE(result.fit.orbit);
	EXPECT_EQ(preparations, 1);
	ASSERT_EQ(result.fit.allResiduals.count(), prepared.size());
	for (std::size_t i = 0; i < prepared.size(); ++i) {
		const auto& weight = result.fit.allResiduals.data()[i].astrometricWeight;
		ASSERT_TRUE(weight);
		EXPECT_NEAR(weight->sigmaRaArcsec, prepared[i].sigmaRa / Constants::ARCSEC_RAD, 1e-14);
		EXPECT_NEAR(weight->sigmaDecArcsec, prepared[i].sigmaDec / Constants::ARCSEC_RAD, 1e-14);
	}
}

TEST_F(FitSessionTest, OutputEpochBelongsToRequestAndRetainsFitStateAndCovariance)
{
	auto data = fitContractCase();
	data.config.sigmaClip.enabled = true;
	const oc::FitSession session{data.config};
	for (bool evaluateOnly : {false, true}) {
		oc::FitRequest request{.initialOrbit = data.initial,
		                       .observations = data.observations.all,
		                       .options = {.evaluateOnly = evaluateOnly}};
		const auto reference = session.optimize(request);
		ASSERT_TRUE(reference.orbit);
		ASSERT_TRUE(reference.fitEpochOrbit);
		ASSERT_TRUE(reference.uncertainty.hasStateCovariance());
		EXPECT_DOUBLE_EQ(reference.orbit->getEpoch(), reference.fitEpochOrbit->getEpoch());

		request.outputEpochJdTdb = reference.fitEpochOrbit->getEpoch() + 5.0;
		const auto transported = session.optimize(request);
		ASSERT_TRUE(transported.orbit);
		ASSERT_TRUE(transported.fitEpochOrbit);
		EXPECT_DOUBLE_EQ(transported.orbit->getEpoch(), *request.outputEpochJdTdb);
		EXPECT_NE(transported.orbit->getR().x, transported.fitEpochOrbit->getR().x);

		request.outputEpochJdTdb.reset();
		const auto repeated = session.optimize(request);
		ASSERT_TRUE(repeated.orbit);
		ASSERT_TRUE(repeated.fitEpochOrbit);
		EXPECT_DOUBLE_EQ(repeated.orbit->getEpoch(), reference.orbit->getEpoch());
		for (const auto* result : {&transported, &repeated}) {
			EXPECT_DOUBLE_EQ(result->fitEpochOrbit->getEpoch(), reference.fitEpochOrbit->getEpoch());
			EXPECT_DOUBLE_EQ(result->fitEpochOrbit->getR().x, reference.fitEpochOrbit->getR().x);
			EXPECT_DOUBLE_EQ(result->fitEpochOrbit->getR().y, reference.fitEpochOrbit->getR().y);
			EXPECT_DOUBLE_EQ(result->fitEpochOrbit->getR().z, reference.fitEpochOrbit->getR().z);
			EXPECT_DOUBLE_EQ(result->fitEpochOrbit->getV().x, reference.fitEpochOrbit->getV().x);
			EXPECT_DOUBLE_EQ(result->fitEpochOrbit->getV().y, reference.fitEpochOrbit->getV().y);
			EXPECT_DOUBLE_EQ(result->fitEpochOrbit->getV().z, reference.fitEpochOrbit->getV().z);
			EXPECT_TRUE((result->uncertainty.cartesianStateCovariance.array() ==
			             reference.uncertainty.cartesianStateCovariance.array())
			                .all());
			EXPECT_TRUE((result->uncertainty.fitParameterCovarianceFormal.array() ==
			             reference.uncertainty.fitParameterCovarianceFormal.array())
			                .all());
			EXPECT_DOUBLE_EQ(result->cost, reference.cost);
			EXPECT_DOUBLE_EQ(result->normalizedRms, reference.normalizedRms);
			EXPECT_EQ(result->keptObsIndices, reference.keptObsIndices);
			EXPECT_EQ(result->sigmaClipPassesRun, reference.sigmaClipPassesRun);
		}
	}
}

TEST_F(FitSessionTest, EvaluateOnlyFreezesConfiguredSelectionAndNonGrav)
{
	auto data = fitContractCase();
	data.observations.all.back().RA += 0.2;
	std::vector<std::uint8_t> mask(data.observations.all.size(), 1);
	mask[1] = 0;
	data.config.sigmaClip.enabled = true;
	data.config.nonGrav.enabled = true;
	data.config.nonGrav.model = oc::NonGravModelKind::Yarkovsky;
	data.config.nonGrav.useInitialValues = true;
	data.config.nonGrav.initialA2 = -3e-14;
	oc::FitRequest request{.initialOrbit = data.initial,
	                       .observations = data.observations.all,
	                       .initialIncluded = mask,
	                       .options = {.selection = oc::FitSelectionMode::Fixed,
	                                   .uncertainty = oc::FitUncertaintyMode::Skip,
	                                   .evaluateOnly = true}};
	int lmStarts = 0;
	int bootstrapStarts = 0;
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		lmStarts += message.starts_with(solverStart);
		bootstrapStarts += message.starts_with("=== Non-grav bootstrap");
	});
	const oc::FitSession session{data.config};
	const auto fixed = session.optimize(request);
	request.options.selection = oc::FitSelectionMode::Configured;
	const auto evaluated = session.optimize(request);
	ASSERT_TRUE(evaluated.fitEpochOrbit);
	ASSERT_TRUE(fixed.fitEpochOrbit);
	EXPECT_EQ(evaluated.keptObsIndices, fixed.keptObsIndices);
	EXPECT_EQ(evaluated.residuals.count(), data.observations.all.size() - 1);
	EXPECT_EQ(evaluated.allResiduals.count(), data.observations.all.size());
	EXPECT_EQ(evaluated.sigmaClipPassesRun, 0);
	EXPECT_EQ(evaluated.sigmaClipRejected, 1);
	EXPECT_EQ(evaluated.iterations, 0);
	EXPECT_FALSE(evaluated.converged);
	EXPECT_FALSE(evaluated.cancelled);
	EXPECT_DOUBLE_EQ(evaluated.cost, fixed.cost);
	EXPECT_DOUBLE_EQ(evaluated.fitEpochOrbit->getR().x, fixed.fitEpochOrbit->getR().x);
	EXPECT_DOUBLE_EQ(evaluated.fitEpochOrbit->getV().y, fixed.fitEpochOrbit->getV().y);
	EXPECT_DOUBLE_EQ(evaluated.finalNonGrav.A2(), data.config.nonGrav.initialA2);
	EXPECT_EQ(lmStarts, 0);
	EXPECT_EQ(bootstrapStarts, 0);
}

TEST_F(FitSessionTest, OptimizerCostDoesNotDependOnResidualReport)
{
	const auto data = fitContractCase();
	for (bool evaluateOnly : {false, true}) {
		OrbitOptimizerConfig config;
		config.nbody = data.config.nbody;
		config.useWeights = data.config.useWeights;
		config.weightProfile = data.config.weightProfile;
		config.obsWeighting = data.config.obsWeighting;
		config.computeUncertainty = false;
		config.logDiagnostics = false;
		config.logCacheStats = false;
		config.evaluateOnly = evaluateOnly;
		const auto prepared = oc::fit::prepareFitObservations(data.observations.all, config);
		OrbitOptimizer withRows(data.initial, prepared, config);
		const auto reference = withRows.optimize();
		config.populateResiduals = false;
		OrbitOptimizer withoutRows(data.initial, prepared, config);
		const auto result = withoutRows.optimize();
		ASSERT_GT(reference.cost, 0.0);
		ASSERT_GT(reference.normalizedRms, 0.0);
		EXPECT_EQ(result.residuals.count(), 0u);
		EXPECT_DOUBLE_EQ(result.cost, reference.cost);
		EXPECT_DOUBLE_EQ(result.normalizedRms, reference.normalizedRms);
		EXPECT_EQ(result.iterations, reference.iterations);
		EXPECT_DOUBLE_EQ(result.fitEpochOrbit->getR().x, reference.fitEpochOrbit->getR().x);
	}
}

TEST_F(FitSessionTest, NormalizedRmsIncludesCorrelatedOpticalAndEachRadarChannel)
{
	const auto data = fitContractCase();
	OrbitOptimizerConfig config;
	config.nbody = data.config.nbody;
	config.useWeights = true;
	config.computeUncertainty = false;
	config.evaluateOnly = true;
	config.logDiagnostics = false;
	config.logCacheStats = false;
	auto prepared = oc::fit::prepareFitObservations(data.observations.all, config);
	ASSERT_GE(prepared.size(), 4u);
	prepared.resize(4);
	for (auto& observation : prepared) {
		observation.sigmaRa = 0.5 * Constants::ARCSEC_RAD;
		observation.sigmaDec = 0.8 * Constants::ARCSEC_RAD;
		observation.sigmaCorr = 0.6;
		observation.timeSigmaDays = 0.0;
		observation.raObs += 0.7 * Constants::ARCSEC_RAD;
		observation.decObs -= 0.4 * Constants::ARCSEC_RAD;
	}
	for (bool delay : {false, true}) {
		auto radar = prepared.front();
		radar.transmitter = radar.observer;
		radar.radar = oc::obs::RadarObservation{.dopplerHz = oc::obs::RadarMeasurement{0.0, 200.0},
		                                        .frequencyHz = 1e9,
		                                        .transmitter = radar.observer.obsCode};
		if (delay)
			radar.radar->delaySeconds = oc::obs::RadarMeasurement{1.0, 0.25};
		prepared.push_back(radar);
	}
	OrbitOptimizer optimizer(data.initial, prepared, config);
	const auto result = optimizer.optimize();
	ASSERT_EQ(result.residuals.count(), 6u);
	double chiSquare = 0.0;
	int dimensions = 0;
	for (const auto& residual : result.residuals.data()) {
		if (residual.radar) {
			for (const auto& channel : {residual.radar->delaySeconds, residual.radar->dopplerHz}) {
				if (!channel)
					continue;
				const double r = channel->residual / channel->sigma;
				chiSquare += r * r;
				++dimensions;
			}
		} else {
			const double ra = residual.dRA_arcsec / 0.5;
			const double dec = residual.dDEC_arcsec / 0.8;
			chiSquare += (ra * ra - 2.0 * 0.6 * ra * dec + dec * dec) / (1.0 - 0.6 * 0.6);
			dimensions += 2;
		}
	}
	ASSERT_EQ(dimensions, 11);
	const double expected = std::sqrt(chiSquare / dimensions);
	EXPECT_NEAR(result.normalizedRms, expected, expected * 1e-12);

	std::vector<std::uint8_t> included(prepared.size(), 1);
	included.back() = 0;
	const auto masked = optimizer.optimizeWithSelection(included, {});
	ASSERT_EQ(masked.residuals.count(), 5u);
	ASSERT_EQ(masked.allResiduals.count(), 6u);
	const auto& excludedRadar = result.residuals.data().back().radar;
	ASSERT_TRUE(excludedRadar);
	for (const auto& channel : {excludedRadar->delaySeconds, excludedRadar->dopplerHz}) {
		if (channel) {
			chiSquare -= std::pow(channel->residual / channel->sigma, 2);
			--dimensions;
		}
	}
	EXPECT_NEAR(masked.normalizedRms, std::sqrt(chiSquare / dimensions), 1e-12);
	const auto reset = optimizer.optimize();
	EXPECT_DOUBLE_EQ(reset.normalizedRms, result.normalizedRms);

	config.computeUncertainty = true;
	OrbitOptimizer covarianceOptimizer(data.initial, prepared, config);
	const auto covarianceFit = covarianceOptimizer.optimizeWithSelection(included, {});
	ASSERT_TRUE(covarianceFit.fitEpochOrbit);
	const auto r = data.initial->getR();
	const auto v = data.initial->getV() * Constants::GAUSS_K;
	const std::array physical{r.x, r.y, r.z, v.x, v.y, v.z};
	const auto linear = covarianceOptimizer.weightedObservationLinearizationPhysical(physical);
	ASSERT_TRUE(linear.valid);
	Eigen::MatrixXd activeJ = linear.jacobian.topRows(dimensions);
	for (int i = 3; i < 6; ++i)
		activeJ.col(i) *= Constants::GAUSS_K;
	const Eigen::JacobiSVD<Eigen::MatrixXd> svd(activeJ, Eigen::ComputeThinV);
	ASSERT_GT(svd.singularValues().minCoeff(), 0.0);
	Eigen::MatrixXd expectedCovariance = svd.matrixV() *
	                                     svd.singularValues().array().square().inverse().matrix().asDiagonal() *
	                                     svd.matrixV().transpose();
	for (int i = 3; i < 6; ++i) {
		expectedCovariance.row(i) *= Constants::GAUSS_K;
		expectedCovariance.col(i) *= Constants::GAUSS_K;
	}
	EXPECT_TRUE(covarianceFit.uncertainty.fitParameterCovarianceFormal.isApprox(expectedCovariance, 1e-6))
	    << "actual:\n"
	    << covarianceFit.uncertainty.fitParameterCovarianceFormal << "\nexpected:\n"
	    << expectedCovariance;

	config.populateResiduals = false;
	config.computeUncertainty = false;
	OrbitOptimizer withoutRows(data.initial, prepared, config);
	EXPECT_DOUBLE_EQ(withoutRows.optimize().normalizedRms, result.normalizedRms);

	prepared.resize(4);
	config.useWeights = false;
	OrbitOptimizer unweighted(data.initial, prepared, config);
	EXPECT_TRUE(std::isnan(unweighted.optimize().normalizedRms));
}

TEST_F(FitSessionTest, RankDeficientFitDoesNotPublishCovariance)
{
	const auto data = fitContractCase();
	OrbitOptimizerConfig config;
	config.nbody = data.config.nbody;
	config.useWeights = true;
	config.evaluateOnly = true;
	config.computeUncertainty = true;
	config.logDiagnostics = false;
	config.logCacheStats = false;
	const auto prepared = oc::fit::prepareFitObservations(data.observations.all, config);
	ASSERT_FALSE(prepared.empty());
	std::vector<oc::obs::PreparedObservation> repeated(8, prepared.front());
	OrbitOptimizer optimizer(data.initial, repeated, config);
	const auto result = optimizer.optimize();
	ASSERT_TRUE(result.orbit);
	EXPECT_EQ(result.uncertainty.fitParameterCovariance.size(), 0);
	EXPECT_FALSE(result.uncertainty.hasElementCovariance);
}

TEST_F(FitSessionTest, GaiaFprMergePreservesProvenanceAndModelSettings)
{
	const auto data = fitContractCase();
	const oc::GaiaFprObservation observation{59000.5, 120.0, 20.0, 0.002, 0.003, 0.4, Vec3D(0.01, 0.002, 0.001)};
	const auto merged = oc::withGaiaFprObservations(data.observations, std::span(&observation, 1), 1.5);
	ASSERT_EQ(merged.all.size(), data.observations.all.size() + 1);
	ASSERT_EQ(merged.iod.size(), data.observations.iod.size() + 1);
	const auto it = std::ranges::find_if(merged.all, [](const auto& o) { return o.obsCode == "258"; });
	ASSERT_NE(it, merged.all.end());
	EXPECT_EQ(it->astrometricConvention, oc::obs::AstrometricConvention::GaiaFpr);
	EXPECT_TRUE(it->isSatellite);
	EXPECT_TRUE(it->hasObserverSigma);
	OrbitOptimizerConfig config;
	config.useWeights = true;
	config.gaiaPhotocenterRadiusKm = merged.gaiaPhotocenterRadiusKm;
	const auto prepared = oc::fit::prepareFitObservations(merged.all, config);
	const auto p = std::ranges::find_if(prepared, [](const auto& o) { return o.observer.obsCode == "258"; });
	ASSERT_NE(p, prepared.end());
	EXPECT_TRUE(p->observerSigma);
	EXPECT_DOUBLE_EQ(p->photocenterRadiusAu, 1.5 / (Constants::AU_IN_METERS / 1000.0));
	EXPECT_THROW((void)oc::withGaiaFprObservations(merged, std::span(&observation, 1), 1.5), std::invalid_argument);
}

TEST_F(FitSessionTest, GaiaFprMergeChecksExistingGaiaCoverage)
{
	auto data = fitContractCase();
	const oc::GaiaFprObservation fpr{59000.5, 120.0, 20.0, 0.002, 0.003, 0.4, Vec3D(0.01, 0.002, 0.001)};
	auto existing = data.observations.all.front();
	existing.obsCode = "258";
	existing.astrometricConvention = oc::obs::AstrometricConvention::CatalogReferenced;
	existing.hasSplitUtc = false;
	existing.hasTimeAccuracy = false;
	existing.timeScale = oc::time::Scale::UTC;
	existing.JD = 2400000.5 + 58998.0;
	data.observations.all.push_back(existing);
	EXPECT_NO_THROW((void)oc::withGaiaFprObservations(data.observations, std::span(&fpr, 1)));

	data.observations.all.back().JD = 2400000.5 + 59000.5;
	EXPECT_THROW((void)oc::withGaiaFprObservations(data.observations, std::span(&fpr, 1)), std::invalid_argument);
	data.observations.all.back().JD = 2400000.5 + 58999.6;
	EXPECT_THROW((void)oc::withGaiaFprObservations(data.observations, std::span(&fpr, 1)), std::invalid_argument);
	data.observations.all.back().JD = 2400000.5 + 58998.0;
	data.observations.all.back().timeScale = oc::time::Scale::TDB;
	EXPECT_THROW((void)oc::withGaiaFprObservations(data.observations, std::span(&fpr, 1)), std::invalid_argument);
}

TEST_F(FitSessionTest, ClippingReusesPreparedRowsAndPreservesOriginalIndices)
{
	auto data = fitContractCase();
	const auto outlier = data.observations.all.size() / 2;
	data.observations.all[outlier].RA += 20.0 * Constants::ARCSEC_RAD;
	std::vector<std::uint8_t> mask(data.observations.all.size(), 1);
	mask[1] = 0;
	data.config.sigmaClip.enabled = true;
	int preparations = 0;
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		preparations += message.starts_with("ObservationPipeline: ready");
	});
	const auto result = oc::FitSession{data.config}.optimize(
	    {.initialOrbit = data.initial, .observations = data.observations.all, .initialIncluded = mask});
	ASSERT_TRUE(result.fitEpochOrbit);
	EXPECT_TRUE(result.uncertainty.hasStateCovariance());
	EXPECT_GT(result.sigmaClipPassesRun, 0);
	EXPECT_EQ(preparations, 1);
	EXPECT_TRUE(result.converged);
	EXPECT_NE(std::ranges::find(result.keptObsIndices, 1), result.keptObsIndices.end());
	EXPECT_EQ(std::ranges::find(result.keptObsIndices, outlier), result.keptObsIndices.end());
	EXPECT_TRUE(std::ranges::is_sorted(result.keptObsIndices));
	EXPECT_EQ(result.residuals.count(), result.keptObsIndices.size());
	EXPECT_EQ(result.allResiduals.count(), data.observations.all.size());
	EXPECT_EQ(result.sigmaClipRejected, data.observations.all.size() - result.keptObsIndices.size());
	EXPECT_EQ(mask[1], 0);
	EXPECT_EQ(mask[outlier], 1);
	for (const auto& residual : result.allResiduals.data()) {
		ASSERT_TRUE(residual.astrometricWeight);
		EXPECT_TRUE(residual.astrometricWeight->postFitChiSquared);
	}
	double chiSquare = 0.0;
	for (const auto& residual : result.residuals.data()) {
		ASSERT_TRUE(residual.astrometricWeight);
		chiSquare += residual.astrometricWeight->chiSquared;
	}
	EXPECT_NEAR(result.normalizedRms, std::sqrt(chiSquare / (2.0 * result.keptObsIndices.size())), 1e-12);
}

TEST_F(FitSessionTest, SelectionModeBelongsToRequest)
{
	auto data = fitContractCase();
	std::vector<std::uint8_t> mask(data.observations.all.size(), 1);
	mask[1] = 0;
	data.config.sigmaClip.enabled = true;
	oc::FitRequest request{.initialOrbit = data.initial,
	                       .observations = data.observations.all,
	                       .initialIncluded = mask,
	                       .options = {.selection = oc::FitSelectionMode::Fixed}};
	const oc::FitSession session{data.config};
	const auto reference = session.optimize(request);
	request.options.selection = oc::FitSelectionMode::Configured;
	const auto refined = session.optimize(request);
	ASSERT_TRUE(refined.fitEpochOrbit);
	EXPECT_EQ(refined.keptObsIndices.size(), data.observations.all.size());
	EXPECT_GT(refined.sigmaClipPassesRun, 0);
	request.options.selection = oc::FitSelectionMode::Fixed;
	const auto result = session.optimize(request);
	ASSERT_TRUE(reference.fitEpochOrbit);
	ASSERT_TRUE(result.fitEpochOrbit);
	EXPECT_DOUBLE_EQ(result.cost, reference.cost);
	EXPECT_EQ(result.iterations, reference.iterations);
	EXPECT_EQ(result.sigmaClipPassesRun, 0);
	EXPECT_EQ(result.sigmaClipRejected, 1);
	EXPECT_EQ(result.keptObsIndices, reference.keptObsIndices);
	EXPECT_EQ(result.allResiduals.count(), data.observations.all.size());
	EXPECT_DOUBLE_EQ(result.fitEpochOrbit->getEpoch(), reference.fitEpochOrbit->getEpoch());
	EXPECT_DOUBLE_EQ(result.fitEpochOrbit->getR().x, reference.fitEpochOrbit->getR().x);
	EXPECT_DOUBLE_EQ(result.fitEpochOrbit->getV().z, reference.fitEpochOrbit->getV().z);
	EXPECT_TRUE((result.uncertainty.fitParameterCovarianceFormal.array() ==
	             reference.uncertainty.fitParameterCovarianceFormal.array())
	                .all());
}

TEST_F(FitSessionTest, SelectionStartsAfterConvergedWeightedAndUnweightedFits)
{
	for (const auto [useNBody, weighted] : {std::pair{false, false}, {false, true}, {true, false}, {true, true}}) {
		SCOPED_TRACE(::testing::Message() << "nbody=" << useNBody << " weighted=" << weighted);
		auto data = fitContractCase();
		data.config.nbody.enabled = useNBody;
		data.config.useWeights = weighted;
		data.config.sigmaClip.enabled = true;
		int checks = 0;
		bool converged = false;
		const FitLogCallback log([&](LogLevel, const std::string& message) {
			if (message.starts_with(solverFinish))
				converged = converged || message.find("converged=1") != std::string::npos;
			if (message.starts_with("[HystereticRejector]")) {
				++checks;
				EXPECT_TRUE(converged);
			}
		});
		const auto result = oc::FitSession{data.config}.optimize({data.initial, data.observations.all});
		ASSERT_TRUE(result.orbit);
		if (!result.converged) {
			EXPECT_FALSE(useNBody);
			EXPECT_TRUE(weighted);
			EXPECT_EQ(result.selectionStatus, oc::FitSelectionStatus::PassLimit);
			EXPECT_EQ(result.sigmaClipPassesRun, data.config.sigmaClip.maxPasses);
			EXPECT_GT(checks, 0);
			continue;
		}
		EXPECT_TRUE(result.converged);
		EXPECT_GT(checks, 0);
	}
}

TEST_F(FitSessionTest, FinalizedSelectionDiagnosticsMatchFreshEvaluation)
{
	for (bool nonGrav : {false, true}) {
		SCOPED_TRACE(nonGrav);
		auto data = fitContractCase();
		data.config.sigmaClip.enabled = true;
		data.config.nonGrav.enabled = nonGrav;
		data.config.nonGrav.model = oc::NonGravModelKind::Yarkovsky;
		data.observations.all[2].RA += 20.0 * Constants::ARCSEC_RAD;
		const std::array overrides{
		    oc::ObservationOverride{data.observations.all[1].id, oc::ObservationDisposition::Exclude}};
		int checks = 0;
		const FitLogCallback log([&](LogLevel, const std::string& message) {
			checks += message.starts_with("[HystereticRejector] rejectChi2=");
		});
		const auto result = oc::FitSession{data.config}.optimize({.initialOrbit = data.initial,
		                                                          .observations = data.observations.all,
		                                                          .outputEpochJdTdb = data.initial->getEpoch() + 30.0,
		                                                          .observationOverrides = overrides});
		ASSERT_TRUE(result.converged);
		ASSERT_TRUE(result.fitEpochOrbit);
		ASSERT_EQ(result.allResiduals.count(), data.observations.all.size());
		EXPECT_GT(checks, 0);
		EXPECT_NE(result.orbit->getEpoch(), result.fitEpochOrbit->getEpoch());

		OrbitOptimizerConfig config;
		config.nbody = data.config.nbody;
		config.useWeights = data.config.useWeights;
		config.weightProfile = data.config.weightProfile;
		config.obsWeighting = data.config.obsWeighting;
		config.evaluateOnly = true;
		config.computeUncertainty = false;
		config.nonGrav = data.config.nonGrav;
		config.nonGrav.useInitialValues = true;
		config.nonGrav.initialA2 = result.finalNonGrav.A2();
		const auto prepared = oc::fit::prepareFitObservations(data.observations.all, config);
		OrbitOptimizer verifier(result.fitEpochOrbit, prepared, config);
		const auto evaluated = verifier.optimize();
		ASSERT_TRUE(evaluated.fitEpochOrbit);
		std::vector<std::uint8_t> mask(prepared.size(), 0);
		for (const auto index : result.keptObsIndices)
			mask[index] = 1;
		const auto fresh = verifier.postFitDiagnostics(result.uncertainty.fitParameterCovarianceFormal, mask);
		ASSERT_TRUE(fresh.valid);
		ASSERT_EQ(fresh.residuals.count(), result.allResiduals.count());
		for (std::size_t i = 0; i < prepared.size(); ++i) {
			const auto& actual = result.allResiduals.data()[i];
			const auto& expected = fresh.residuals.data()[i];
			EXPECT_DOUBLE_EQ(actual.dRA_arcsec, expected.dRA_arcsec);
			EXPECT_DOUBLE_EQ(actual.dDEC_arcsec, expected.dDEC_arcsec);
			ASSERT_TRUE(actual.astrometricWeight);
			ASSERT_TRUE(expected.astrometricWeight);
			EXPECT_EQ(actual.astrometricWeight->postFitChiSquared, expected.astrometricWeight->postFitChiSquared);
		}
	}
}

TEST_F(FitSessionTest, ExecutionGraphMatchesOptionCombinations)
{
	auto data = fitContractCase();
	data.observations.all[2].RA += 20.0 * Constants::ARCSEC_RAD;
	std::vector<std::uint8_t> mask(data.observations.all.size(), 1);
	mask[0] = 0;
	const std::array overrides{
	    oc::ObservationOverride{data.observations.all[0].id, oc::ObservationDisposition::Include},
	    oc::ObservationOverride{data.observations.all[1].id, oc::ObservationDisposition::Exclude}};
	const double outputEpoch = data.initial->getEpoch() + 30.0;
	for (const auto model :
	     {oc::NonGravModelKind::None, oc::NonGravModelKind::Yarkovsky, oc::NonGravModelKind::MarsdenH2O}) {
		for (bool clipping : {false, true}) {
			for (int mode : {0, 1, 2}) {
				for (bool warm : {false, true}) {
					for (const auto uncertainty : {oc::FitUncertaintyMode::Skip, oc::FitUncertaintyMode::Compute}) {
						SCOPED_TRACE(::testing::Message()
						             << "model=" << static_cast<int>(model) << " clip=" << clipping << " mode=" << mode
						             << " warm=" << warm << " covariance=" << static_cast<int>(uncertainty));
						auto cfg = data.config;
						cfg.sigmaClip.enabled = clipping;
						cfg.nonGrav.enabled = model != oc::NonGravModelKind::None;
						cfg.nonGrav.model = model;
						NonGravModel ng;
						ng.setEnabled(cfg.nonGrav.enabled);
						ng.setModel(static_cast<NgModel>(model));
						ng.setParams(1e-13, -3e-14, 1e-14);
						const bool evaluate = mode == 2;
						const bool selects = clipping && mode == 0;
						int preparations = 0, lmStarts = 0, bootstraps = 0;
						int transports = 0, selectionChecks = 0;
						std::string covarianceWarning;
						const FitLogCallback log([&](LogLevel, const std::string& message) {
							if (message.starts_with("[Covariance] Jacobian"))
								covarianceWarning = message;
							preparations += message.starts_with("ObservationPipeline: ready");
							lmStarts += message.starts_with(solverStart);
							bootstraps += message.starts_with("=== Non-grav bootstrap");
							transports += message.starts_with("==== ORBIT AT TARGET EPOCH");
							selectionChecks += message.starts_with("[HystereticRejector] rejectChi2=");
						});
						const auto result = oc::FitSession{cfg}.optimize(
						    {.initialOrbit = data.initial,
						     .observations = data.observations.all,
						     .initialIncluded = mask,
						     .warmStartNonGrav = warm ? &ng : nullptr,
						     .options = {.selection =
						                     mode == 1 ? oc::FitSelectionMode::Fixed : oc::FitSelectionMode::Configured,
						                 .uncertainty = uncertainty,
						                 .evaluateOnly = evaluate},
						     .outputEpochJdTdb = outputEpoch,
						     .observationOverrides = overrides});
						ASSERT_TRUE(result.orbit);
						EXPECT_DOUBLE_EQ(result.orbit->getEpoch(), outputEpoch);
						EXPECT_EQ(preparations, 1);
						EXPECT_EQ(transports, 1);
						EXPECT_EQ(lmStarts > 0, !evaluate);
						EXPECT_EQ(selectionChecks > 0, selects);
						const bool needsBootstrap = cfg.nonGrav.enabled && !warm && !evaluate;
						EXPECT_EQ(bootstraps, needsBootstrap ? 1 : 0);
						// With zero Marsden accelerations, DT has no observable effect.
						const bool unobservableDt = evaluate && !warm && model == oc::NonGravModelKind::MarsdenH2O;
						EXPECT_EQ(result.uncertainty.hasStateCovariance(),
						          uncertainty == oc::FitUncertaintyMode::Compute && !unobservableDt)
						    << covarianceWarning;
						EXPECT_NE(std::ranges::find(result.keptObsIndices, 0), result.keptObsIndices.end());
						EXPECT_EQ(std::ranges::find(result.keptObsIndices, 1), result.keptObsIndices.end());
						EXPECT_EQ(result.allResiduals.count(), data.observations.all.size());
						if (warm && model == oc::NonGravModelKind::MarsdenH2O && mode == 1) {
							EXPECT_TRUE(result.converged);
							EXPECT_TRUE(std::isfinite(result.finalNonGrav.DT()));
							if (uncertainty == oc::FitUncertaintyMode::Compute)
								EXPECT_TRUE(result.uncertainty.fitParameterCovarianceFormal.allFinite());
						}
						if (evaluate) {
							EXPECT_EQ(lmStarts, 0);
							EXPECT_FALSE(result.converged);
						}
					}
				}
			}
		}
	}
}

TEST_F(FitSessionTest, SelectionPassLimitChecksTheLastCorrection)
{
	auto data = fitContractCase();
	data.observations.all[data.observations.all.size() / 2].RA += 20 * Constants::ARCSEC_RAD;
	data.config.sigmaClip.enabled = true;
	data.config.sigmaClip.maxPasses = 1;
	const auto result = oc::FitSession{data.config}.optimize({data.initial, data.observations.all});
	ASSERT_TRUE(result.orbit);
	EXPECT_EQ(result.sigmaClipPassesRun, 1);
	EXPECT_TRUE(result.converged);
	EXPECT_EQ(result.selectionStatus, oc::FitSelectionStatus::Stable);
	EXPECT_TRUE(result.uncertainty.hasStateCovariance());
}

TEST_F(FitSessionTest, SelectionPassLimitPreservesTheLastFittedMaskWhenUnstable)
{
	auto data = fitContractCase();
	data.observations.all[data.observations.all.size() / 2].RA += 200 * Constants::ARCSEC_RAD;
	data.observations.all[data.observations.all.size() / 3].RA += 20 * Constants::ARCSEC_RAD;
	data.config.sigmaClip.enabled = true;
	data.config.sigmaClip.maxPasses = 1;
	const auto result = oc::FitSession{data.config}.optimize({data.initial, data.observations.all});
	ASSERT_TRUE(result.orbit);
	EXPECT_EQ(result.sigmaClipPassesRun, 1);
	EXPECT_FALSE(result.converged);
	EXPECT_EQ(result.selectionStatus, oc::FitSelectionStatus::PassLimit);
	EXPECT_EQ(result.residuals.count(), result.keptObsIndices.size());
	EXPECT_EQ(result.sigmaClipRejected, data.observations.all.size() - result.keptObsIndices.size());
	auto fixed = oc::FitRequest{result.fitEpochOrbit, data.observations.all};
	std::vector<std::uint8_t> mask(data.observations.all.size(), 0);
	for (const auto index : result.keptObsIndices)
		mask[index] = 1;
	fixed.initialIncluded = mask;
	fixed.options.evaluateOnly = true;
	const auto evaluated = oc::FitSession{data.config}.optimize(fixed);
	EXPECT_NEAR(evaluated.cost, result.cost, 1e-10);
}

TEST_F(FitSessionTest, ShortArcRejectionDoesNotRequireTwentyObservations)
{
	auto observations = oc::loadObservations(getExecutableDir() / "examples/NEO/impacts/2018_LA.xml");
	ASSERT_EQ(observations.all.size(), 18u);
	oc::prepareObservations(observations, {.debias = true});
	oc::FitConfig config;
	config.useWeights = true;
	config.obsWeighting.useObserverSigmas = true;
	config.nbody.enabled = true;
	config.nbody.planetMask = 4095;
	config.nbody.useSb441 = true;
	config.sigmaClip.enabled = true;

	// Reproduce the old stopping condition without changing the fit or its data.
	config.sigmaClip.minObsKept = 20;
	const auto blocked = oc::FitPipeline{config}.run({observations});
	ASSERT_TRUE(blocked.fit.fitEpochOrbit);
	EXPECT_FALSE(blocked.fit.converged);
	EXPECT_EQ(blocked.fit.selectionStatus, oc::FitSelectionStatus::InsufficientObservations);
	EXPECT_EQ(blocked.fit.sigmaClipPassesRun, 0);
	EXPECT_EQ(blocked.fit.keptObsIndices.size(), 17u);
	EXPECT_LT(blocked.fit.iterations, 100);

	config.sigmaClip.minObsKept = oc::SigmaClipSettings{}.minObsKept;
	const auto fitted = oc::FitPipeline{config}.run({observations});
	ASSERT_TRUE(fitted.fit.fitEpochOrbit);
	EXPECT_TRUE(fitted.fit.converged);
	EXPECT_EQ(fitted.fit.selectionStatus, oc::FitSelectionStatus::Stable);
	EXPECT_GT(fitted.fit.sigmaClipPassesRun, 0);
	EXPECT_LT(fitted.fit.keptObsIndices.size(), blocked.fit.keptObsIndices.size());
	EXPECT_GE(fitted.fit.keptObsIndices.size(), 11u);
	EXPECT_TRUE(fitted.fit.uncertainty.hasStateCovariance());
	EXPECT_EQ(fitted.fit.residuals.count(), fitted.fit.keptObsIndices.size());

	std::vector<std::uint8_t> mask(observations.all.size(), 0);
	for (const auto index : fitted.fit.keptObsIndices)
		mask[index] = 1;
	const auto evaluated =
	    oc::FitSession{config}.optimize({.initialOrbit = fitted.fit.fitEpochOrbit,
	                                     .observations = observations.all,
	                                     .initialIncluded = mask,
	                                     .options = {.selection = oc::FitSelectionMode::Fixed, .evaluateOnly = true}});
	EXPECT_NEAR(evaluated.cost, fitted.fit.cost, 1e-10);
	EXPECT_EQ(evaluated.keptObsIndices, fitted.fit.keptObsIndices);
	EXPECT_EQ(evaluated.selectionStatus, oc::FitSelectionStatus::NotApplied);
	EXPECT_TRUE(evaluated.uncertainty.fitParameterCovarianceFormal.isApprox(
	    fitted.fit.uncertainty.fitParameterCovarianceFormal, 1e-9));
}

TEST_F(FitSessionTest, SunOnlyNonGravRetainsNumericalBootstrap)
{
	auto data = fitContractCase();
	data.config.nbody.enabled = false;
	data.config.nonGrav.enabled = true;
	data.config.nonGrav.model = oc::NonGravModelKind::Yarkovsky;
	data.config.sigmaClip.enabled = true;
	int bootstraps = 0, transports = 0;
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		bootstraps += message.starts_with("=== Non-grav bootstrap: Sun-only numerical");
		transports += message.starts_with("==== ORBIT AT TARGET EPOCH");
	});
	const auto result = oc::FitSession{data.config}.optimize({.initialOrbit = data.initial,
	                                                          .observations = data.observations.all,
	                                                          .outputEpochJdTdb = data.initial->getEpoch() + 30.0});
	ASSERT_TRUE(result.orbit);
	EXPECT_EQ(bootstraps, 1);
	EXPECT_EQ(transports, 1);
	EXPECT_TRUE(result.finalNonGrav.enabled());
	EXPECT_TRUE(result.uncertainty.hasStateCovariance());
}

TEST_F(FitSessionTest, OutputFinalizationReusesAcceptedLinearization)
{
	const auto data = fitContractCase();
	OrbitOptimizerConfig config;
	config.nbody = data.config.nbody;
	config.useWeights = true;
	config.weightProfile = data.config.weightProfile;
	config.obsWeighting = data.config.obsWeighting;
	config.computeUncertainty = false;
	OrbitOptimizer optimizer(data.initial, oc::fit::prepareFitObservations(data.observations.all, config), config);
	const auto fitted = optimizer.optimize();
	ASSERT_TRUE(fitted.fitEpochOrbit);
	int lmStarts = 0;
	const FitLogCallback log(
	    [&](LogLevel, const std::string& message) { lmStarts += message.starts_with(solverStart); });
	const double epoch = data.initial->getEpoch() + 30.0;
	const auto first = optimizer.finalizeOutput(epoch, true);
	ASSERT_TRUE(first.uncertainty.hasStateCovariance());
	ASSERT_TRUE(optimizer.finalizeOutput(epoch + 30.0, true).orbit);
	const auto returned = optimizer.finalizeOutput(epoch, true);
	EXPECT_EQ(lmStarts, 0);
	EXPECT_EQ(returned.fitEpochOrbit->getEpoch(), fitted.fitEpochOrbit->getEpoch());
	EXPECT_EQ(returned.fitEpochOrbit->getR().x, fitted.fitEpochOrbit->getR().x);
	EXPECT_EQ(returned.fitEpochOrbit->getV().z, fitted.fitEpochOrbit->getV().z);
	EXPECT_EQ(returned.cost, fitted.cost);
	EXPECT_EQ(returned.orbit->getR().x, first.orbit->getR().x);
	EXPECT_TRUE((returned.uncertainty.elementCovariance.array() == first.uncertainty.elementCovariance.array()).all());
	const auto withoutCovariance = optimizer.finalizeOutput(0.0, false);
	EXPECT_TRUE(withoutCovariance.uncertainty.empty());
	EXPECT_EQ(withoutCovariance.orbit->getEpoch(), fitted.fitEpochOrbit->getEpoch());
	std::stop_source cancel;
	cancel.request_stop();
	const auto stopped = optimizer.finalizeOutput(epoch, true, {.cancel = cancel.get_token()});
	EXPECT_TRUE(stopped.cancelled);
	EXPECT_FALSE(stopped.orbit);
	EXPECT_THROW((void)optimizer.finalizeOutput(epoch, true), std::logic_error);
}

TEST_F(FitSessionTest, CancellationAtLmBoundariesAndCovariancePublishesNoSolution)
{
	auto data = fitContractCase();
	for (const std::string marker : {solverStart, solverFinish, "[Covariance]"}) {
		std::stop_source cancel;
		bool reached = false;
		int progressAfterCancellation = 0;
		const FitLogCallback log([&](LogLevel, const std::string& message) {
			if (message.starts_with(marker)) {
				reached = true;
				cancel.request_stop();
			}
		});
		oc::FitSession session{data.config};
		const auto result = session.optimize({.initialOrbit = data.initial, .observations = data.observations.all},
		                                     {.cancel = cancel.get_token(), .progress = [&](const auto&) {
			                                      if (cancel.stop_requested())
				                                      ++progressAfterCancellation;
		                                      }});
		EXPECT_TRUE(reached) << marker;
		EXPECT_TRUE(result.cancelled) << marker;
		EXPECT_FALSE(result.converged);
		EXPECT_FALSE(result.orbit);
		EXPECT_FALSE(result.fitEpochOrbit);
		EXPECT_EQ(result.residuals.count(), 0u);
		EXPECT_FALSE(result.uncertainty.hasElementCovariance);
		EXPECT_EQ(progressAfterCancellation, 0);
	}
}

TEST_F(FitSessionTest, CancellationDuringBootstrapRestoresNonGravForNextCall)
{
	const auto data = fitContractCase();
	OrbitOptimizerConfig config;
	config.nbody = data.config.nbody;
	config.nonGrav.enabled = true;
	config.nonGrav.model = oc::NonGravModelKind::Yarkovsky;
	config.nonGrav.useInitialValues = true;
	config.nonGrav.initialA2 = -3e-14;
	config.computeUncertainty = false;
	config.logCacheStats = false;
	OrbitOptimizer optimizer(data.initial, oc::fit::prepareFitObservations(data.observations.all, config), config);
	std::stop_source cancel;
	int bootstrapStarts = 0;
	int extendedStarts = 0;
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		bootstrapStarts += message.starts_with("=== Non-grav bootstrap");
		extendedStarts += message.starts_with("=== Non-grav extended");
		if (message.starts_with(solverFinish))
			cancel.request_stop();
	});
	const auto first = optimizer.optimize({.cancel = cancel.get_token()});
	ASSERT_TRUE(first.cancelled);
	EXPECT_FALSE(optimizer.postFitDiagnostics().valid);
	cancel = std::stop_source{};
	const auto second = optimizer.optimize({.cancel = cancel.get_token()});
	EXPECT_TRUE(second.cancelled);
	EXPECT_EQ(bootstrapStarts, 2);
	EXPECT_EQ(extendedStarts, 0);
}

TEST_F(FitSessionTest, CancellationAfterSelectionSkipsFinalEvaluation)
{
	auto data = fitContractCase();
	data.config.sigmaClip.enabled = true;
	std::stop_source cancel;
	int preparations = 0;
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		preparations += message.starts_with("ObservationPipeline: ready");
		if (message.starts_with("[FitSession]") && message.find("status=NO_CHANGE") != std::string::npos)
			cancel.request_stop();
	});
	const auto result = oc::FitSession{data.config}.optimize(
	    {.initialOrbit = data.initial, .observations = data.observations.all}, {.cancel = cancel.get_token()});
	EXPECT_TRUE(cancel.stop_requested());
	EXPECT_TRUE(result.cancelled);
	EXPECT_FALSE(result.orbit);
	EXPECT_EQ(preparations, 1);
}

TEST_F(FitSessionTest, PipelineDoesNotPublishCancelledFinalFit)
{
	const auto data = fitContractCase();
	std::stop_source cancel;
	bool finalFit = false;
	oc::FitRunControl control{.cancel = cancel.get_token()};
	control.progress = [&](const std::string& message) {
		if (message.starts_with("Optimizing "))
			finalFit = true;
	};
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		if (finalFit && message.starts_with(solverFinish))
			cancel.request_stop();
	});
	const auto result = oc::FitPipeline{data.config}.run({data.observations, {}}, control);
	EXPECT_TRUE(result.fit.cancelled);
	EXPECT_FALSE(result.fit.orbit);
	EXPECT_FALSE(result.fit.fitEpochOrbit);
	EXPECT_TRUE(result.fit.uncertainty.empty());
	EXPECT_EQ(result.fit.residuals.count(), 0u);
	EXPECT_TRUE(result.observations.all.empty());
	EXPECT_TRUE(result.observations.iod.empty());
	EXPECT_FALSE(result.initial.orbit);
	EXPECT_FALSE(result.extension);
	EXPECT_EQ(result.fitInputObservations, 0u);
	EXPECT_TRUE(finalFit);
	EXPECT_TRUE(cancel.stop_requested());
}

TEST_F(FitSessionTest, RunControlIsNotRetainedBetweenCalls)
{
	const auto data = fitContractCase();
	const oc::FitSession session{data.config};
	const oc::FitRequest request{.initialOrbit = data.initial,
	                             .observations = data.observations.all,
	                             .options = {.uncertainty = oc::FitUncertaintyMode::Skip, .evaluateOnly = true}};
	int progressCalls = 0;
	{
		std::stop_source cancel;
		const auto cancelled = session.optimize(request, {.cancel = cancel.get_token(), .progress = [&](const auto&) {
			                                                  ++progressCalls;
			                                                  cancel.request_stop();
		                                                  }});
		EXPECT_TRUE(cancelled.cancelled);
		EXPECT_FALSE(cancelled.orbit);
		EXPECT_TRUE(cancelled.uncertainty.empty());
	}
	EXPECT_EQ(progressCalls, 1);
	const auto result = session.optimize(request);
	EXPECT_FALSE(result.cancelled);
	ASSERT_TRUE(result.orbit);
	EXPECT_EQ(progressCalls, 1);
	EXPECT_THROW(
	    (void)session.optimize(request, {.progress = [](const auto&) { throw std::logic_error("callback failure"); }}),
	    std::logic_error);
	const auto repeated = session.optimize(request);
	ASSERT_TRUE(repeated.orbit);
	EXPECT_DOUBLE_EQ(repeated.cost, result.cost);
	EXPECT_DOUBLE_EQ(repeated.orbit->getR().x, result.orbit->getR().x);
}

TEST_F(FitSessionTest, OptimizerRunControlIsNotRetainedBetweenCalls)
{
	const auto data = fitContractCase();
	OrbitOptimizerConfig config;
	config.nbody = data.config.nbody;
	config.evaluateOnly = true;
	config.computeUncertainty = false;
	OrbitOptimizer optimizer(data.initial, oc::fit::prepareFitObservations(data.observations.all, config), config);
	int calls = 0;
	{
		std::stop_source cancel;
		const auto& cancelled = optimizer.optimize({.cancel = cancel.get_token(), .progress = [&](const auto&) {
			                                            ++calls;
			                                            cancel.request_stop();
		                                            }});
		EXPECT_TRUE(cancelled.cancelled);
		EXPECT_FALSE(cancelled.orbit);
	}
	const auto& result = optimizer.optimize();
	EXPECT_FALSE(result.cancelled);
	EXPECT_TRUE(result.orbit);
	EXPECT_EQ(calls, 1);
}

TEST_F(FitSessionTest, InitialOrbitCallbackCanCancelBeforeExtension)
{
	const auto data = fitContractCase();
	const oc::FitPipeline pipeline{data.config};
	std::stop_source cancel;
	std::vector<std::string> messages;
	int initialCalls = 0;
	const auto result = pipeline.run(
	    {data.observations, {}},
	    {.cancel = cancel.get_token(), .progress = [&](const auto& message) { messages.push_back(message); }},
	    [&](const oc::InitialOrbitInfo& initial) {
		    EXPECT_TRUE(initial.orbit);
		    ++initialCalls;
		    cancel.request_stop();
	    });
	EXPECT_TRUE(result.fit.cancelled);
	EXPECT_FALSE(result.initial.orbit);
	EXPECT_FALSE(result.extension);
	EXPECT_EQ(initialCalls, 1);
	EXPECT_EQ(messages, std::vector<std::string>{"Starting IOD"});
	EXPECT_THROW((void)pipeline.run({data.observations, {}}, {},
	                                [](const auto&) { throw std::logic_error("initial orbit callback failure"); }),
	             std::logic_error);
	const auto repeated = pipeline.run({data.observations, {}});
	EXPECT_FALSE(repeated.fit.cancelled);
	EXPECT_TRUE(repeated.fit.orbit);
	EXPECT_EQ(initialCalls, 1);
	EXPECT_EQ(messages.size(), 1u);
}

TEST_F(FitSessionTest, RejectsNullInitialOrbit)
{
	oc::FitSession session{oc::FitConfig{}};
	std::vector<ObsNode> obs(3);

	EXPECT_THROW((void)session.optimize({.initialOrbit = nullptr, .observations = obs}), std::invalid_argument);
}

TEST_F(FitSessionTest, RejectsTooFewObservations)
{
	oc::FitSession session{oc::FitConfig{}};
	auto orbit = std::make_shared<Orbit>();
	std::vector<ObsNode> obs(2);
	EXPECT_THROW((void)session.optimize({.initialOrbit = orbit, .observations = obs}), std::invalid_argument);
}

TEST_F(FitSessionTest, Loads2020VvAdesUsingNativeCatalogPath)
{
	ASSERT_TRUE(std::filesystem::is_regular_file(getExecutableDir() / "OBSCODE.dat"));
	const auto input =
	    std::filesystem::path(__FILE__).parent_path() / "data" / "2020VV_ADES_2022.ESA11_table_replay.xml";
	ASSERT_TRUE(std::filesystem::is_regular_file(input));
	const auto loaded = oc::loadObservations(input);
	EXPECT_EQ(loaded.all.size(), 78u);
	EXPECT_GE(loaded.iod.size(), 3u);
	EXPECT_FALSE(loaded.objectName.empty());
	EXPECT_EQ(loaded.source, input);
	EXPECT_EQ(loaded.format, "ADES");
	EXPECT_EQ(loaded.all.front().astCat, "Gaia2");
	EXPECT_EQ(loaded.iod.front().astCat, "Gaia2");
	EXPECT_FALSE(loaded.preparation);
	EXPECT_TRUE(
	    std::is_sorted(loaded.all.begin(), loaded.all.end(), [](const auto& a, const auto& b) { return a.JD < b.JD; }));
	EXPECT_EQ(std::count_if(loaded.all.begin(), loaded.all.end(),
	                        [](const auto& observation) { return observation.hasObserverSigma; }),
	          11);
}

TEST_F(FitSessionTest, PreparationWithoutDebiasPreservesBothPools)
{
	const auto input = std::filesystem::path(__FILE__).parent_path() / "data/2020VV_ADES_2022.ESA11_table_replay.xml";
	auto loaded = oc::loadObservations(input);
	const auto original = loaded;
	oc::prepareObservations(loaded);
	ASSERT_TRUE(loaded.preparation);
	EXPECT_FALSE(loaded.preparation->debiasRequested);
	EXPECT_TRUE(loaded.preparation->biasTable.empty());
	EXPECT_EQ(loaded.preparation->correctedAll, 0u);
	EXPECT_EQ(loaded.preparation->correctedIod, 0u);
	const auto check = [](const auto& actual, const auto& expected) {
		ASSERT_EQ(actual.size(), expected.size());
		for (std::size_t i = 0; i < actual.size(); ++i) {
			EXPECT_EQ(actual[i].JD, expected[i].JD);
			EXPECT_EQ(actual[i].RA, expected[i].RA);
			EXPECT_EQ(actual[i].DEC, expected[i].DEC);
			EXPECT_EQ(actual[i].astCat, expected[i].astCat);
			EXPECT_FALSE(actual[i].submittedCoordinates);
			EXPECT_EQ(actual[i].biasRaArcsec, 0.0);
			EXPECT_EQ(actual[i].biasDecArcsec, 0.0);
		}
	};
	check(loaded.all, original.all);
	check(loaded.iod, original.iod);
}

TEST_F(FitSessionTest, PreparationMatchesSelectedDebiasForBothPools)
{
	const auto input = std::filesystem::path(__FILE__).parent_path() / "data/2020VV_ADES_2022.ESA11_table_replay.xml";
	auto loaded = oc::loadObservations(input);
	auto expected = loaded;
	const auto table =
	    Calculations::astrometryBiasPathForProfile(getExecutableDir(), oc::AstrometricWeightProfile::Vfcc17);
	ASSERT_TRUE(std::filesystem::is_regular_file(table));
	const auto allCount = Calculations::debiasObservations(expected.all, table);
	const auto iodCount = Calculations::debiasObservations(expected.iod, table);
	ASSERT_GT(allCount, 0u);
	oc::prepareObservations(loaded, {.debias = true,
	                                 .weightProfile = oc::AstrometricWeightProfile::Vfcc17,
	                                 .runtimeDirectory = getExecutableDir()});
	ASSERT_TRUE(loaded.preparation);
	EXPECT_TRUE(loaded.preparation->debiasRequested);
	EXPECT_EQ(loaded.preparation->biasTable, table);
	EXPECT_EQ(loaded.preparation->correctedAll, allCount);
	EXPECT_EQ(loaded.preparation->correctedIod, iodCount);
	const auto check = [](const auto& actual, const auto& reference) {
		ASSERT_EQ(actual.size(), reference.size());
		for (std::size_t i = 0; i < actual.size(); ++i) {
			EXPECT_EQ(actual[i].RA, reference[i].RA);
			EXPECT_EQ(actual[i].DEC, reference[i].DEC);
			EXPECT_EQ(actual[i].submittedRaRad(), reference[i].submittedRaRad());
			EXPECT_EQ(actual[i].submittedDecRad(), reference[i].submittedDecRad());
			EXPECT_EQ(actual[i].biasRaArcsec, reference[i].biasRaArcsec);
			EXPECT_EQ(actual[i].biasDecArcsec, reference[i].biasDecArcsec);
			EXPECT_EQ(actual[i].vector.x, reference[i].vector.x);
			EXPECT_EQ(actual[i].vector.y, reference[i].vector.y);
			EXPECT_EQ(actual[i].vector.z, reference[i].vector.z);
			EXPECT_EQ(actual[i].utcMjdDay, reference[i].utcMjdDay);
			EXPECT_EQ(actual[i].utcSecondsOfDay, reference[i].utcSecondsOfDay);
			EXPECT_EQ(actual[i].rmsRaArcsec, reference[i].rmsRaArcsec);
			EXPECT_EQ(actual[i].rmsDecArcsec, reference[i].rmsDecArcsec);
		}
	};
	check(loaded.all, expected.all);
	check(loaded.iod, expected.iod);
	EXPECT_THROW(oc::prepareObservations(loaded, {.debias = true}), std::logic_error);
	check(loaded.all, expected.all);
	check(loaded.iod, expected.iod);
}

TEST_F(FitSessionTest, EvaluateOnlyPreservesAdesPhotometryThroughThePipeline)
{
	const auto input =
	    std::filesystem::path(__FILE__).parent_path() / "data" / "2020VV_ADES_2022.ESA11_table_replay.xml";
	auto loaded = oc::loadObservations(input);
	ASSERT_GE(loaded.all.size(), 2u);
	loaded.all[0].nuclearMagnitude = true;
	loaded.all[1].nuclearMagnitude = false;
	const auto& all = loaded.all;
	ASSERT_EQ(all.size(), 78u);
	std::size_t reported = 0;
	for (const auto& obs : all) {
		if (obs.rmsMag > 0.0)
			++reported;
	}
	EXPECT_EQ(reported, 18u);

	auto initial = make_orbit_from_states_vectors(StateVectors(Vec3D(1.0, 0.1, 0.02), Vec3D(0.01, 0.9, 0.1)),
	                                              all.front().epochJdTdb());
	oc::FitConfig config;
	config.weightProfile = oc::AstrometricWeightProfile::Vfcc17;
	oc::FitRunOptions options;
	options.selection = oc::FitSelectionMode::Fixed;
	options.uncertainty = oc::FitUncertaintyMode::Skip;
	options.evaluateOnly = true;
	const std::vector<std::uint8_t> included(all.size(), 1);
	const auto result = oc::FitSession(config).optimize(
	    {.initialOrbit = initial, .observations = all, .initialIncluded = included, .options = options});
	ASSERT_EQ(result.residuals.count(), all.size());
	for (std::size_t i = 0; i < all.size(); ++i) {
		const auto& row = result.residuals.data()[i];
		EXPECT_DOUBLE_EQ(row.mag, all[i].MAG);
		EXPECT_EQ(row.magBand, all[i].magBand);
		EXPECT_EQ(row.nuclearMagnitude, all[i].nuclearMagnitude);
		EXPECT_DOUBLE_EQ(row.magSigma, oc::photometry::resolveSigma(all[i].magString, all[i].rmsMag));
	}
}

TEST_F(FitSessionTest, FixedRequestPreservesMaskAndNonGravWarmStart)
{
	std::vector<ObsNode> observations(8);
	for (std::size_t i = 0; i < observations.size(); ++i) {
		observations[i].JD = 2459165.5 + 0.3 * i;
		observations[i].RA = 1.0 + 0.01 * i;
		observations[i].DEC = 0.1;
		observations[i].obsCode = "500";
	}
	const auto initial =
	    make_orbit_from_states_vectors(StateVectors(Vec3D(1.4, .4, .1), Vec3D(-.004, .012, .001) / Constants::GAUSS_K),
	                                   observations.front().epochJdTdb());
	ASSERT_TRUE(initial);
	NonGravModel nonGrav;
	nonGrav.setEnabled(true);
	nonGrav.setModel(NgModel::Yarkovsky);
	nonGrav.setParams(0.0, -3.0e-14, 0.0);
	const std::array<std::uint8_t, 8> mask{2, 0, 1, 0, 255, 0, 1, 0};
	oc::FitConfig config;
	config.nbody.enabled = true;
	config.nonGrav.enabled = true;
	config.nonGrav.model = oc::NonGravModelKind::Yarkovsky;
	config.sigmaClip.enabled = true;
	const oc::FitRequest request{.initialOrbit = initial,
	                             .observations = observations,
	                             .initialIncluded = mask,
	                             .warmStartNonGrav = &nonGrav,
	                             .options = {.selection = oc::FitSelectionMode::Fixed,
	                                         .uncertainty = oc::FitUncertaintyMode::Skip,
	                                         .evaluateOnly = true}};
	const auto result = oc::FitSession(config).optimize(request);
	ASSERT_TRUE(result.orbit);
	EXPECT_EQ(result.keptObsIndices, (std::vector<std::size_t>{0, 2, 4, 6}));
	EXPECT_EQ(result.sigmaClipPassesRun, 0);
	EXPECT_EQ(result.sigmaClipRejected, 4);
	EXPECT_EQ(result.residuals.count(), 4u);
	EXPECT_EQ(result.allResiduals.count(), observations.size());
	EXPECT_EQ(result.iterations, 0);
	EXPECT_FALSE(result.uncertainty.hasElementCovariance);
	EXPECT_TRUE(result.finalNonGrav.enabled());
	EXPECT_EQ(result.finalNonGrav.model(), NgModel::Yarkovsky);
	EXPECT_DOUBLE_EQ(result.finalNonGrav.A2(), nonGrav.A2());
	EXPECT_DOUBLE_EQ(nonGrav.A2(), -3.0e-14);
	EXPECT_EQ(mask, (std::array<std::uint8_t, 8>{2, 0, 1, 0, 255, 0, 1, 0}));
}

TEST_F(FitSessionTest, BorisovPipelineStartsFreshNonGravSelection)
{
	const auto path = getExecutableDir() / "examples/comets/2I_Borisov.xml";
	ASSERT_TRUE(std::filesystem::exists(path));
	auto loaded = oc::loadObservations(path);
	oc::FitConfig cfg;
	cfg.useWeights = true;
	cfg.weightProfile = oc::AstrometricWeightProfile::Vfcc17;
	cfg.obsWeighting.useObserverSigmas = true;
	cfg.nbody.enabled = true;
	cfg.nbody.planetMask = 2303;
	cfg.nbody.useSb441 = true;
	cfg.sigmaClip.enabled = true;
	cfg.nonGrav.enabled = true;
	cfg.nonGrav.model = oc::NonGravModelKind::MarsdenH2O;
	cfg.nonGrav.solveMask = 15;
	oc::prepareObservations(loaded, {.debias = true, .weightProfile = cfg.weightProfile});
	bool fittedNonGrav = false;
	int checksAfterNonGrav = 0;
	const FitLogCallback log([&](LogLevel, const std::string& message) {
		if (message.starts_with("Non-grav done."))
			fittedNonGrav = true;
		if (fittedNonGrav && message.starts_with("[HystereticRejector] rejectChi2="))
			++checksAfterNonGrav;
	});
	const auto result = oc::FitPipeline{cfg}.run({loaded});
	ASSERT_TRUE(result.fit.orbit);
	EXPECT_EQ(result.fitInputObservations, loaded.all.size());
	EXPECT_EQ(result.selection, oc::FitSelectionMode::Configured);
	EXPECT_TRUE(result.fit.converged);
	EXPECT_GT(result.fit.sigmaClipPassesRun, 0);
	EXPECT_GT(checksAfterNonGrav, 0);
	ASSERT_EQ(result.fit.allResiduals.count(), loaded.all.size());
	for (std::size_t i = 0; i < loaded.all.size(); ++i) {
		const auto& metric = result.fit.allResiduals.data()[i].astrometricWeight;
		ASSERT_TRUE(metric && metric->postFitChiSquared);
		const bool included = std::ranges::find(result.fit.keptObsIndices, i) != result.fit.keptObsIndices.end();
		const double sigma = cfg.sigmaClip.thresholdNSigma;
		if (included)
			EXPECT_LE(*metric->postFitChiSquared, sigma * sigma);
		else
			EXPECT_GT(*metric->postFitChiSquared, sigma * sigma * Calculations::HystereticRejectConfig::recoveryRatio);
	}
}

TEST_F(FitSessionTest, BorisovSelectionChecksEveryCompletedCorrection)
{
	const auto path = getExecutableDir() / "examples/comets/2I_Borisov.xml";
	ASSERT_TRUE(std::filesystem::exists(path));
	for (const auto [sigma, nonGrav] : {std::pair{2.0, true}, {3.0, true}, {3.0, false}}) {
		SCOPED_TRACE(::testing::Message() << "sigma=" << sigma << " nonGrav=" << nonGrav);
		auto loaded = oc::loadObservations(path);
		oc::FitConfig cfg;
		cfg.useWeights = true;
		cfg.weightProfile = oc::AstrometricWeightProfile::Vfcc17;
		cfg.obsWeighting.useObserverSigmas = true;
		cfg.nbody.enabled = true;
		cfg.nbody.planetMask = 2303;
		cfg.nbody.useSb441 = true;
		cfg.sigmaClip.enabled = true;
		cfg.sigmaClip.thresholdNSigma = sigma;
		cfg.sigmaClip.maxPasses = 15;
		cfg.nonGrav.enabled = nonGrav;
		cfg.nonGrav.model = oc::NonGravModelKind::MarsdenH2O;
		cfg.nonGrav.solveMask = 15;
		oc::prepareObservations(loaded, {.debias = true, .weightProfile = cfg.weightProfile});
		const FitLogCallback log([&](LogLevel, const std::string& message) {
			EXPECT_EQ(message.find("selection did not stabilize"), std::string::npos);
		});
		const auto result = oc::FitPipeline{cfg}.run({loaded});
		ASSERT_TRUE(result.fit.converged);
		EXPECT_EQ(result.fit.selectionStatus, oc::FitSelectionStatus::Stable);
		EXPECT_GT(result.fit.sigmaClipPassesRun, 0);
		EXPECT_LE(result.fit.sigmaClipPassesRun, cfg.sigmaClip.maxPasses);
		ASSERT_EQ(result.fit.allResiduals.count(), loaded.all.size());
		std::vector<std::uint8_t> mask(loaded.all.size(), 0);
		for (const auto index : result.fit.keptObsIndices)
			mask[index] = 1;
		std::vector<double> chiSquared;
		chiSquared.reserve(mask.size());
		for (std::size_t i = 0; i < mask.size(); ++i) {
			const auto& metric = result.fit.allResiduals.data()[i].astrometricWeight;
			ASSERT_TRUE(metric && metric->postFitChiSquared);
			const double chi2 = *metric->postFitChiSquared;
			chiSquared.push_back(chi2);
			if (mask[i])
				EXPECT_LE(chi2, sigma * sigma);
			else
				EXPECT_GT(chi2, sigma == 2.0 ? 3.684 : 8.289);
		}
		const auto original = mask;
		Calculations::HystereticRejectConfig rejection;
		rejection.rejectChiSquared = sigma * sigma;
		rejection.recoverChiSquared = sigma == 2.0 ? 3.684 : 8.289;
		const auto check = Calculations::applyHystereticReject(loaded.all, mask, chiSquared, rejection);
		EXPECT_EQ(check.status, Calculations::RejectStatus::NoChange);
		EXPECT_EQ(mask, original);
		EXPECT_DOUBLE_EQ(check.appliedRejectChiSquared, sigma * sigma);
		const std::string label = !nonGrav ? "gravitySigma3" : sigma == 2.0 ? "sigma2" : "sigma3";
		RecordProperty(label + "Kept", result.fit.keptObsIndices.size());
		RecordProperty(label + "Refits", result.fit.sigmaClipPassesRun);
	}
}

TEST_F(FitSessionTest, ExternalFormalCovarianceUsesMaskAndPhysicalVelocityUnits)
{
	std::vector<ObsNode> observations(8);
	for (std::size_t i = 0; i < observations.size(); ++i) {
		observations[i].JD = 2459165.5 + 0.3 * i;
		observations[i].RA = 1.0 + 0.01 * i;
		observations[i].DEC = 0.1;
		observations[i].obsCode = "500";
	}
	const Vec3D r(1.4, .4, .1), v(-.004, .012, .001);
	const double epoch = observations.front().epochJdTdb();
	auto orbit = make_orbit_from_states_vectors(StateVectors(r, v / Constants::GAUSS_K), epoch);
	ASSERT_TRUE(orbit);
	const std::array<double, 6> physical{r.x, r.y, r.z, v.x, v.y, v.z};
	OrbitOptimizerConfig config;
	config.evaluateOnly = true;
	config.useWeights = true;
	config.nbody.enabled = true;
	config.computeUncertainty = false;
	config.logDiagnostics = false;
	config.weightProfile = oc::AstrometricWeightProfile::Uniform;
	OrbitOptimizer optimizer(orbit, oc::fit::prepareFitObservations(observations, config), config);
	const auto& direct = optimizer.optimize();
	ASSERT_TRUE(direct.orbit);
	EXPECT_FALSE(direct.converged);
	EXPECT_EQ(direct.iterations, 0);
	oc::FitConfig fitConfig;
	fitConfig.nbody.enabled = true;
	fitConfig.weightProfile = oc::AstrometricWeightProfile::Uniform;
	oc::FitRunOptions options;
	options.selection = oc::FitSelectionMode::Fixed;
	options.uncertainty = oc::FitUncertaintyMode::Skip;
	options.evaluateOnly = true;
	const auto evaluated =
	    oc::FitSession(fitConfig).optimize({.initialOrbit = orbit, .observations = observations, .options = options});
	ASSERT_TRUE(evaluated.orbit);
	EXPECT_FALSE(evaluated.converged);
	EXPECT_EQ(evaluated.iterations, 0);
	oc::FitParameterCovariance formal = oc::FitParameterCovariance::Zero(6, 6);
	formal.diagonal() << 1e-16, 2e-16, 3e-16, 1e-17, 2e-17, 3e-17;
	const std::array<std::uint8_t, 8> mask{1, 0, 1, 0, 1, 0, 1, 0};
	const auto diagnostics = optimizer.postFitDiagnostics(formal, mask);
	ASSERT_TRUE(diagnostics.valid);
	EXPECT_TRUE(diagnostics.formalCovariance.isApprox(formal, 0.0));
	const auto linear = optimizer.weightedObservationLinearizationPhysical(physical);
	ASSERT_TRUE(linear.valid);
	for (std::size_t i = 0; i < observations.size(); ++i) {
		const auto& metric = diagnostics.residuals.data()[i].astrometricWeight;
		ASSERT_TRUE(metric);
		ASSERT_TRUE(metric->postFitChiSquared);
		const Eigen::MatrixXd jac = linear.jacobian.middleRows(2 * i, 2);
		const Eigen::Matrix2d covariance =
		    Eigen::Matrix2d::Identity() + (mask[i] ? -1.0 : 1.0) * jac * formal * jac.transpose();
		const Eigen::Vector2d residual(linear.residuals[2 * i], linear.residuals[2 * i + 1]);
		const double expected =
		    (covariance(1, 1) * residual.x() * residual.x() - 2.0 * covariance(0, 1) * residual.x() * residual.y() +
		     covariance(0, 0) * residual.y() * residual.y()) /
		    (covariance(0, 0) * covariance(1, 1) - covariance(0, 1) * covariance(1, 0));
		EXPECT_NEAR(*metric->postFitChiSquared, expected, 1e-8 * expected);
		if (mask[i])
			EXPECT_GT(*metric->postFitChiSquared, metric->chiSquared);
		else
			EXPECT_LT(*metric->postFitChiSquared, metric->chiSquared);
	}
}

TEST_F(FitSessionTest, ObservationLinearizationRetainsSplitUtcPrecision)
{
	std::vector<ObsNode> observations(3);
	for (std::size_t i = 0; i < observations.size(); ++i) {
		auto& o = observations[i];
		o.hasSplitUtc = true;
		o.utcMjdDay = 59139;
		o.utcSecondsOfDay = 21600.0 + static_cast<double>(i) * 0.000005;
		o.JD = 2400000.5 + (o.utcMjdDay + o.utcSecondsOfDay / 86400.0);
		o.obsCode = "500";
		o.RA = 1.0;
		o.DEC = .1;
		ASSERT_EQ(o.JD, observations.front().JD);
	}
	const Vec3D r(1.4, .4, .1), v(-.004, .012, .001);
	const double epoch = observations.front().epochJdTdb();
	auto orbit = make_orbit_from_states_vectors(StateVectors(r, v / Constants::GAUSS_K), epoch);
	const std::array<double, 6> physical{r.x, r.y, r.z, v.x, v.y, v.z};
	for (const bool useNBody : {false, true}) {
		OrbitOptimizerConfig config;
		config.nbody.enabled = useNBody;
		config.logDiagnostics = false;
		config.weightProfile = oc::AstrometricWeightProfile::Uniform;
		OrbitOptimizer optimizer(orbit, oc::fit::prepareFitObservations(observations, config), config);
		const auto linear = optimizer.weightedObservationLinearizationPhysical(physical);
		ASSERT_TRUE(linear.valid);
		ASSERT_EQ(linear.residuals.size(), 6u);
		const double difference =
		    std::hypot(linear.residuals[4] - linear.residuals[0], linear.residuals[5] - linear.residuals[1]);
		EXPECT_GT(difference, 1e-14) << "Distinct split UTC times collapsed into full JD";
	}
}

TEST_F(FitSessionTest, OptimizerResultOwnsDataAndSeparatesFitAndOutputEpochs)
{
	const auto input =
	    std::filesystem::path(__FILE__).parent_path() / "data" / "2020VV_ADES_2022.ESA11_table_replay.xml";
	const auto loaded = oc::loadObservations(input);
	constexpr double fitEpoch = 2459170.2715327661;
	const auto initial = make_orbit_from_states_vectors(
	    StateVectors(Vec3D(.62434444840234238, .83392068676583608, .0034709928450530553),
	                 Vec3D(-.012853569243756456, .011750580852295998, .000093051552130828674) / Constants::GAUSS_K),
	    fitEpoch);
	ASSERT_TRUE(initial);
	const Vec3D originalPosition = initial->getR();

	for (const bool numerical : {false, true}) {
		oc::FitResult saved;
		{
			OrbitOptimizerConfig config;
			config.nbody.enabled = numerical;
			config.evaluateOnly = true;
			config.useWeights = true;
			config.logDiagnostics = false;
			config.logCacheStats = false;
			config.outputEpochTdb = fitEpoch + 5.0;
			OrbitOptimizer optimizer(initial, oc::fit::prepareFitObservations(loaded.all, config), config);
			saved = optimizer.optimize();
			const auto& repeated = optimizer.optimize();
			ASSERT_TRUE(saved.orbit);
			ASSERT_TRUE(saved.fitEpochOrbit);
			EXPECT_EQ(saved.iterations, 0);
			EXPECT_FALSE(saved.converged);
			EXPECT_FALSE(saved.cancelled);
			EXPECT_DOUBLE_EQ(saved.orbit->getEpoch(), config.outputEpochTdb);
			EXPECT_DOUBLE_EQ(saved.fitEpochOrbit->getEpoch(), fitEpoch);
			EXPECT_DOUBLE_EQ(saved.fitEpochOrbit->getR().x, originalPosition.x);
			ASSERT_TRUE(saved.uncertainty.hasElementCovariance);
			EXPECT_TRUE(saved.uncertainty.cartesianStateCovariance.allFinite());
			EXPECT_TRUE(saved.uncertainty.cartesianStateCovariance.isApprox(
			    repeated.uncertainty.cartesianStateCovariance, 0.0));
			EXPECT_DOUBLE_EQ(saved.cost, repeated.cost);
			EXPECT_NE(saved.fitEpochOrbit.get(), repeated.fitEpochOrbit.get());
		}
		EXPECT_EQ(saved.residuals.count(), loaded.all.size());
		EXPECT_DOUBLE_EQ(saved.orbit->getEpoch(), fitEpoch + 5.0);
		EXPECT_DOUBLE_EQ(initial->getR().x, originalPosition.x);
	}
}

TEST_F(FitSessionTest, OptimizerRejectsMissingInitialOrbit)
{
	EXPECT_THROW(OrbitOptimizer(nullptr, {}), std::invalid_argument);
}

TEST(ObservationSelectionTest, IdsSurviveCopiesAndOverridesFollowReorderedRows)
{
	std::vector<ObsNode> observations(6);
	const auto includedId = observations[1].id;
	const auto excludedId = observations[4].id;
	ASSERT_NE(includedId, excludedId);
	auto copy = observations;
	std::ranges::reverse(copy);
	const std::vector<oc::ObservationOverride> overrides{{includedId, oc::ObservationDisposition::Include},
	                                                     {excludedId, oc::ObservationDisposition::Exclude}};
	const auto selection = oc::fit::resolveObservationSelection(copy, {}, overrides);
	EXPECT_EQ(selection.included, (std::vector<std::uint8_t>{1, 0, 1, 1, 1, 1}));
	EXPECT_EQ(selection.locked, (std::vector<std::uint8_t>{0, 1, 0, 0, 1, 0}));
	EXPECT_EQ(observations[1].id, includedId);
	auto invalid = overrides;
	invalid.push_back(overrides.front());
	EXPECT_THROW((void)oc::fit::resolveObservationSelection(copy, {}, invalid), std::invalid_argument);
	invalid = {{"missing", oc::ObservationDisposition::Exclude}};
	EXPECT_THROW((void)oc::fit::resolveObservationSelection(copy, {}, invalid), std::invalid_argument);
	copy[0].id = copy[1].id;
	EXPECT_THROW((void)oc::fit::resolveObservationSelection(copy, {}, overrides), std::invalid_argument);
}

TEST(ObservationSelectionTest, RejectionAndRecoveryRespectManualChoices)
{
	std::vector<ObsNode> observations(8);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000.0 + i;
	const std::vector<oc::ObservationOverride> overrides{{observations[0].id, oc::ObservationDisposition::Exclude},
	                                                     {observations[4].id, oc::ObservationDisposition::Include}};
	const std::array<std::uint8_t, 8> initial{1, 1, 1, 0, 1, 1, 1, 1};
	auto selection = oc::fit::resolveObservationSelection(observations, initial, overrides);
	const std::array<double, 8> chi{0, 0, 100, 0, 10000, 0, 0, 0};
	const auto rejected =
	    Calculations::applyHystereticReject(observations, selection.included, chi, {}, selection.locked);
	EXPECT_EQ(rejected.status, Calculations::RejectStatus::Changed);
	EXPECT_EQ(selection.included, (std::vector<std::uint8_t>{0, 1, 0, 1, 1, 1, 1, 1}));
	auto automatic = oc::fit::resolveObservationSelection(observations, selection.included, {});
	(void)Calculations::applyHystereticReject(observations, automatic.included, chi);
	EXPECT_EQ(automatic.included[0], 1);
	EXPECT_EQ(automatic.included[4], 0);
}

TEST(ObservationSelectionTest, ReplacedRecordsStayExcludedUnlessManuallyIncluded)
{
	std::vector<ObsNode> observations(8);
	for (std::size_t i = 0; i < observations.size(); ++i)
		observations[i].JD = 2459000.0 + i;
	observations[2].note2 = 'X';
	observations[5].note2 = 'x';
	auto selection = oc::fit::resolveObservationSelection(observations, {}, {});
	EXPECT_EQ(selection.included, (std::vector<std::uint8_t>{1, 1, 0, 1, 1, 0, 1, 1}));
	EXPECT_EQ(selection.locked, (std::vector<std::uint8_t>{0, 0, 1, 0, 0, 1, 0, 0}));
	const std::array<double, 8> chi{};
	const auto recovered =
	    Calculations::applyHystereticReject(observations, selection.included, chi, {}, selection.locked);
	EXPECT_EQ(recovered.status, Calculations::RejectStatus::NoChange);
	EXPECT_EQ(selection.included[2], 0);
	EXPECT_EQ(selection.included[5], 0);

	const std::array<std::uint8_t, 8> initial{1, 1, 1, 1, 1, 1, 1, 1};
	const std::vector<oc::ObservationOverride> overrides{{observations[2].id, oc::ObservationDisposition::Include},
	                                                     {observations[5].id, oc::ObservationDisposition::Auto}};
	selection = oc::fit::resolveObservationSelection(observations, initial, overrides);
	EXPECT_EQ(selection.included, (std::vector<std::uint8_t>{1, 1, 1, 1, 1, 0, 1, 1}));
	EXPECT_EQ(selection.locked, (std::vector<std::uint8_t>{0, 0, 1, 0, 0, 1, 0, 0}));
}

TEST_F(FitSessionTest, LoadedIdentitiesMatchAcrossPools)
{
	const auto data = fitContractCase();
	for (const auto& o : data.observations.iod) {
		const auto found = std::ranges::find(data.observations.all, o.id, &ObsNode::id);
		ASSERT_NE(found, data.observations.all.end());
		EXPECT_EQ(found->JD, o.JD);
		EXPECT_EQ(found->submittedRaRad(), o.submittedRaRad());
	}
}

TEST_F(FitSessionTest, RefinementOwnsSelectionAndDoesNotReloadOrRunIod)
{
	auto data = fitContractCase();
	oc::FitPipelineResult previous;
	previous.observations = data.observations;
	previous.observations.source = "this-file-does-not-exist.xml";
	previous.effectiveConfig = data.config;
	previous.fit = oc::FitSession{data.config}.optimize({data.initial, data.observations.all});
	ASSERT_TRUE(previous.fit.orbit);
	const auto oldState = previous.fit.orbit->getR();
	const auto oldMask = previous.fit.keptObsIndices;
	const auto oldCov = previous.fit.uncertainty.fitParameterCovariance;
	data.config.sigmaClip.enabled = true;
	const std::vector<oc::ObservationOverride> overrides{
	    {data.observations.all[1].id, oc::ObservationDisposition::Exclude},
	    {data.observations.all[2].id, oc::ObservationDisposition::Include}};
	const oc::SolutionRefinementRequest request{
	    previous, data.config, {.debias = true, .weightProfile = data.config.weightProfile}, overrides};
	const auto refined = oc::refineSolution(request, {.progress = [](const std::string& message) {
		                                        EXPECT_EQ(message.find("IOD"), std::string::npos);
		                                        EXPECT_EQ(message.find("ATT-EXT"), std::string::npos);
	                                        }});
	ASSERT_TRUE(refined.fit.orbit);
	EXPECT_EQ(refined.observationOverrides, overrides);
	EXPECT_EQ(refined.observations.all.size(), previous.observations.all.size());
	EXPECT_EQ(refined.observations.all[1].id, previous.observations.all[1].id);
	EXPECT_EQ(std::ranges::count(refined.fit.keptObsIndices, 1u), 0);
	EXPECT_EQ(std::ranges::count(refined.fit.keptObsIndices, 2u), 1);
	EXPECT_EQ(refined.fit.allResiduals.count(), previous.observations.all.size());
	EXPECT_EQ(refined.fit.uncertainty.fitParameterCovariance.rows(), 6);
	EXPECT_EQ(previous.fit.keptObsIndices, oldMask);
	EXPECT_EQ(previous.fit.orbit->getR().x, oldState.x);
	EXPECT_TRUE((previous.fit.uncertainty.fitParameterCovariance.array() == oldCov.array()).all());
	std::stop_source stop;
	const auto cancelled = oc::refineSolution(
	    request, {.cancel = stop.get_token(), .progress = [&](const auto&) { stop.request_stop(); }});
	EXPECT_TRUE(cancelled.fit.cancelled);
	EXPECT_FALSE(cancelled.fit.orbit);
	EXPECT_EQ(previous.fit.keptObsIndices, oldMask);
	auto invalid = overrides;
	invalid.front().id = "unknown";
	EXPECT_THROW((void)oc::refineSolution({previous, data.config, request.preparation, invalid}),
	             std::invalid_argument);
}

TEST_F(FitSessionTest, RefinementCanEnableAndDisableNonGravWithManualExclusion)
{
	auto data = fitContractCase();
	oc::FitPipelineResult previous;
	previous.observations = data.observations;
	previous.effectiveConfig = data.config;
	previous.fit = oc::FitSession{data.config}.optimize({data.initial, data.observations.all});
	const std::vector<oc::ObservationOverride> overrides{
	    {data.observations.all[1].id, oc::ObservationDisposition::Exclude}};
	const oc::ObservationPreparationOptions preparation{.debias = true, .weightProfile = data.config.weightProfile};
	data.config.nonGrav.enabled = true;
	data.config.nonGrav.model = oc::NonGravModelKind::Yarkovsky;
	auto enabled = oc::refineSolution({previous, data.config, preparation, overrides});
	ASSERT_TRUE(enabled.fit.orbit);
	EXPECT_TRUE(enabled.fit.finalNonGrav.enabled());
	EXPECT_EQ(enabled.fit.finalNonGrav.model(), NgModel::Yarkovsky);
	EXPECT_EQ(enabled.fit.uncertainty.fitParameterCovariance.rows(), 7);
	EXPECT_EQ(std::ranges::count(enabled.fit.keptObsIndices, 1u), 0);
	data.config.nonGrav.enabled = false;
	const auto disabled = oc::refineSolution({enabled, data.config, preparation, overrides});
	ASSERT_TRUE(disabled.fit.orbit);
	EXPECT_FALSE(disabled.fit.finalNonGrav.enabled());
	EXPECT_EQ(disabled.fit.uncertainty.fitParameterCovariance.rows(), 6);
	EXPECT_EQ(std::ranges::count(disabled.fit.keptObsIndices, 1u), 0);
	EXPECT_TRUE(enabled.fit.finalNonGrav.enabled());

	data.config.nonGrav.enabled = true;
	data.config.nonGrav.model = oc::NonGravModelKind::MarsdenH2O;
	data.config.sigmaClip.enabled = true;
	const auto marsden = oc::refineSolution({disabled, data.config, preparation, overrides});
	ASSERT_TRUE(marsden.fit.orbit);
	EXPECT_EQ(marsden.fit.finalNonGrav.model(), NgModel::MarsdenH2O);
	EXPECT_EQ(marsden.fit.uncertainty.fitParameterCovariance.rows(), 10);
	EXPECT_EQ(std::ranges::count(marsden.fit.keptObsIndices, 1u), 0);
	data.config.nonGrav.solveMask = static_cast<std::uint8_t>(NgParameterMask::A1 | NgParameterMask::A2);
	data.config.nonGrav.useInitialValues = true;
	data.config.nonGrav.initialA3 = 3e-10;
	const auto changed = oc::refineSolution({marsden, data.config, preparation, overrides});
	ASSERT_TRUE(changed.fit.orbit);
	EXPECT_EQ(changed.fit.finalNonGrav.fitParameterMask(), NgParameterMask::A1 | NgParameterMask::A2);
	EXPECT_EQ(changed.fit.finalNonGrav.A3(), 3e-10);
	EXPECT_EQ(changed.fit.uncertainty.fitParameterCovariance.rows(), 8);
}

TEST_F(FitSessionTest, NonGravWarmStartCannotOverrideRequestedModelOrSolveMask)
{
	auto data = fitContractCase();
	NonGravModel previous;
	previous.setEnabled(true);
	previous.setModel(NgModel::MarsdenH2O);
	previous.setParams(1e-10, 2e-10, 3e-10);
	data.config.nonGrav.enabled = true;
	data.config.nonGrav.model = oc::NonGravModelKind::MarsdenH2O;
	data.config.nonGrav.solveMask = static_cast<std::uint8_t>(NgParameterMask::A1 | NgParameterMask::A2);
	oc::FitRequest request{.initialOrbit = data.initial,
	                       .observations = data.observations.all,
	                       .warmStartNonGrav = &previous,
	                       .options = {.evaluateOnly = true}};
	auto result = oc::FitSession{data.config}.optimize(request);
	ASSERT_TRUE(result.orbit);
	EXPECT_EQ(result.finalNonGrav.fitParameterMask(), NgParameterMask::A1 | NgParameterMask::A2);
	EXPECT_EQ(result.finalNonGrav.A1(), previous.A1());
	EXPECT_EQ(result.finalNonGrav.A2(), previous.A2());
	EXPECT_EQ(result.finalNonGrav.A3(), previous.A3());
	EXPECT_EQ(result.uncertainty.fitParameterCovariance.rows(), 8);
	data.config.nonGrav.model = oc::NonGravModelKind::Yarkovsky;
	data.config.nonGrav.useInitialValues = true;
	data.config.nonGrav.initialA2 = -5e-14;
	result = oc::FitSession{data.config}.optimize(request);
	EXPECT_EQ(result.finalNonGrav.model(), NgModel::Yarkovsky);
	EXPECT_EQ(result.finalNonGrav.A2(), -5e-14);
}

TEST_F(FitSessionTest, RefinementRepreparesFromSubmittedCoordinatesWithoutAccumulatingBias)
{
	auto data = fitContractCase();
	oc::FitPipelineResult previous;
	previous.observations = data.observations;
	previous.effectiveConfig = data.config;
	previous.fit = oc::FitSession{data.config}.optimize({data.initial, data.observations.all});
	const auto off =
	    oc::refineSolution({previous, data.config, {.debias = false, .weightProfile = data.config.weightProfile}});
	for (std::size_t i = 0; i < off.observations.all.size(); ++i) {
		EXPECT_EQ(off.observations.all[i].RA, previous.observations.all[i].submittedRaRad());
		EXPECT_EQ(off.observations.all[i].id, previous.observations.all[i].id);
		EXPECT_EQ(off.observations.all[i].biasRaArcsec, 0.0);
	}
	const auto on =
	    oc::refineSolution({off, data.config, {.debias = true, .weightProfile = data.config.weightProfile}});
	for (std::size_t i = 0; i < on.observations.all.size(); ++i) {
		EXPECT_EQ(on.observations.all[i].RA, previous.observations.all[i].RA);
		EXPECT_EQ(on.observations.all[i].DEC, previous.observations.all[i].DEC);
	}
}
