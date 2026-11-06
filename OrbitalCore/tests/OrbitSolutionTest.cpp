#include "orbital_core/OrbitSolution.h"
#include "orbital_core/Exporter.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/TimeScale.h"
#ifdef ORBITAL_CORE_HAZARD_AVAILABLE
#include "orbital_core/hazard/IntervalMonteCarlo.h"
#endif
#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {
	using namespace orbital_core;

	class OrbitSolutionTest : public testing::Test {
	   protected:
		std::filesystem::path m_root;
		FitResult m_fit;

		void SetUp() override
		{
			m_root =
			    std::filesystem::temp_directory_path() /
			    ("orbit-solution-test-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
			ASSERT_TRUE(std::filesystem::create_directory(m_root));
			m_fit.fitEpochOrbit =
			    make_orbit_from_states_vectors(StateVectors(Vec3D(.9, .1, .05), Vec3D(-.1, 1.05, .02)), 2460000.5);
			ASSERT_TRUE(m_fit.fitEpochOrbit);
			m_fit.fitEpochOrbit->setEpochMjd(std::nextafter(60000.125, 60001.0));
			m_fit.orbit = m_fit.fitEpochOrbit;
			m_fit.converged = true;
			m_fit.uncertainty.cartesianStateCovariance = StateCov6x6::Identity() * 3e-16;
			m_fit.uncertainty.cartesianStateCovariance(0, 3) = 2e-17;
			m_fit.uncertainty.cartesianStateCovariance(3, 0) = 2e-17;
			m_fit.uncertainty.cartesianStateCovarianceFormal = m_fit.uncertainty.cartesianStateCovariance / 2.75;
			m_fit.uncertainty.covarianceObjectiveScale = 2.75;
		}

		void TearDown() override
		{
			if (std::filesystem::canonical(m_root).parent_path() ==
			        std::filesystem::canonical(std::filesystem::temp_directory_path()) &&
			    m_root.filename().string().starts_with("orbit-solution-test-"))
				std::filesystem::remove_all(m_root);
		}

		std::filesystem::path save(bool full = false)
		{
			ExportContext context;
			context.objectName = "Test Orbit";
			context.fitResult = &m_fit;
			context.rootDirectory = m_root;
			return full ? Exporter::save(context, {ExportProduct::StateVectors, ExportProduct::Covariance,
			                                       ExportProduct::FitSummary})
			                  .directory
			            : Exporter::save(context, {ExportProduct::StateVectors}).directory;
		}

		void replace(const std::filesystem::path& file, const std::string& from, const std::string& to)
		{
			std::ifstream input(file, std::ios::binary);
			std::string contents(std::istreambuf_iterator<char>(input), {});
			input.close();
			const auto pos = contents.find(from);
			ASSERT_NE(pos, std::string::npos);
			contents.replace(pos, from.size(), to);
			std::ofstream(file, std::ios::binary) << contents;
		}

		void nonGrav(NgModel model, bool fitDt = false)
		{
			m_fit.finalNonGrav.setModel(model);
			m_fit.finalNonGrav.setEnabled(true);
			if (fitDt)
				m_fit.finalNonGrav.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2 |
				                                         NgParameterMask::A3 | NgParameterMask::DT);
			m_fit.finalNonGrav.setParams(model == NgModel::Yarkovsky ? 0 : 1e-10, -2e-12,
			                             model == NgModel::Yarkovsky ? 0 : 3e-11);
			m_fit.finalNonGrav.setDT(fitDt ? .5 : 0);
			const int dimension = 6 + m_fit.finalNonGrav.fitParameterCount();
			m_fit.uncertainty.fitParameterCovariance = FitParameterCovariance::Identity(dimension, dimension) * 1e-24;
			m_fit.uncertainty.fitParameterCovariance.topLeftCorner<6, 6>() = m_fit.uncertainty.cartesianStateCovariance;
			m_fit.uncertainty.fitParameterCovariance(0, 6) = m_fit.uncertainty.fitParameterCovariance(6, 0) = 1e-22;
			m_fit.uncertainty.fitParameterCovarianceFormal = m_fit.uncertainty.fitParameterCovariance / 2.75;
		}
	};

	TEST_F(OrbitSolutionTest, LoadsExistingCarWithoutFitOrRescaling)
	{
		const auto directory = save();
		const auto solution = loadOrbitSolution(directory);
		EXPECT_EQ(solution.objectName, "Test Orbit");
		EXPECT_DOUBLE_EQ(solution.nominal.epochMjd(), m_fit.fitEpochOrbit->epochMjd());
		EXPECT_NE(solution.nominal.epochMjd(), time::mjdFromJd(solution.nominal.getEpoch()));
		EXPECT_DOUBLE_EQ(solution.nominal.getR().x, m_fit.fitEpochOrbit->getR().x);
		EXPECT_LT((solution.nominal.getV() - m_fit.fitEpochOrbit->getV()).magnitude(), 1e-15);
		EXPECT_TRUE(solution.covariance.isApprox(m_fit.uncertainty.cartesianStateCovariance, 1e-15));
		EXPECT_FALSE(solution.nonGrav.enabled());
		EXPECT_FALSE(solution.absoluteMagnitude.has_value());
	}

	TEST_F(OrbitSolutionTest, CometExportDoesNotMislabelItsMagnitudeAsAsteroidH)
	{
		Residual row{};
		row.mag = 15.0;
		row.magBand = 'T';
		row.nuclearMagnitude = false;
		row.magSigma = 0.5;
		row.helioDistAu = row.geoDistAu = 1.0;
		m_fit.residuals.add(row);
		ExportContext context;
		context.objectName = "24P";
		context.fitResult = &m_fit;
		context.rootDirectory = m_root;
		context.cometPhotometry = true;
		const auto directory =
		    Exporter::save(context, {ExportProduct::StateVectors, ExportProduct::FitSummary}).directory;
		const auto solution = loadOrbitSolution(directory);
		EXPECT_FALSE(solution.absoluteMagnitude);
		std::ifstream oef(directory / "orbit_cartesian_oef.txt");
		const std::string content(std::istreambuf_iterator<char>(oef), {});
		EXPECT_NE(content.find(" MAG  99.999  0.150"), std::string::npos);
		std::ifstream summary(directory / "fit_summary.txt");
		const std::string details(std::istreambuf_iterator<char>(summary), {});
		EXPECT_NE(details.find("M1 15"), std::string::npos);
		EXPECT_NE(details.find("K1_fixed 10"), std::string::npos);
	}

	TEST_F(OrbitSolutionTest, ExportedAsteroidGMatchesTheMagnitudeModel)
	{
		Residual row{};
		row.mag = 20.0;
		row.magBand = 'V';
		row.magSigma = 0.5;
		row.helioDistAu = row.geoDistAu = 1.0;
		row.phaseAngleRad = 0.6;
		m_fit.residuals.add(row);
		ExportContext context;
		context.objectName = "Test Orbit";
		context.fitResult = &m_fit;
		context.rootDirectory = m_root;
		context.photometry.asteroidG = 0.42;
		const auto directory =
		    Exporter::save(context, {ExportProduct::StateVectors, ExportProduct::FitSummary}).directory;
		const auto solution = loadOrbitSolution(directory);
		const auto expected = photometry::fit(m_fit.residuals.data(), 0.42);
		ASSERT_TRUE(solution.absoluteMagnitude);
		ASSERT_TRUE(expected.absoluteMagnitude);
		EXPECT_DOUBLE_EQ(solution.slopeG, 0.42);
		EXPECT_NEAR(*solution.absoluteMagnitude, *expected.absoluteMagnitude, 0.00051);
	}

	TEST_F(OrbitSolutionTest, TtAndTdbEpochsReferToTheSameState)
	{
		const auto directory = save();
		const auto tdb = loadOrbitSolution(directory);
		std::ostringstream source, target;
		source << std::fixed << std::setprecision(12) << std::setw(23) << m_fit.fitEpochOrbit->epochMjd() << " TDB";
		target << std::fixed << std::setprecision(12) << std::setw(23)
		       << time::fromMjdTdb(m_fit.fitEpochOrbit->epochMjd(), time::Scale::TT) << " TDT";
		replace(directory / "orbit_cartesian_oef.txt", source.str(), target.str());
		const auto tt = loadOrbitSolution(directory);
		EXPECT_DOUBLE_EQ(tt.nominal.epochMjd(), tdb.nominal.epochMjd());
		EXPECT_DOUBLE_EQ(tt.nominal.getR().x, tdb.nominal.getR().x);
		EXPECT_DOUBLE_EQ(tt.nominal.getV().x, tdb.nominal.getV().x);
		EXPECT_TRUE(tt.covariance.isApprox(tdb.covariance, 1e-15));
	}

	TEST_F(OrbitSolutionTest, LoadsFullNonGravCovarianceAndCrossTerms)
	{
		for (const auto model : {NgModel::Yarkovsky, NgModel::MarsdenH2O, NgModel::MarsdenCO, NgModel::InverseSquare}) {
			nonGrav(model);
			const auto solution = loadOrbitSolution(save(true));
			EXPECT_EQ(solution.nonGrav.model(), model);
			EXPECT_EQ(solution.nonGrav.fitParameterMask(), m_fit.finalNonGrav.fitParameterMask());
			EXPECT_DOUBLE_EQ(solution.nonGrav.A2(), -2e-12);
			EXPECT_DOUBLE_EQ(solution.nonGrav.DT(), m_fit.finalNonGrav.DT());
			EXPECT_TRUE(solution.covariance.isApprox(m_fit.uncertainty.fitParameterCovariance, 1e-14));
		}
	}

	TEST_F(OrbitSolutionTest, MigratesLegacyCoCoefficientsAndCovariance)
	{
		nonGrav(NgModel::MarsdenCO);
		m_fit.finalNonGrav.setDT(2.0);
		const auto directory = save();
		replace(directory / "orbit_cartesian_oef.txt", "! OC_NG 3 ", "! OC_NG 1 ");
		const auto solution = loadOrbitSolution(directory);
		constexpr double scale = 0.01003 / 0.0408373333128795;
		EXPECT_DOUBLE_EQ(solution.nonGrav.A1(), m_fit.finalNonGrav.A1() * scale);
		EXPECT_DOUBLE_EQ(solution.nonGrav.A2(), m_fit.finalNonGrav.A2() * scale);
		EXPECT_DOUBLE_EQ(solution.nonGrav.A3(), m_fit.finalNonGrav.A3() * scale);
		EXPECT_DOUBLE_EQ(solution.nonGrav.DT(), 2.0);
		Eigen::MatrixXd conversion = Eigen::MatrixXd::Identity(10, 10);
		conversion(6, 6) = conversion(7, 7) = conversion(8, 8) = scale;
		const Eigen::MatrixXd expected = conversion * m_fit.uncertainty.fitParameterCovariance * conversion;
		EXPECT_TRUE(solution.covariance.isApprox(expected, 1e-14));
		for (double r : {0.5, 1.0, 5.0})
			EXPECT_NEAR(solution.nonGrav.A1() * solution.nonGrav.gr(r),
			            m_fit.finalNonGrav.A1() * scale * m_fit.finalNonGrav.gr(r), 1e-24);
	}

	TEST_F(OrbitSolutionTest, PreservesInverseSquareWithAndWithoutFittedDelay)
	{
		nonGrav(NgModel::InverseSquare, true);
		const auto withDelay = save();
		const auto delayed = loadOrbitSolution(withDelay).nonGrav;
		EXPECT_EQ(delayed.model(), NgModel::InverseSquare);
		EXPECT_DOUBLE_EQ(delayed.DT(), .5);
		EXPECT_TRUE(hasNgParameter(delayed.fitParameterMask(), NgParameterMask::DT));

		nonGrav(NgModel::InverseSquare);
		m_fit.finalNonGrav.setDT(.5);
		const auto withoutDelay = save();
		const auto original = loadOrbitSolution(withoutDelay).nonGrav;
		EXPECT_EQ(original.model(), NgModel::InverseSquare);
		EXPECT_DOUBLE_EQ(original.DT(), .5);
		EXPECT_FALSE(hasNgParameter(original.fitParameterMask(), NgParameterMask::DT));
	}

	TEST_F(OrbitSolutionTest, DoesNotSilentlyDropNonGrav)
	{
		nonGrav(NgModel::Yarkovsky);
		const auto solution = loadOrbitSolution(save());
		EXPECT_TRUE(solution.nonGrav.enabled());
		EXPECT_EQ(solution.covariance.rows(), 7);
	}

	TEST_F(OrbitSolutionTest, StandaloneOefIgnoresSidecarsAndRejectsMissingNgr)
	{
		nonGrav(NgModel::Yarkovsky);
		const auto directory = save(true);
		replace(directory / "fit_summary.txt", "model 4", "model 3");
		EXPECT_NO_THROW((void)loadOrbitSolution(directory));
		replace(directory / "orbit_cartesian_oef.txt", " NGR", "! NGR");
		EXPECT_THROW((void)loadOrbitSolution(directory), std::invalid_argument);
	}

	TEST_F(OrbitSolutionTest, CoreExtensionCannotBeMistakenForStandardOef)
	{
		for (const auto model : {NgModel::InverseSquare, NgModel::MarsdenCO}) {
			nonGrav(model);
			const auto directory = save();
			EXPECT_EQ(loadOrbitSolution(directory).nonGrav.model(), model);
			replace(directory / "orbit_cartesian_oef.txt", "'OC-OEF1.0'", "'OEF2.0'");
			EXPECT_THROW((void)loadOrbitSolution(directory), std::invalid_argument);
		}
	}

	TEST_F(OrbitSolutionTest, PreservesDisabledModelCoefficientsAndConfiguredMask)
	{
		nonGrav(NgModel::InverseSquare, true);
		m_fit.finalNonGrav.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2);
		m_fit.finalNonGrav.setEnabled(false);
		m_fit.uncertainty.fitParameterCovariance = m_fit.uncertainty.cartesianStateCovariance;
		auto restored = loadOrbitSolution(save()).nonGrav;
		EXPECT_FALSE(restored.enabled());
		EXPECT_EQ(restored.model(), NgModel::InverseSquare);
		EXPECT_DOUBLE_EQ(restored.A1(), m_fit.finalNonGrav.A1());
		EXPECT_DOUBLE_EQ(restored.A2(), m_fit.finalNonGrav.A2());
		EXPECT_DOUBLE_EQ(restored.DT(), .5);
		restored.setEnabled(true);
		EXPECT_EQ(restored.fitParameterMask(), NgParameterMask::A1 | NgParameterMask::A2);
	}

	TEST_F(OrbitSolutionTest, LoadsEachExportedCoordinateChart)
	{
		nonGrav(NgModel::Yarkovsky);
		ExportContext context;
		context.objectName = "2020VV";
		context.fitResult = &m_fit;
		context.rootDirectory = m_root;
		const auto directory = Exporter::save(context, {ExportProduct::OrbitElements}).directory;
		for (const char* file : {"orbit_equinoctial_oef.txt", "orbit_keplerian_oef.txt"}) {
			const auto solution = loadOrbitSolution(directory / file);
			EXPECT_LT((solution.nominal.getR() - m_fit.fitEpochOrbit->getR()).magnitude(), 1e-13);
			const auto sigma = m_fit.uncertainty.fitParameterCovariance.diagonal().cwiseSqrt().eval();
			EXPECT_LT((solution.covariance - m_fit.uncertainty.fitParameterCovariance)
			              .cwiseQuotient(sigma * sigma.transpose())
			              .cwiseAbs()
			              .maxCoeff(),
			          1e-6);
			EXPECT_DOUBLE_EQ(solution.nonGrav.A2(), m_fit.finalNonGrav.A2());
		}
	}

	TEST_F(OrbitSolutionTest, PreservesFixedNonGravWithSixDimensionalCovariance)
	{
		nonGrav(NgModel::MarsdenH2O);
		m_fit.finalNonGrav.setDT(-12.3);
		m_fit.finalNonGrav.setFitParameterMask(NgParameterMask::None);
		m_fit.uncertainty.fitParameterCovariance = m_fit.uncertainty.cartesianStateCovariance;
		const auto solution = loadOrbitSolution(save());
		EXPECT_TRUE(solution.nonGrav.enabled());
		EXPECT_EQ(solution.nonGrav.numFitParams(), 0);
		EXPECT_DOUBLE_EQ(solution.nonGrav.DT(), -12.3);
		EXPECT_DOUBLE_EQ(solution.nonGrav.A1(), 1e-10);
	}

	TEST_F(OrbitSolutionTest, RejectsWrongFrameTimeAndRecords)
	{
		for (const auto& [from, to] :
		     std::vector<std::pair<std::string, std::string>>{{"ECLM J2000", "EQUM J2000"},
		                                                      {" TDB", " UTC"},
		                                                      {" CAR", " COM"},
		                                                      {" LSP   0  0    6", " LSP   0  0    7"},
		                                                      {" COV", " BAD"},
		                                                      {" MAG", " MJD 60000 TDT\n MAG"},
		                                                      {" MAG", " MAG 0 0.15\n MAG"},
		                                                      {"format  =", "format  !"}}) {
			const auto directory = save();
			replace(directory / "orbit_cartesian_oef.txt", from, to);
			EXPECT_THROW((void)loadOrbitSolution(directory), std::invalid_argument) << from;
		}
	}

	TEST_F(OrbitSolutionTest, RejectsIndefiniteOrMissingCovariance)
	{
		m_fit.uncertainty.cartesianStateCovariance(0, 1) = 9e-16;
		m_fit.uncertainty.cartesianStateCovariance(1, 0) = 9e-16;
		EXPECT_THROW((void)loadOrbitSolution(save()), std::invalid_argument);
		m_fit.uncertainty.cartesianStateCovariance.setZero();
		EXPECT_THROW((void)loadOrbitSolution(save()), std::invalid_argument);
	}

#ifdef ORBITAL_CORE_HAZARD_AVAILABLE
	TEST_F(OrbitSolutionTest, CancellationPreventsPropagationAndCallbacks)
	{
		const auto solution = loadOrbitSolution(save());
		NBodySettings forces;
		forces.enabled = true;
		hazard::IntervalMcConfig config;
		config.startMjdTdb = solution.nominal.epochMjd();
		config.endMjdTdb = config.startMjdTdb + 1;
		std::stop_source stop;
		stop.request_stop();
		const auto result = hazard::runIntervalMonteCarlo(solution, forces, config, stop.get_token(), [](int, int) {
			ADD_FAILURE() << "Cancelled propagation called progress";
		});
		EXPECT_TRUE(result.cancelled);
		EXPECT_EQ(result.completed, 0);
		EXPECT_EQ(result.config.altitudeKm, 100);
		EXPECT_TRUE(result.events.empty());
	}

	TEST_F(OrbitSolutionTest, RejectsPropagationSettingsBeforeStarting)
	{
		auto solution = loadOrbitSolution(save());
		NBodySettings forces;
		forces.enabled = true;
		hazard::IntervalMcConfig config;
		config.startMjdTdb = solution.nominal.epochMjd();
		config.endMjdTdb = config.startMjdTdb + 1;
		config.workers = 0;
		EXPECT_THROW((void)hazard::runIntervalMonteCarlo(solution, forces, config), std::invalid_argument);
		config.workers = 1;
		config.integrator = static_cast<IntegratorKind>(99);
		EXPECT_THROW((void)hazard::runIntervalMonteCarlo(solution, forces, config), std::invalid_argument);
		config.integrator = IntegratorKind::RA15;
		solution.covariance.resize(7, 7);
		EXPECT_THROW((void)hazard::runIntervalMonteCarlo(solution, forces, config), std::invalid_argument);
	}
#endif
}  // namespace
