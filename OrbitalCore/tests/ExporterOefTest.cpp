#include <gtest/gtest.h>

#include "orbital_core/Exporter.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/Constants.h"
#include "orbital_core/StateVectors.h"
#include "hazard/LovCoordinateChart.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace orbital_core;

namespace {
	struct OefData {
		std::map<std::string, std::vector<double>> fields;

		Eigen::VectorXd vector(const std::string& tag) const
		{
			const auto& values = fields.at(tag);
			return Eigen::Map<const Eigen::VectorXd>(values.data(), values.size());
		}

		Eigen::MatrixXd matrix(const std::string& tag, int n = 6) const
		{
			const auto& values = fields.at(tag);
			if (values.size() != static_cast<std::size_t>(n * (n + 1) / 2))
				throw std::runtime_error("Unexpected OEF matrix dimension");
			Eigen::MatrixXd result(n, n);
			std::size_t k = 0;
			for (int i = 0; i < n; ++i) {
				for (int j = i; j < n; ++j)
					result(i, j) = result(j, i) = values[k++];
			}
			return result;
		}
	};

	OefData readOef(const std::filesystem::path& path)
	{
		std::ifstream input(path);
		if (!input)
			throw std::runtime_error("Cannot read OEF: " + path.string());
		OefData result;
		std::string line;
		while (std::getline(input, line)) {
			std::istringstream row(line);
			std::string tag;
			row >> tag;
			if (tag == "!")
				row >> tag;
			double value;
			while (row >> value)
				result.fields[tag].push_back(value);
		}
		return result;
	}

	std::string readFileBytes(const std::filesystem::path& path)
	{
		std::ifstream input(path, std::ios::binary);
		if (!input)
			throw std::runtime_error("Cannot read export: " + path.string());
		std::ostringstream contents;
		contents << input.rdbuf();
		return contents.str();
	}

	double standardizedInverseError(const OefData& data, int n = 6)
	{
		const Eigen::MatrixXd cov = data.matrix("COV", n);
		const Eigen::MatrixXd normal = data.matrix("NOR", n);
		const Eigen::VectorXd sigma = cov.diagonal().cwiseSqrt();
		const Eigen::MatrixXd scale = sigma * sigma.transpose();
		const Eigen::MatrixXd correlation = cov.cwiseQuotient(scale);
		const Eigen::MatrixXd scaledNormal = normal.cwiseProduct(scale);
		return (correlation * scaledNormal - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff();
	}

	class ExporterOefTest : public ::testing::Test {
	   protected:
		FitResult fit;
		std::filesystem::path root;
		bool ownsRoot = false;

		void SetUp() override
		{
			const auto input =
			    readOef(std::filesystem::path(__FILE__).parent_path() / "data/2020VV_oef_covariance.txt");
			const auto state = input.vector("CAR");
			fit.fitEpochOrbit = make_orbit_from_states_vectors(
			    StateVectors(Vec3D(state(0), state(1), state(2)),
			                 Vec3D(state(3), state(4), state(5)) * (1.0 / Constants::GAUSS_K)),
			    input.vector("MJD")(0) + 2400000.5);
			ASSERT_TRUE(fit.fitEpochOrbit);
			fit.orbit = fit.fitEpochOrbit;
			fit.uncertainty.cartesianStateCovariance = input.matrix("COV");
			const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
			root = std::filesystem::temp_directory_path() / ("orbital_oef_test_" + std::to_string(stamp));
			ownsRoot = std::filesystem::create_directory(root);
			ASSERT_TRUE(ownsRoot);
		}

		void TearDown() override
		{
			if (ownsRoot) {
				std::error_code ec;
				std::filesystem::remove_all(root, ec);
			}
		}

		ExportReport save(std::initializer_list<ExportProduct> products)
		{
			ExportContext context;
			context.objectName = "2020VV";
			context.fitResult = &fit;
			context.rootDirectory = root;
			return Exporter::save(context, products);
		}

		void setCometaryOrbit(double eccentricity, double trueAnomaly = 0.6)
		{
			cometary::Vector6 elements;
			elements << 0.8, eccentricity, 0.4, 0.7, 1.1, trueAnomaly;
			const auto state = cometary::toCartesian(elements);
			ASSERT_TRUE(state);
			const double epoch = fit.fitEpochOrbit->getEpoch();
			fit.fitEpochOrbit = make_orbit_from_states_vectors(
			    StateVectors(Vec3D((*state)(0), (*state)(1), (*state)(2)),
			                 Vec3D((*state)(3), (*state)(4), (*state)(5)) * (1.0 / Constants::GAUSS_K)),
			    epoch);
			ASSERT_TRUE(fit.fitEpochOrbit);
			fit.orbit = fit.fitEpochOrbit;
		}
	};
}  // namespace

TEST_F(ExporterOefTest, NormalInvertsPrinted2020VvCovariance)
{
	const auto original = fit.uncertainty.cartesianStateCovariance;
	const auto report = save({ExportProduct::StateVectors, ExportProduct::OrbitElements});
	const auto car = readOef(report.directory / "orbit_cartesian_oef.txt");
	const auto equ = readOef(report.directory / "orbit_equinoctial_oef.txt");
	EXPECT_LT(standardizedInverseError(car), 2e-7);
	EXPECT_LT(standardizedInverseError(equ), 2e-7);
	EXPECT_TRUE(car.matrix("COV").isApprox(original, 1e-15));
	EXPECT_EQ((fit.uncertainty.cartesianStateCovariance - original).norm(), 0.0);
	EXPECT_EQ(car.vector("MJD")(0), equ.vector("MJD")(0));
	EXPECT_FALSE(std::filesystem::exists(report.directory / "orbit_cometary_oef.txt"));
	EXPECT_TRUE(std::filesystem::exists(report.directory / "orbit_keplerian_oef.txt"));
	EXPECT_EQ(report.files.size(), 3u);
}

TEST_F(ExporterOefTest, EigenCommentsUseRadianBasisAndWeakAxisSign)
{
	const auto report = save({ExportProduct::OrbitElements});
	const auto equ = readOef(report.directory / "orbit_equinoctial_oef.txt");
	Eigen::MatrixXd radCov = equ.matrix("COV");
	radCov.row(5) *= Constants::DEG_TO_RAD;
	radCov.col(5) *= Constants::DEG_TO_RAD;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(radCov);
	ASSERT_EQ(solver.info(), Eigen::Success);
	const auto eig = equ.vector("EIG");
	for (int i = 0; i < 6; ++i)
		EXPECT_NEAR(eig(i), std::sqrt(solver.eigenvalues()(i)), 6e-6 * eig(i));
	EXPECT_NEAR(eig(5), 2.83631e-5, 1e-10);

	const auto weak = equ.vector("WEA");
	EXPECT_GE(weak(0), 0.0);
	Eigen::VectorXd expected = solver.eigenvectors().col(5);
	if (expected(0) < 0.0)
		expected = -expected;
	EXPECT_LT((weak - expected).cwiseAbs().maxCoeff(), 5.1e-6);
	EXPECT_NEAR(weak.norm(), 1.0, 1e-5);
	const Eigen::Vector<double, 6> referenceWeakDirection(0.45757, -0.11198, 0.40334, 0.00291, 0.00859, -0.78443);
	EXPECT_LT((weak - referenceWeakDirection).cwiseAbs().maxCoeff(), 1e-5);
}

TEST_F(ExporterOefTest, NormalUsesTheSameTriangleAsCovarianceOutput)
{
	auto& cov = fit.uncertainty.cartesianStateCovariance;
	cov.setIdentity();
	cov(0, 1) = 0.25;
	cov(1, 0) = 0.5;
	const auto report = save({ExportProduct::StateVectors});
	const auto car = readOef(report.directory / "orbit_cartesian_oef.txt");
	EXPECT_EQ(car.matrix("COV")(1, 0), 0.25);
	EXPECT_LT(standardizedInverseError(car), 1e-13);
}

TEST_F(ExporterOefTest, NormalHandlesMixedUnits)
{
	const Eigen::Vector<double, 6> sigma(1e-12, 1e-8, 1e-3, 1.0, 1e4, 1e8);
	Eigen::Matrix<double, 6, 6> corr = Eigen::Matrix<double, 6, 6>::Constant(0.2);
	corr.diagonal().setOnes();
	fit.uncertainty.cartesianStateCovariance = sigma.asDiagonal() * corr * sigma.asDiagonal();
	const auto report = save({ExportProduct::StateVectors});
	EXPECT_LT(standardizedInverseError(readOef(report.directory / "orbit_cartesian_oef.txt")), 1e-13);
}

TEST_F(ExporterOefTest, OmitsNormalForSingularOrIndefiniteCovariance)
{
	auto& cov = fit.uncertainty.cartesianStateCovariance;
	for (int variant = 0; variant < 3; ++variant) {
		SCOPED_TRACE(variant);
		cov.setIdentity();
		if (variant == 0)
			cov(0, 0) = 0.0;
		else if (variant == 1)
			cov(0, 0) = -1.0;
		else
			cov(0, 1) = cov(1, 0) = 2.0;
		const auto report = save({ExportProduct::StateVectors});
		const auto car = readOef(report.directory / "orbit_cartesian_oef.txt");
		EXPECT_TRUE(car.fields.contains("COV"));
		EXPECT_FALSE(car.fields.contains("NOR"));
	}
}

TEST_F(ExporterOefTest, EquinoctialNormalIncludesFittedNonGravParameter)
{
	fit.finalNonGrav.setModel(NgModel::Yarkovsky);
	fit.finalNonGrav.setEnabled(true);
	auto& cov = fit.uncertainty.fitParameterCovariance;
	cov = FitParameterCovariance::Zero(7, 7);
	cov.topLeftCorner(6, 6) = fit.uncertainty.cartesianStateCovariance;
	cov(6, 6) = 1e-26;
	for (int i = 0; i < 6; ++i)
		cov(i, 6) = cov(6, i) = 0.2 * cov(i, 0) / std::sqrt(cov(0, 0)) * 1e-13;
	const auto report = save({ExportProduct::OrbitElements});
	const auto equ = readOef(report.directory / "orbit_equinoctial_oef.txt");
	ASSERT_EQ(equ.vector("EIG").size(), 7);
	EXPECT_NEAR(equ.matrix("COV", 7)(6, 6), 1e-6, 1e-20);
	EXPECT_LT(standardizedInverseError(equ, 7), 2e-7);
}

TEST_F(ExporterOefTest, EquinoctialExportUsesTheLovCoordinatesAndCovariance)
{
	const auto chart =
	    hazard::makeEquinoctialLovChart(*fit.fitEpochOrbit, fit.uncertainty.cartesianStateCovariance, fit.finalNonGrav);
	ASSERT_TRUE(chart.valid);
	const auto report = save({ExportProduct::OrbitElements});
	const auto output = readOef(report.directory / "orbit_equinoctial_oef.txt");
	Eigen::VectorXd expected = chart.nominalCoordinates;
	expected(5) *= Constants::RAD_TO_DEG;
	EXPECT_LT((output.vector("EQU") - expected).norm(), 1.0e-13);
	Eigen::MatrixXd covariance = chart.coordinateCovariance;
	covariance.row(5) *= Constants::RAD_TO_DEG;
	covariance.col(5) *= Constants::RAD_TO_DEG;
	EXPECT_LT((output.matrix("COV") - covariance).norm() / covariance.norm(), 2.0e-15);
}

TEST_F(ExporterOefTest, CometaryExportPreservesConicsAcrossParabola)
{
	for (double e : {0.990001, 1.0 - 1e-8, 1.0, 1.0 + 1e-8, 1.2}) {
		for (double nu : {-0.6, 0.6}) {
			SCOPED_TRACE(e);
			SCOPED_TRACE(nu);
			ASSERT_NO_FATAL_FAILURE(setCometaryOrbit(e, nu));
			const auto original = fit.uncertainty.cartesianStateCovariance;
			const auto report = save({ExportProduct::OrbitElements, ExportProduct::StateVectors});
			const auto cot = readOef(report.directory / "orbit_cometary_oef.txt");
			const auto car = readOef(report.directory / "orbit_cartesian_oef.txt");
			ASSERT_EQ(cot.vector("COT").size(), 6);
			EXPECT_NEAR(cot.vector("COT")(0), 0.8, 3e-14);
			EXPECT_NEAR(cot.vector("COT")(1), e, 3e-14);
			const double expectedNu = nu < 0.0 ? nu + Constants::TWOPI : nu;
			EXPECT_NEAR(cot.vector("COT")(5), expectedNu * Constants::RAD_TO_DEG, 3e-12);
			EXPECT_TRUE(cot.matrix("COV").allFinite());
			EXPECT_TRUE(car.matrix("COV").isApprox(original, 1e-15));
			EXPECT_EQ(cot.vector("MJD")(0), car.vector("MJD")(0));
			EXPECT_FALSE(std::filesystem::exists(report.directory / "orbit_keplerian_oef.txt"));
			EXPECT_FALSE(std::filesystem::exists(report.directory / "orbit_equinoctial_oef.txt"));
			EXPECT_EQ(report.files.size(), 2u);
			EXPECT_EQ((fit.uncertainty.cartesianStateCovariance - original).norm(), 0.0);
		}
	}
}

TEST_F(ExporterOefTest, CometaryAllAndReorderedProductsReportCartesianExactlyOnce)
{
	ASSERT_NO_FATAL_FAILURE(setCometaryOrbit(1.2));
	auto checkFiles = [](const ExportReport& report, std::size_t expectedCount) {
		const std::set<std::filesystem::path> unique(report.files.begin(), report.files.end());
		EXPECT_EQ(report.files.size(), expectedCount);
		EXPECT_EQ(unique.size(), report.files.size());
		EXPECT_EQ(std::count(report.files.begin(), report.files.end(), report.directory / "orbit_cartesian_oef.txt"),
		          1);
		EXPECT_EQ(std::count(report.files.begin(), report.files.end(), report.directory / "orbit_cometary_oef.txt"), 1);
	};
	const auto all = save({ExportProduct::All});
	checkFiles(all, 3u);
	const std::string car = readFileBytes(all.directory / "orbit_cartesian_oef.txt");
	const auto reordered = save({ExportProduct::StateVectors, ExportProduct::OrbitElements});
	checkFiles(reordered, 2u);
	EXPECT_EQ(readFileBytes(reordered.directory / "orbit_cartesian_oef.txt"), car);
}

TEST_F(ExporterOefTest, ReusedDirectoryDoesNotReportPreviousKeplerianOrEquinoctialFiles)
{
	const auto previous = save({ExportProduct::StateVectors, ExportProduct::OrbitElements});
	const auto equFile = previous.directory / "orbit_equinoctial_oef.txt";
	const auto kepFile = previous.directory / "orbit_keplerian_oef.txt";
	const auto carFile = previous.directory / "orbit_cartesian_oef.txt";
	const std::string equ = readFileBytes(equFile);
	const std::string kep = readFileBytes(kepFile);
	const std::string car = readFileBytes(carFile);
	ASSERT_NO_FATAL_FAILURE(setCometaryOrbit(1.2));
	const auto current = save({ExportProduct::All});
	EXPECT_EQ(current.directory, previous.directory);
	EXPECT_EQ(std::count(current.files.begin(), current.files.end(), equFile), 0);
	EXPECT_EQ(std::count(current.files.begin(), current.files.end(), kepFile), 0);
	EXPECT_EQ(std::count(current.files.begin(), current.files.end(), carFile), 1);
	EXPECT_EQ(readFileBytes(equFile), equ);
	EXPECT_EQ(readFileBytes(kepFile), kep);
	EXPECT_NE(readFileBytes(carFile), car);
	EXPECT_TRUE(readOef(current.directory / "orbit_cometary_oef.txt").vector("COT").allFinite());
}

TEST_F(ExporterOefTest, CometaryCovarianceUsesAnalyticJacobianAndBothDegreeAxes)
{
	ASSERT_NO_FATAL_FAILURE(setCometaryOrbit(1.0));
	fit.orbit = std::make_shared<Orbit>(*fit.fitEpochOrbit);
	fit.orbit->setEpoch(fit.fitEpochOrbit->getEpoch() + 20.0);
	const Eigen::Vector<double, 6> sigma(1e-6, 2e-6, 3e-6, 1e-8, 2e-8, 3e-8);
	cometary::Matrix6 correlation = cometary::Matrix6::Constant(0.1);
	correlation.diagonal().setOnes();
	fit.uncertainty.cartesianStateCovariance = sigma.asDiagonal() * correlation * sigma.asDiagonal();
	const auto coordinates = hazard::cometaryCoordinates(*fit.fitEpochOrbit, fit.finalNonGrav);
	ASSERT_TRUE(coordinates);
	const auto jacobian = cometary::cartesianToElementsJacobian(coordinates->head<6>());
	ASSERT_TRUE(jacobian);
	cometary::Matrix6 expected = *jacobian * fit.uncertainty.cartesianStateCovariance * jacobian->transpose();
	expected = (0.5 * (expected + expected.transpose())).eval();
	for (int i = 2; i < 6; ++i) {
		expected.row(i) *= Constants::RAD_TO_DEG;
		expected.col(i) *= Constants::RAD_TO_DEG;
	}
	const auto report = save({ExportProduct::OrbitElements});
	const auto output = readOef(report.directory / "orbit_cometary_oef.txt");
	EXPECT_LT((output.matrix("COV") - expected).norm() / expected.norm(), 2e-14);
	EXPECT_LT(standardizedInverseError(output), 2e-7);
	EXPECT_EQ(output.vector("MJD")(0), fit.fitEpochOrbit->epochMjd());
	EXPECT_NE(output.vector("MJD")(0), fit.orbit->epochMjd());
	EXPECT_EQ(report.files.size(), 2u);
}

TEST_F(ExporterOefTest, CometaryCovariancePreservesNonGravCrossTerms)
{
	ASSERT_NO_FATAL_FAILURE(setCometaryOrbit(1.2));
	fit.finalNonGrav.setModel(NgModel::Yarkovsky);
	fit.finalNonGrav.setEnabled(true);
	auto& cov = fit.uncertainty.fitParameterCovariance;
	cov = FitParameterCovariance::Zero(7, 7);
	cov.topLeftCorner(6, 6) = fit.uncertainty.cartesianStateCovariance;
	cov(6, 6) = 1e-26;
	for (int i = 0; i < 6; ++i)
		cov(i, 6) = cov(6, i) = 0.2 * cov(i, 0) / std::sqrt(cov(0, 0)) * 1e-13;
	const Eigen::MatrixXd original = cov;
	const auto coordinates = hazard::cometaryCoordinates(*fit.fitEpochOrbit, fit.finalNonGrav);
	ASSERT_TRUE(coordinates);
	const auto orbitalJacobian = cometary::cartesianToElementsJacobian(coordinates->head<6>());
	ASSERT_TRUE(orbitalJacobian);
	Eigen::MatrixXd jacobian = Eigen::MatrixXd::Identity(7, 7);
	jacobian.topLeftCorner<6, 6>() = *orbitalJacobian;
	Eigen::MatrixXd expected = jacobian * original * jacobian.transpose();
	expected = (0.5 * (expected + expected.transpose())).eval();
	for (int i = 2; i < 6; ++i) {
		expected.row(i) *= Constants::RAD_TO_DEG;
		expected.col(i) *= Constants::RAD_TO_DEG;
	}
	const auto report = save({ExportProduct::OrbitElements});
	const auto output = readOef(report.directory / "orbit_cometary_oef.txt");
	const auto actual = output.matrix("COV", 7);
	expected.row(6) *= 1e10;
	expected.col(6) *= 1e10;
	EXPECT_LT((actual - expected).norm() / expected.norm(), 2e-14);
	for (int i = 0; i < 6; ++i)
		EXPECT_NEAR(actual(i, 6), expected(i, 6), 2e-14 * std::abs(expected(i, 6)));
	EXPECT_NEAR(actual(6, 6), 1e-6, 1e-20);
	EXPECT_EQ(output.vector("LSP")(2), 7.0);
	EXPECT_EQ((cov - original).norm(), 0.0);
}

TEST_F(ExporterOefTest, CometaryExportMarksUnavailableCovarianceWithoutZeros)
{
	ASSERT_NO_FATAL_FAILURE(setCometaryOrbit(1.2));
	fit.uncertainty.cartesianStateCovariance.setZero();
	const auto report = save({ExportProduct::OrbitElements});
	const auto output = readOef(report.directory / "orbit_cometary_oef.txt");
	EXPECT_TRUE(output.vector("COT").allFinite());
	EXPECT_FALSE(output.fields.contains("COV"));
	EXPECT_FALSE(output.fields.contains("RMS"));
	EXPECT_FALSE(output.fields.contains("NOR"));
	std::ifstream input(report.directory / "orbit_cometary_oef.txt");
	std::ostringstream contents;
	contents << input.rdbuf();
	EXPECT_NE(contents.str().find("! COVARIANCE_UNAVAILABLE\n"), std::string::npos);
	EXPECT_TRUE(std::filesystem::exists(report.directory / "orbit_cartesian_oef.txt"));
}

TEST_F(ExporterOefTest, UnsupportedCometaryChartStillExportsCartesianState)
{
	ASSERT_NO_FATAL_FAILURE(setCometaryOrbit(1.2));
	fit.fitEpochOrbit->setR(Vec3D(0.8, 0.0, 0.0));
	fit.fitEpochOrbit->setV(Vec3D(std::sqrt(2.2 / 0.8), 0.0, 0.0));
	ExportContext context;
	context.objectName = "2020VV";
	context.fitResult = &fit;
	context.rootDirectory = root;
	EXPECT_THROW((void)Exporter::save(context, {ExportProduct::OrbitElements}), std::runtime_error);
	std::filesystem::path carFile;
	for (const auto& entry : std::filesystem::recursive_directory_iterator(root)) {
		if (entry.path().filename() == "orbit_cartesian_oef.txt")
			carFile = entry.path();
	}
	ASSERT_FALSE(carFile.empty());
	const auto output = readOef(carFile);
	EXPECT_TRUE(output.vector("CAR").allFinite());
	EXPECT_FALSE(std::filesystem::exists(carFile.parent_path() / "orbit_cometary_oef.txt"));
}

TEST_F(ExporterOefTest, CovarianceTextAppendsOptionalCometaryBlockWithItsOwnEpoch)
{
	auto saveText = [&]() {
		const auto report = save({ExportProduct::Covariance});
		std::ifstream input(report.directory / "covariance.txt", std::ios::binary);
		std::ostringstream contents;
		contents << input.rdbuf();
		return contents.str();
	};
	fit.uncertainty.hasElementCovariance = false;
	fit.uncertainty.cometaryElementCovariance.reset();
	const std::string legacy = saveText();
	ASSERT_FALSE(legacy.empty());
	EXPECT_EQ(legacy.find("cometary_element_covariance"), std::string::npos);

	CometaryElementCovariance cot;
	cot.epochJdTdb = fit.fitEpochOrbit->getEpoch() + 37.125;
	for (int r = 0; r < 6; ++r) {
		for (int c = 0; c < 6; ++c) {
			cot.formalCovariance(r, c) = (r == c ? 10.0 : 1.0) * (r + 1) * (c + 1) / 128.0;
			cot.covariance(r, c) = 2.0 * cot.formalCovariance(r, c);
		}
	}
	fit.uncertainty.cometaryElementCovariance = cot;
	const std::string output = saveText();
	EXPECT_EQ(output.substr(0, legacy.size()), legacy);
	const std::string epochTag = "cometary_element_covariance_epoch_jd_tdb ";
	const auto epochPos = output.find(epochTag);
	ASSERT_NE(epochPos, std::string::npos);
	std::istringstream epochInput(output.substr(epochPos + epochTag.size()));
	double epoch = 0.0;
	ASSERT_TRUE(static_cast<bool>(epochInput >> epoch));
	EXPECT_DOUBLE_EQ(epoch, cot.epochJdTdb);
	EXPECT_NE(output.find("cometary_element_covariance_units AU 1 rad rad rad rad"), std::string::npos);

	auto checkMatrix = [&](const char* tag, const auto& expected) {
		const std::string section = '[' + std::string(tag) + ']';
		const auto pos = output.find(section);
		ASSERT_NE(pos, std::string::npos);
		std::istringstream input(output.substr(pos + section.size()));
		std::string rowsTag, colsTag;
		int rows = 0, cols = 0;
		ASSERT_TRUE(static_cast<bool>(input >> rowsTag >> rows >> colsTag >> cols));
		ASSERT_EQ(rowsTag, "rows");
		ASSERT_EQ(colsTag, "cols");
		ASSERT_EQ(rows, 6);
		ASSERT_EQ(cols, 6);
		for (int r = 0; r < 6; ++r) {
			for (int c = 0; c < 6; ++c) {
				double actual = 0.0;
				ASSERT_TRUE(static_cast<bool>(input >> actual));
				EXPECT_DOUBLE_EQ(actual, expected(r, c));
			}
		}
	};
	checkMatrix("cometary_element_covariance_q_e_i_node_peri_nu", cot.covariance);
	checkMatrix("cometary_element_covariance_formal_q_e_i_node_peri_nu", cot.formalCovariance);
	fit.uncertainty.cometaryElementCovariance.reset();
	EXPECT_EQ(saveText(), legacy);
}
