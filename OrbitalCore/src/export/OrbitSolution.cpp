#include "orbital_core/OrbitSolution.h"
#include "orbital_core/TimeScale.h"
#include "hazard/LovCoordinateChart.h"
#include "OefLayout.h"
#include "OefCoordinates.h"
#include <Eigen/Cholesky>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace orbital_core {
	namespace {
		using Rows = std::vector<std::vector<std::string>>;

		Rows readRows(const std::filesystem::path& file)
		{
			std::ifstream stream(file);
			if (!stream)
				throw std::runtime_error("Cannot read solution file: " + file.string());
			Rows rows;
			for (std::string line; std::getline(stream, line);) {
				if (line.starts_with("! OC_NG "))
					line.erase(0, 2);
				line.erase(line.find('!') == std::string::npos ? line.size() : line.find('!'));
				std::istringstream input(line);
				input.imbue(std::locale::classic());
				std::vector<std::string> row;
				for (std::string token; input >> token;)
					row.push_back(std::move(token));
				if (!row.empty())
					rows.push_back(std::move(row));
			}
			if (stream.bad())
				throw std::runtime_error("Cannot finish reading solution file");
			return rows;
		}

		double number(std::string token)
		{
			std::replace(token.begin(), token.end(), 'D', 'E');
			std::replace(token.begin(), token.end(), 'd', 'e');
			std::istringstream input(token);
			input.imbue(std::locale::classic());
			double value;
			if (!(input >> value) || !input.eof() || !std::isfinite(value))
				throw std::invalid_argument("Invalid solution number: " + token);
			return value;
		}

		int integer(const std::string& token, int low, int high)
		{
			const double value = number(token);
			if (value < low || value > high || value != std::floor(value))
				throw std::invalid_argument("Invalid solution integer: " + token);
			return static_cast<int>(value);
		}

		const std::vector<std::string>& required(const Rows& rows, const std::string& tag, std::size_t size = 0)
		{
			const std::vector<std::string>* found = nullptr;
			for (const auto& row : rows) {
				if (row[0] != tag)
					continue;
				if (found || (size && row.size() != size))
					throw std::invalid_argument("Duplicate or malformed solution field: " + tag);
				found = &row;
			}
			if (!found)
				throw std::invalid_argument("Missing solution field: " + tag);
			return *found;
		}

		void validateCovariance(const FitParameterCovariance& covariance)
		{
			if (covariance.rows() < 6 || covariance.rows() > 10 || covariance.rows() != covariance.cols() ||
			    !covariance.allFinite() || (covariance.diagonal().array() <= 0).any())
				throw std::invalid_argument("Invalid solution covariance dimensions or diagonal");
			const Eigen::VectorXd sigma = covariance.diagonal().cwiseSqrt();
			const Eigen::MatrixXd correlation = covariance.cwiseQuotient(sigma * sigma.transpose());
			if ((correlation - correlation.transpose()).cwiseAbs().maxCoeff() > 1e-12 ||
			    Eigen::LLT<Eigen::MatrixXd>(correlation).info() != Eigen::Success)
				throw std::invalid_argument("Solution covariance is not symmetric positive definite");
		}

		bool contains(const Rows& rows, const std::string& tag)
		{
			return std::ranges::any_of(rows, [&](const auto& row) { return row[0] == tag; });
		}

		std::vector<int> readNonGrav(const Rows& body, NonGravModel& ng, bool extended)
		{
			ng.setModel(NgModel::None);
			if (!contains(body, "LSP")) {
				if (contains(body, "NGR"))
					throw std::invalid_argument("NGR without LSP");
				return {};
			}
			const auto& lsp = required(body, "LSP");
			if (lsp.size() < 4)
				throw std::invalid_argument("Incomplete LSP");
			const int model = integer(lsp[1], 0, extended ? 103 : 2);
			if ((model > 2 && model < 101) || model == 102)
				throw std::invalid_argument("Unknown extended OEF model");
			const int count = integer(lsp[2], 0, 4);
			const int fitted = integer(lsp[3], 6, 10) - 6;
			if (lsp.size() != 4 + fitted || count != (model == 0 ? 0 : model == 1 ? 2 : 4))
				throw std::invalid_argument("Inconsistent OEF LSP layout");
			if (model == 0) {
				if (fitted || contains(body, "NGR"))
					throw std::invalid_argument("Non-grav values without a model");
				return {};
			}
			const auto& ngr = required(body, "NGR", count + 1);
			if (model == 1) {
				if (number(ngr[1]) != 0.0)
					throw std::invalid_argument("OEF area/mass radiation pressure is not supported by Core");
				ng.setModel(NgModel::Yarkovsky);
				ng.setParams(0.0, number(ngr[2]) * 1e-10, 0.0);
			} else {
				ng.setModel(model == 2     ? NgModel::MarsdenH2O
				            : model == 101 ? NgModel::MarsdenCO
				                           : NgModel::InverseSquare);
				ng.setParams(number(ngr[1]), number(ngr[2]), number(ngr[3]));
				ng.setDT(number(ngr[4]));
			}
			ng.setEnabled(true);
			std::vector<int> slots;
			unsigned mask = 0;
			for (int i = 0; i < fitted; ++i) {
				const int slot = integer(lsp[4 + i], 1, count);
				if ((model == 1 && slot != 2) || std::ranges::find(slots, slot) != slots.end())
					throw std::invalid_argument("Unsupported or duplicate fitted non-grav slot");
				slots.push_back(slot);
				mask |= 1u << (slot - 1);
			}
			ng.setFitParameterMask(static_cast<NgParameterMask>(mask));
			if (ng.numFitParams() != fitted)
				throw std::invalid_argument("Unsupported non-grav solve mask");
			return slots;
		}

		void restoreNonGravMetadata(const Rows& header, NonGravModel& ng, bool extended)
		{
			if (!contains(header, "OC_NG")) {
				if (extended)
					throw std::invalid_argument("Core OEF extension requires model metadata");
				return;
			}
			const auto& row = required(header, "OC_NG", 9);
			(void)integer(row[1], 1, 3);
			const auto originalModel = static_cast<NgModel>(integer(row[2], 0, 4));
			NonGravModel restored;
			restored.setModel(originalModel);
			const bool enabled = integer(row[3], 0, 1) != 0;
			const auto mask = static_cast<NgParameterMask>(integer(row[4], 0, 15));
			if (NonGravModel::sanitizeFitParameterMaskFor(restored.model(), mask) != mask)
				throw std::invalid_argument("Invalid Core OEF fitted parameter mask");
			restored.setFitParameterMask(mask);
			restored.setParams(number(row[5]), number(row[6]), number(row[7]));
			restored.setDT(number(row[8]));
			restored.setEnabled(enabled);
			const auto saved = oef::layout(restored), native = oef::layout(ng);
			if (oef::requiresExtension(restored) != extended || restored.enabled() != enabled ||
			    saved.model != native.model || saved.solved != native.solved ||
			    saved.values.size() != native.values.size())
				throw std::invalid_argument("OEF model metadata disagrees with LSP/NGR");
			for (std::size_t i = 0; i < saved.values.size(); ++i)
				if (std::abs(saved.values[i] - native.values[i]) >
				    1e-14 * std::max(std::abs(saved.values[i]), std::abs(native.values[i])))
					throw std::invalid_argument("OEF model metadata disagrees with NGR values");
			ng = restored;
		}
	}  // namespace

	OrbitSolution loadOrbitSolution(const std::filesystem::path& path)
	{
		const auto rows = readRows(std::filesystem::is_directory(path) ? path / "orbit_cartesian_oef.txt" : path);
		const auto header = std::find(rows.begin(), rows.end(), std::vector<std::string>{"END_OF_HEADER"});
		if (header == rows.end() || header + 1 == rows.end())
			throw std::invalid_argument("Missing OEF header or object name");
		const Rows headerRows(rows.begin(), header);
		const auto& format = required(headerRows, "format", 3);
		const bool extended = format[2] == "'OC-OEF1.0'";
		if (required(headerRows, "format", 3)[1] != "=" || required(headerRows, "rectype", 3)[1] != "=" ||
		    required(headerRows, "refsys", 4)[1] != "=" || (format[2] != "'OEF2.0'" && !extended) ||
		    required(headerRows, "rectype", 3)[2] != "'ML'" || required(headerRows, "refsys", 4)[2] != "ECLM" ||
		    required(headerRows, "refsys", 4)[3] != "J2000")
			throw std::invalid_argument("Expected OEF2.0 ML, heliocentric ECLM J2000");
		OrbitSolution result;
		for (const auto& word : *(header + 1)) {
			if (!result.objectName.empty())
				result.objectName += ' ';
			result.objectName += word;
		}
		const Rows body(header + 2, rows.end());
		const auto& epoch = required(body, "MJD", 3);
		if (epoch[2] != "TDT" && epoch[2] != "TT" && epoch[2] != "TDB")
			throw std::invalid_argument("Solution epoch must be TDT/TT/TDB");
		const double mjd = time::toMjdTdb(number(epoch[1]), epoch[2] == "TDB" ? time::Scale::TDB : time::Scale::TT);
		if (!(time::jdFromMjd(mjd) > 0))
			throw std::invalid_argument("Invalid solution epoch");
		const auto slots = readNonGrav(body, result.nonGrav, extended);
		restoreNonGravMetadata(headerRows, result.nonGrav, extended);
		const int n = 6 + result.nonGrav.numFitParams();
		std::string coordinates;
		std::vector<double> cov;
		for (const auto& row : body) {
			if (row[0] == "CAR" || row[0] == "EQU" || row[0] == "KEP" || row[0] == "COT") {
				if (!coordinates.empty())
					throw std::invalid_argument("Multiple orbital states in OEF");
				coordinates = row[0];
			} else if (row[0] == "COV") {
				if (row.size() < 2)
					throw std::invalid_argument("Empty OEF covariance row");
				for (std::size_t i = 1; i < row.size(); ++i)
					cov.push_back(number(row[i]));
			} else if (row[0] != "MJD" && row[0] != "LSP" && row[0] != "NGR" && row[0] != "MAG" && row[0] != "NOR") {
				throw std::invalid_argument("Unsupported OEF record: " + row[0]);
			}
		}
		const auto& elements = required(body, coordinates, 7);
		Eigen::VectorXd values(n);
		for (int i = 0; i < 6; ++i)
			values(i) = number(elements[i + 1]);
		for (int i = 6; i < n; ++i)
			values(i) = result.nonGrav.fitParameterValue(i - 6);
		if (cov.size() != static_cast<std::size_t>(n * (n + 1) / 2))
			throw std::invalid_argument("OEF covariance is incomplete or has the wrong dimension");
		Eigen::MatrixXd fileCov(n, n);
		std::size_t k = 0;
		for (int i = 0; i < n; ++i)
			for (int j = i; j < n; ++j)
				fileCov(i, j) = fileCov(j, i) = cov[k++];
		validateCovariance(fileCov);
		const auto spec = oef::layout(result.nonGrav);
		Eigen::MatrixXd conversion = Eigen::MatrixXd::Zero(n, n);
		conversion.topLeftCorner<6, 6>().setIdentity();
		for (int i = 6; i < n; ++i) {
			const auto at = std::ranges::find(slots, spec.solved[i - 6]);
			conversion(i, 6 + std::distance(slots.begin(), at)) = 1.0 / spec.scales(i);
		}
		if (coordinates != "CAR") {
			const int first = coordinates == "EQU" ? 5 : 2;
			for (int i = first; i < 6; ++i) {
				values(i) *= Constants::DEG_TO_RAD;
				conversion.row(i) *= Constants::DEG_TO_RAD;
			}
		}
		result.covariance = conversion * fileCov * conversion.transpose();
		if (coordinates != "CAR") {
			hazard::LovCoordinateChart chart;
			chart.valid = true;
			chart.nonGrav = result.nonGrav;
			chart.epochTdb = time::jdFromMjd(mjd);
			chart.coordinateSystem = coordinates == "EQU"   ? hazard::LovCoordinateSystem::Equinoctial
			                         : coordinates == "KEP" ? hazard::LovCoordinateSystem::Keplerian
			                                                : hazard::LovCoordinateSystem::Cometary;
			chart.nominalCoordinates = values;
			const auto physical = hazard::physicalParametersFromChart(chart, values);
			const std::optional<Eigen::MatrixXd> jacobian =
			    coordinates == "KEP" ? std::optional<Eigen::MatrixXd>{oef::keplerianToPhysicalJacobian(values)}
			                         : hazard::physicalJacobianFromChart(chart, values);
			if (!physical || !jacobian)
				throw std::invalid_argument("Invalid OEF elements or singular coordinate conversion");
			values = *physical;
			result.covariance = (*jacobian * result.covariance * jacobian->transpose()).eval();
			result.covariance = (0.5 * (result.covariance + result.covariance.transpose())).eval();
		}
		if (contains(headerRows, "OC_NG") && required(headerRows, "OC_NG", 9)[1] == "1" &&
		    result.nonGrav.model() == NgModel::MarsdenCO) {
			// Version 1 used the r0=10 normalizer with the r0=5 activity law.
			constexpr double scale = 0.01003 / 0.0408373333128795;
			result.nonGrav.setParams(result.nonGrav.A1() * scale, result.nonGrav.A2() * scale,
			                         result.nonGrav.A3() * scale);
			for (int i = 0; i < result.nonGrav.fitParameterCount(); ++i) {
				if (result.nonGrav.fitParameterAt(i) != NgParameter::DT) {
					result.covariance.row(6 + i) *= scale;
					result.covariance.col(6 + i) *= scale;
				}
			}
		}
		auto orbit =
		    make_orbit_from_states_vectors(StateVectors(Vec3D(values(0), values(1), values(2)),
		                                                Vec3D(values(3), values(4), values(5)) / Constants::GAUSS_K),
		                                   time::jdFromMjd(mjd));
		if (!orbit)
			throw std::invalid_argument("Invalid Cartesian state in solution");
		orbit->setEpochMjd(mjd);
		result.nominal = *orbit;
		validateCovariance(result.covariance);
		if (contains(body, "MAG")) {
			const auto& mag = required(body, "MAG", 3);
			const double h = number(mag[1]);
			result.slopeG = number(mag[2]);
			if (h != 0 && h != 99.999)
				result.absoluteMagnitude = h;
		}
		return result;
	}
}  // namespace orbital_core
