#include "orbital_core/Exporter.h"
#include "OefLayout.h"
#include "obs/RwoMetadata.h"
#include "OefCoordinates.h"
#include "EnvironmentVariable.h"
#include <filesystem>
#include <span>
#include <string>
#include <utility>
#include <cstdint>

#include "orbital_core/Constants.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/StateVectors.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/Moid.h"

#include "hazard/LovCoordinateChart.h"
#include "hazard/CloseApproachKinematics.h"
#include "orbital_core/Photometry.h"
#include "orbital_core/uncertainty/CovarianceDiagnostics.h"

#include <Eigen/Cholesky>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <unordered_set>
#include <vector>

namespace orbital_core {
	namespace {

		[[nodiscard]] std::string localDateStamp()
		{
			const auto now = std::chrono::system_clock::now();
			const std::time_t tt = std::chrono::system_clock::to_time_t(now);
			std::tm local{};
#ifdef _WIN32
			localtime_s(&local, &tt);
#else
			localtime_r(&tt, &local);
#endif

			std::ostringstream os;
			os << std::put_time(&local, "%Y-%m-%d");
			return os.str();
		}

		[[nodiscard]] std::string sanitizePathPart(std::string text)
		{
			if (text.empty())
				return "object";

			for (char& ch : text) {
				const unsigned char uch = static_cast<unsigned char>(ch);
				if (uch < 32 || ch == '<' || ch == '>' || ch == ':' || ch == '"' || ch == '/' || ch == '\\' ||
				    ch == '|' || ch == '?' || ch == '*') {
					ch = '_';
				}
			}

			while (!text.empty() && (text.back() == ' ' || text.back() == '.'))
				text.pop_back();
			while (!text.empty() && text.front() == ' ')
				text.erase(text.begin());

			return text.empty() ? std::string("object") : text;
		}

		[[nodiscard]] std::filesystem::path appendWorkDirectory(std::filesystem::path root)
		{
			if (root.filename() == "work")
				return root;
			return root / "work";
		}

		[[nodiscard]] std::filesystem::path localAppDataWorkDirectory()
		{
#ifdef _WIN32
			if (const auto localAppData = environmentVariable(L"LOCALAPPDATA"); localAppData && !localAppData->empty()) {
				return std::filesystem::path(*localAppData) / "Trajecta" / "work";
			}
			const auto home = environmentVariable(L"HOME");
#elif defined(__linux__)
			if (const auto data = environmentVariable("XDG_DATA_HOME"); data && !data->empty()) {
				return std::filesystem::path(*data) / "Trajecta" / "work";
			}
			const auto home = environmentVariable("HOME");
#else
			const auto home = environmentVariable("HOME");
#endif
			if (home && !home->empty()) {
				return std::filesystem::path(*home) / ".local" / "share" / "Trajecta" / "work";
			}
			return {};
		}

		[[nodiscard]] bool ensureDirectory(const std::filesystem::path& dir, std::error_code& ec)
		{
			ec.clear();
			std::filesystem::create_directories(dir, ec);
			if (ec)
				return false;

			std::error_code existsEc;
			return std::filesystem::is_directory(dir, existsEc) && !existsEc;
		}

		[[nodiscard]] std::filesystem::path resolveBaseDirectory(const ExportContext& context)
		{
#if defined(__linux__)
			if (context.rootDirectory.empty()) {
				const auto data = localAppDataWorkDirectory();
				if (!data.empty())
					return data;
			}
#endif
			std::filesystem::path root = context.rootDirectory.empty() ? getExecutableDir() : context.rootDirectory;
			if (root.empty())
				throw std::runtime_error("Exporter root directory is empty");

			return appendWorkDirectory(std::move(root));
		}

		[[nodiscard]] std::filesystem::path resolveDirectory(const ExportContext& context)
		{
			const std::filesystem::path leaf =
			    std::filesystem::path(localDateStamp()) / sanitizePathPart(context.objectName);

			const std::filesystem::path primary = resolveBaseDirectory(context);
			std::error_code primaryEc;
			const std::filesystem::path primaryDir = primary / leaf;
			if (ensureDirectory(primaryDir, primaryEc))
				return primaryDir;

			const std::filesystem::path fallback = localAppDataWorkDirectory();
			std::error_code fallbackEc;
			if (!fallback.empty() && fallback != primary) {
				const std::filesystem::path fallbackDir = fallback / leaf;
				if (ensureDirectory(fallbackDir, fallbackEc))
					return fallbackDir;
			}

			std::ostringstream os;
			os << "failed to create export directory: " << primaryDir.string();
			if (primaryEc)
				os << " (" << primaryEc.message() << ')';
			if (!fallback.empty()) {
				os << "; fallback failed: " << (fallback / leaf).string();
				if (fallbackEc)
					os << " (" << fallbackEc.message() << ')';
			}
			throw std::runtime_error(os.str());
		}

		void requireFit(const ExportContext& context)
		{
			if (!context.fitResult)
				throw std::invalid_argument("Exporter requires FitResult");
		}

		[[nodiscard]] std::ofstream openTextFile(const std::filesystem::path& file)
		{
			std::ofstream os(file, std::ios::out | std::ios::trunc);
			if (!os)
				throw std::runtime_error("failed to open export file: " + file.string());
			os << std::setprecision(17);
			return os;
		}

		[[nodiscard]] std::string objectDesignator(const ExportContext& context)
		{
			return sanitizePathPart(context.objectName);
		}

		[[nodiscard]] double exportedAbsoluteMagnitude(const FitResult& fit, double slopeG)
		{
			const double H =
			    photometry::fit(photometry::observations(fit).data(), slopeG).absoluteMagnitude.value_or(0.0);
			return std::isfinite(H) && H > 0.0 ? H : 0.0;
		}

		void writeOefHeader(std::ostream& os, const NonGravModel& ng)
		{
			os << "format  = '" << (oef::requiresExtension(ng) ? "OC-OEF1.0" : "OEF2.0") << "'       ! file format\n";
			os << "rectype = 'ML'           ! record type (1L/ML)\n";
			os << "refsys  = ECLM J2000     ! default reference system\n";
			if (oef::requiresExtension(ng))
				os << "! Core extension: this force law is not supported by standard OEF\n";
			auto configured = ng;
			configured.setEnabled(true);
			os << "! OC_NG 3 " << static_cast<int>(ng.model()) << ' ' << ng.enabled() << ' '
			   << static_cast<int>(configured.fitParameterMask()) << std::setprecision(17) << ' ' << ng.A1() << ' '
			   << ng.A2() << ' ' << ng.A3() << ' ' << ng.DT() << '\n';
			os << "END_OF_HEADER\n";
		}

		void writeMagnitudeLine(std::ostream& os, const ExportContext& context)
		{
			const double h = context.cometPhotometry
			                     ? 0.0
			                     : exportedAbsoluteMagnitude(*context.fitResult, context.photometry.asteroidG);
			os << " MAG " << std::fixed << std::setprecision(3) << std::setw(7) << (h > 0.0 ? h : 99.999)
			   << std::setw(7) << context.photometry.asteroidG << '\n'
			   << std::defaultfloat << std::setprecision(17);
		}

		void writeLspLine(std::ostream& os, const FitResult& fit)
		{
			const auto spec = oef::layout(fit.finalNonGrav);
			os << "! Non-gravitational parameters: model used, number of model parameters, dimension\n";
			os << " LSP   " << spec.model << "  " << spec.values.size() << "    " << spec.scales.size();
			for (int slot : spec.solved)
				os << ' ' << slot;
			os << '\n';
			if (!spec.values.empty()) {
				os << " NGR" << std::scientific << std::setprecision(17);
				for (double value : spec.values)
					os << ' ' << value;
				os << '\n' << std::defaultfloat << std::setprecision(17);
			}
		}

		void writeOefScientific(std::ostream& os, double value, int width = 25, int precision = 16);

		struct NodalDistances {
			double ascendingAu{std::numeric_limits<double>::quiet_NaN()};
			double descendingAu{std::numeric_limits<double>::quiet_NaN()};
		};

		[[nodiscard]] bool finiteVector(const Vec3D& v) noexcept
		{
			return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
		}

		[[nodiscard]] std::optional<NodalDistances> nodalDistances(const Vec3D& objectR, const Vec3D& objectVCanonical,
		                                                           const Vec3D& earthR,
		                                                           const Vec3D& earthVCanonical) noexcept
		{
			if (!finiteVector(objectR) || !finiteVector(objectVCanonical) || !finiteVector(earthR) ||
			    !finiteVector(earthVCanonical))
				return std::nullopt;

			const Vec3D cObject = objectR.cross(objectVCanonical);
			const Vec3D cEarth = earthR.cross(earthVCanonical);
			const Vec3D node = cEarth.cross(cObject);
			const double nodeNorm = node.magnitude();
			const double objectRNorm = objectR.magnitude();
			const double earthRNorm = earthR.magnitude();
			if (nodeNorm <= 0.0 || objectRNorm <= 0.0 || earthRNorm <= 0.0)
				return std::nullopt;

			constexpr double kObjectMu = 1.0;
			constexpr double kEarthOrbitMu = 1.0 + Constants::GM_EARTH_PLUS_MOON / Constants::GM_SUN;

			Vec3D objectLenz = objectVCanonical.cross(cObject) / kObjectMu - objectR / objectRNorm;
			Vec3D earthLenz = earthVCanonical.cross(cEarth) / kEarthOrbitMu - earthR / earthRNorm;

			const double objectE = objectLenz.magnitude();
			const double earthE = earthLenz.magnitude();
			if (objectE <= 0.0 || earthE <= 0.0)
				return std::nullopt;

			const double cosObject = std::clamp(node.dot(objectLenz) / (nodeNorm * objectE), -1.0, 1.0);
			const double cosEarth = std::clamp(node.dot(earthLenz) / (nodeNorm * earthE), -1.0, 1.0);

			const double objectC2 = cObject.dot(cObject);
			const double earthC2 = cEarth.dot(cEarth);
			const double objectQ = (objectC2 / kObjectMu) / (1.0 + objectE);
			const double earthQ = (earthC2 / kEarthOrbitMu) / (1.0 + earthE);

			const double objectRp = objectQ * (1.0 + objectE) / (1.0 + objectE * cosObject);
			const double objectRm = objectQ * (1.0 + objectE) / (1.0 - objectE * cosObject);
			const double earthRp = earthQ * (1.0 + earthE) / (1.0 + earthE * cosEarth);
			const double earthRm = earthQ * (1.0 + earthE) / (1.0 - earthE * cosEarth);

			if (!std::isfinite(objectRp) || !std::isfinite(objectRm) || !std::isfinite(earthRp) ||
			    !std::isfinite(earthRm))
				return std::nullopt;

			NodalDistances out;
			out.ascendingAu = objectRp - earthRp;
			out.descendingAu = objectRm - earthRm;
			if (objectRp < 0.0 || earthRp < 0.0)
				out.ascendingAu = out.descendingAu;
			else if (objectRm < 0.0 || earthRm < 0.0)
				out.descendingAu = out.ascendingAu;
			return out;
		}

		struct EarthOrbitUsefulData {
			double moidAu{std::numeric_limits<double>::quiet_NaN()};
			std::optional<NodalDistances> nodes{};
		};

		[[nodiscard]] std::optional<EarthOrbitUsefulData> earthOrbitUsefulData(const Orbit& orbit) noexcept
		{
			IEphemerisSource* eph = getDefaultEphemerisSource();
			if (!eph || !eph->isAvailable())
				return std::nullopt;

			EphemerisRecord earth;
			try {
				earth = eph->queryJD(EphemerisPlanet::Earth, orbit.getEpoch());
			} catch (...) {
				return std::nullopt;
			}

			const Vec3D earthVCanonical = earth.V * (1.0 / Constants::GAUSS_K);
			EarthOrbitUsefulData out;
			out.nodes = nodalDistances(orbit.getR(), orbit.getV(), earth.R, earthVCanonical);

			auto earthOrbit = make_orbit_from_states_vectors(StateVectors(earth.R, earthVCanonical), orbit.getEpoch());
			if (earthOrbit) {
				const auto moid = hazard::computeMoid(orbit, *earthOrbit);
				if (moid.isValid())
					out.moidAu = moid.moidAu;
			}
			return out;
		}

		void writeKeplerianUsefulData(std::ostream& os, const Orbit& orbit, const FitResult& fit, bool comet,
		                              double slopeG)
		{
			os << "! PERIHELION ";
			writeOefScientific(os, orbit.getPerigee(), 23);
			os << '\n';
			os << "! APHELION ";
			writeOefScientific(os, orbit.getApogee(), 23);
			os << '\n';

			const auto earthData = earthOrbitUsefulData(orbit);
			if (earthData && earthData->nodes) {
				os << "! ANODE ";
				writeOefScientific(os, earthData->nodes->ascendingAu, 23);
				os << '\n';
				os << "! DNODE ";
				writeOefScientific(os, earthData->nodes->descendingAu, 23);
				os << '\n';
			}
			if (earthData && std::isfinite(earthData->moidAu)) {
				os << "! MOID ";
				writeOefScientific(os, earthData->moidAu, 23);
				os << '\n';
			}

			os << "! PERIOD ";
			writeOefScientific(os, orbit.getOrbitalPeriod(), 23);
			os << '\n';

			if (!comet && earthData && std::isfinite(earthData->moidAu)) {
				const double H = exportedAbsoluteMagnitude(fit, slopeG);
				const bool pha = H > 0.0 && H < 22.0 && earthData->moidAu <= 0.05;
				os << "! PHA " << (pha ? 'T' : 'F') << '\n';
			}

			const double vInf = hazard::circularEarthVInfinityKmS(orbit);
			if (std::isfinite(vInf)) {
				os << "! VINFTY ";
				writeOefScientific(os, vInf, 23);
				os << '\n';
			}

			if (std::isfinite(fit.uncertainty.orbitQualityU)) {
				std::ostringstream uLine;
				uLine << std::fixed << std::setprecision(1) << fit.uncertainty.orbitQualityU;
				os << "! U_PAR " << uLine.str() << '\n';
			}
		}

		template <class Matrix>
		[[nodiscard]] bool usableCovariance(const Matrix& matrix)
		{
			if (matrix.rows() <= 0 || matrix.cols() <= 0 || matrix.rows() != matrix.cols())
				return false;

			bool any = false;
			for (Eigen::Index r = 0; r < matrix.rows(); ++r) {
				for (Eigen::Index c = 0; c < matrix.cols(); ++c) {
					const double v = matrix(r, c);
					if (!std::isfinite(v))
						return false;
					any = any || v != 0.0;
				}
			}
			return any;
		}

		template <class Matrix>
		void writeUpperTriangular(std::ostream& os, const char* tag, const Matrix& matrix)
		{
			std::vector<double> values;
			values.reserve(static_cast<std::size_t>(matrix.rows() * (matrix.rows() + 1) / 2));
			for (Eigen::Index r = 0; r < matrix.rows(); ++r) {
				for (Eigen::Index c = r; c < matrix.cols(); ++c)
					values.push_back(matrix(r, c));
			}

			const auto oldFlags = os.flags();
			const auto oldPrecision = os.precision();
			os << std::uppercase << std::scientific << std::setprecision(15);
			for (std::size_t i = 0; i < values.size(); i += 3) {
				os << ' ' << tag;
				const std::size_t end = std::min(i + 3, values.size());
				for (std::size_t j = i; j < end; ++j)
					os << ' ' << std::setw(23) << values[j];
				os << '\n';
			}
			os.flags(oldFlags);
			os.precision(oldPrecision);
		}

		template <class Matrix>
		void writeRmsLine(std::ostream& os, const Matrix& cov)
		{
			if (!usableCovariance(cov))
				return;

			const auto oldFlags = os.flags();
			const auto oldPrecision = os.precision();
			os << "! RMS ";
			os << std::uppercase << std::scientific << std::setprecision(5);
			for (Eigen::Index i = 0; i < cov.rows(); ++i) {
				const double variance = cov(i, i);
				os << ' ' << std::setw(13) << (variance > 0.0 ? std::sqrt(variance) : 0.0);
			}
			os << '\n';
			os.flags(oldFlags);
			os.precision(oldPrecision);
		}

		template <class Matrix>
		void writeNormalFromCovariance(std::ostream& os, const Matrix& cov)
		{
			if (!usableCovariance(cov))
				return;

			const Eigen::MatrixXd dense = cov.template selfadjointView<Eigen::Upper>();
			if ((dense.diagonal().array() <= 0.0).any())
				return;
			const Eigen::VectorXd sigma = dense.diagonal().cwiseSqrt();
			Eigen::MatrixXd scaled(dense.rows(), dense.cols());
			for (Eigen::Index r = 0; r < dense.rows(); ++r) {
				for (Eigen::Index c = r; c < dense.cols(); ++c)
					scaled(r, c) = scaled(c, r) = dense(r, c) / sigma(r) / sigma(c);
			}
			Eigen::LLT<Eigen::MatrixXd> cholesky(scaled);
			if (cholesky.info() != Eigen::Success)
				return;
			Eigen::MatrixXd inverse = cholesky.solve(Eigen::MatrixXd::Identity(dense.rows(), dense.cols()));
			inverse = (0.5 * (inverse + inverse.transpose())).eval();
			for (Eigen::Index r = 0; r < dense.rows(); ++r) {
				for (Eigen::Index c = r; c < dense.cols(); ++c)
					inverse(r, c) = inverse(c, r) = inverse(r, c) / sigma(r) / sigma(c);
			}
			if (inverse.allFinite())
				writeUpperTriangular(os, "NOR", inverse);
		}

		template <class Matrix>
		void writeCovarianceAndNormal(std::ostream& os, const Matrix& cov)
		{
			if (!usableCovariance(cov))
				return;

			writeRmsLine(os, cov);
			writeUpperTriangular(os, "COV", cov);
			writeNormalFromCovariance(os, cov);
		}

		void writeOefScientific(std::ostream& os, double value, int width, int precision)
		{
			const auto oldFlags = os.flags();
			const auto oldPrecision = os.precision();
			os << std::uppercase << std::scientific << std::setprecision(precision) << std::setw(width) << value;
			os.flags(oldFlags);
			os.precision(oldPrecision);
		}

		void writeOefFixed(std::ostream& os, double value, int width = 19, int precision = 13)
		{
			const auto oldFlags = os.flags();
			const auto oldPrecision = os.precision();
			os << std::fixed << std::setprecision(precision) << std::setw(width) << value;
			os.flags(oldFlags);
			os.precision(oldPrecision);
		}

		void writeOefMjd(std::ostream& os, double mjdTdb)
		{
			const auto oldFlags = os.flags();
			const auto oldPrecision = os.precision();
			os << " MJD" << std::fixed << std::setprecision(12) << std::setw(23) << mjdTdb << " TDB\n";
			os.flags(oldFlags);
			os.precision(oldPrecision);
		}

		void writeCartesianOefValues(std::ostream& os, const Vec3D& r, const Vec3D& v)
		{
			os << " CAR";
			writeOefScientific(os, r.x);
			writeOefScientific(os, r.y);
			writeOefScientific(os, r.z);
			writeOefScientific(os, v.x);
			writeOefScientific(os, v.y);
			writeOefScientific(os, v.z);
			os << '\n';
		}

		void writeEquinoctialOefValues(std::ostream& os, const Eigen::VectorXd& values)
		{
			os << " EQU";
			writeOefScientific(os, values(0));
			writeOefScientific(os, values(1));
			writeOefScientific(os, values(2));
			writeOefScientific(os, values(3));
			writeOefScientific(os, values(4));
			writeOefScientific(os, values(5));
			os << '\n';
		}

		void writeKeplerianOefValues(std::ostream& os, const Eigen::VectorXd& values)
		{
			os << " KEP";
			writeOefScientific(os, values(0));
			writeOefScientific(os, values(1));
			writeOefFixed(os, values(2), 18, 13);
			writeOefFixed(os, values(3), 18, 13);
			writeOefFixed(os, values(4), 18, 13);
			writeOefScientific(os, values(5), 25, 16);
			os << '\n';
		}

		void writeCometaryOefValues(std::ostream& os, const Eigen::VectorXd& values)
		{
			os << " COT";
			writeOefScientific(os, values(0));
			writeOefScientific(os, values(1));
			for (int i = 2; i < 6; ++i)
				writeOefFixed(os, values(i), 18, 13);
			os << '\n';
		}

		[[nodiscard]] FitParameterCovariance physicalCovarianceForExport(const FitResult& fit)
		{
			if (usableCovariance(fit.uncertainty.fitParameterCovariance))
				return fit.uncertainty.fitParameterCovariance;

			if (fit.finalNonGrav.numFitParams() != 0)
				throw std::invalid_argument("Cannot export non-grav without full parameter covariance");

			FitParameterCovariance cov(6, 6);
			cov = fit.uncertainty.cartesianStateCovariance;
			return cov;
		}

		void scaleAngleCoordinatesToDegrees(Eigen::VectorXd& value, Eigen::MatrixXd* cov, int firstAngle, int lastAngle)
		{
			for (int i = firstAngle; i <= lastAngle && i < value.size(); ++i)
				value(i) = Convertors::radiansToDegrees(value(i));

			if (!cov)
				return;

			for (int i = firstAngle; i <= lastAngle && i < cov->rows(); ++i)
				cov->row(i) *= Constants::RAD_TO_DEG;
			for (int i = firstAngle; i <= lastAngle && i < cov->cols(); ++i)
				cov->col(i) *= Constants::RAD_TO_DEG;
		}

		[[nodiscard]] Eigen::VectorXd equinoctialValuesForOef(const Orbit& orbit, const NonGravModel& ng)
		{
			auto coordinates = hazard::equinoctialCoordinates(orbit, ng);
			if (!coordinates)
				throw std::runtime_error("cannot export EQU for an unsupported Cartesian state");
			Eigen::VectorXd out = *coordinates;
			out(5) = Convertors::radiansToDegrees(out(5));
			return out;
		}

		[[nodiscard]] std::optional<Eigen::MatrixXd> equinoctialCovarianceRadians(
		    const Orbit& orbit, const FitParameterCovariance& physicalCov, const NonGravModel& ng)
		{
			auto chart = hazard::makeEquinoctialLovChart(orbit, physicalCov, ng);
			if (!chart.valid)
				return std::nullopt;

			return chart.coordinateCovariance;
		}

		template <class Matrix>
		void writeEigenAndWeakAxis(std::ostream& os, const Matrix& cov)
		{
			if (!usableCovariance(cov))
				return;

			const auto spectrum = uncertainty::covarianceSpectrum(cov);
			if (!spectrum)
				return;

			const auto oldFlags = os.flags();
			const auto oldPrecision = os.precision();
			os << std::uppercase << std::scientific << std::setprecision(5);
			os << "! EIG ";
			for (Eigen::Index i = 0; i < spectrum->signedSigmas.size(); ++i) {
				os << ' ' << std::setw(13) << spectrum->signedSigmas(i);
			}
			os << '\n';

			os << std::fixed << std::setprecision(5);
			os << "! WEA ";
			for (Eigen::Index i = 0; i < spectrum->weakDirection.size(); ++i)
				os << ' ' << std::setw(9) << spectrum->weakDirection(i);
			os << '\n';
			os.flags(oldFlags);
			os.precision(oldPrecision);
		}

		[[nodiscard]] const Orbit_p& exportStateOrbit(const FitResult& fit)
		{
			return fit.fitEpochOrbit ? fit.fitEpochOrbit : fit.orbit;
		}

		void writeCartesianOef(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			const auto file = dir / "orbit_cartesian_oef.txt";
			if (std::find(report.files.begin(), report.files.end(), file) != report.files.end())
				return;
			auto os = openTextFile(file);
			const FitResult& fit = *context.fitResult;
			const Orbit_p& orbit = exportStateOrbit(fit);
			if (!orbit)
				throw std::runtime_error("cannot export cartesian OEF without orbit");

			const Vec3D r = orbit->getR();
			const Vec3D v = orbit->getV() * Constants::GAUSS_K;

			writeOefHeader(os, fit.finalNonGrav);
			os << objectDesignator(context) << '\n';
			os << "! Cartesian state: x, y, z [AU], vx, vy, vz [AU/day]\n";
			writeCartesianOefValues(os, r, v);
			writeOefMjd(os, orbit->epochMjd());
			writeMagnitudeLine(os, context);
			writeLspLine(os, fit);
			writeCovarianceAndNormal(os, oef::scaledCovariance(physicalCovarianceForExport(fit), fit.finalNonGrav));

			report.files.push_back(file);
		}

		void writeEquinoctialOef(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			const auto file = dir / "orbit_equinoctial_oef.txt";
			auto os = openTextFile(file);
			const FitResult& fit = *context.fitResult;
			const Orbit_p& orbit = exportStateOrbit(fit);
			if (!orbit)
				throw std::runtime_error("cannot export equinoctial OEF without orbit");

			const Eigen::VectorXd values = equinoctialValuesForOef(*orbit, fit.finalNonGrav);
			const FitParameterCovariance physicalCov = physicalCovarianceForExport(fit);
			const auto covRadians = equinoctialCovarianceRadians(*orbit, physicalCov, fit.finalNonGrav);

			writeOefHeader(os, fit.finalNonGrav);
			os << objectDesignator(context) << '\n';
			os << "! Equinoctial elements: a, e*sin(LP), e*cos(LP), tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.\n";
			writeEquinoctialOefValues(os, values);
			writeOefMjd(os, orbit->epochMjd());
			writeMagnitudeLine(os, context);
			writeLspLine(os, fit);
			if (covRadians) {
				Eigen::MatrixXd cov = oef::scaledCovariance(*covRadians, fit.finalNonGrav);
				Eigen::VectorXd dummy = Eigen::VectorXd::Zero(cov.rows());
				scaleAngleCoordinatesToDegrees(dummy, &cov, 5, 5);
				writeRmsLine(os, cov);
				writeEigenAndWeakAxis(os, oef::scaledCovariance(*covRadians, fit.finalNonGrav));
				writeUpperTriangular(os, "COV", cov);
				writeNormalFromCovariance(os, cov);
			}

			report.files.push_back(file);
		}

		void writeKeplerianOef(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			const auto file = dir / "orbit_keplerian_oef.txt";
			auto os = openTextFile(file);
			const FitResult& fit = *context.fitResult;
			const Orbit_p& orbit = exportStateOrbit(fit);
			if (!orbit)
				throw std::runtime_error("cannot export keplerian OEF without orbit");

			const auto radians = hazard::keplerianCoordinates(*orbit, fit.finalNonGrav);
			if (!radians)
				throw std::runtime_error("cannot export keplerian OEF for this orbit");
			Eigen::VectorXd values = *radians;
			scaleAngleCoordinatesToDegrees(values, nullptr, 2, 5);

			const FitParameterCovariance physicalCov = physicalCovarianceForExport(fit);
			std::optional<Eigen::MatrixXd> cov;
			if (usableCovariance(physicalCov)) {
				const auto jacobian = oef::physicalToKeplerianJacobian(*radians);
				cov = jacobian * physicalCov * jacobian.transpose();
				*cov = (0.5 * (*cov + cov->transpose())).eval();
				Eigen::VectorXd dummy = Eigen::VectorXd::Zero(cov->rows());
				scaleAngleCoordinatesToDegrees(dummy, &*cov, 2, 5);
			}

			writeOefHeader(os, fit.finalNonGrav);
			os << objectDesignator(context) << '\n';
			os << "! Keplerian elements: a, e, i, long. node, arg. peric., mean anomaly\n";
			writeKeplerianOefValues(os, values);
			writeOefMjd(os, orbit->epochMjd());
			writeMagnitudeLine(os, context);
			writeLspLine(os, fit);
			writeKeplerianUsefulData(os, *orbit, fit, context.cometPhotometry, context.photometry.asteroidG);
			if (cov)
				writeCovarianceAndNormal(os, oef::scaledCovariance(*cov, fit.finalNonGrav));

			report.files.push_back(file);
		}

		void writeCometaryOef(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			const FitResult& fit = *context.fitResult;
			const Orbit_p& orbit = exportStateOrbit(fit);
			if (!orbit)
				throw std::runtime_error("cannot export cometary OEF without orbit");

			const auto radians = hazard::cometaryCoordinates(*orbit, fit.finalNonGrav);
			if (!radians || radians->size() < 6 || !radians->allFinite())
				throw std::runtime_error("cannot export COT for an unsupported Cartesian state");
			Eigen::VectorXd values = *radians;
			scaleAngleCoordinatesToDegrees(values, nullptr, 2, 5);

			const FitParameterCovariance physicalCov = physicalCovarianceForExport(fit);
			std::optional<Eigen::MatrixXd> cov;
			if (usableCovariance(physicalCov)) {
				auto chart = hazard::makeCometaryLovChart(*orbit, physicalCov, fit.finalNonGrav);
				if (chart.valid && chart.coordinateCovariance.rows() == values.size() &&
				    usableCovariance(chart.coordinateCovariance)) {
					cov = std::move(chart.coordinateCovariance);
					Eigen::VectorXd dummy = Eigen::VectorXd::Zero(cov->rows());
					scaleAngleCoordinatesToDegrees(dummy, &*cov, 2, 5);
				}
			}

			const auto file = dir / "orbit_cometary_oef.txt";
			auto os = openTextFile(file);
			writeOefHeader(os, fit.finalNonGrav);
			os << objectDesignator(context) << '\n';
			os << "! Cometary elements: q [AU], e, i, long. node, arg. peric., true anomaly [deg]\n";
			writeCometaryOefValues(os, values);
			writeOefMjd(os, orbit->epochMjd());
			writeMagnitudeLine(os, context);
			writeLspLine(os, fit);
			if (cov && cov->allFinite()) {
				*cov = oef::scaledCovariance(*cov, fit.finalNonGrav);
				if ((cov->diagonal().array() >= 0.0).all())
					writeRmsLine(os, *cov);
				writeUpperTriangular(os, "COV", *cov);
				writeNormalFromCovariance(os, *cov);
			} else {
				os << "! COVARIANCE_UNAVAILABLE\n";
			}
			report.files.push_back(file);
		}

		[[nodiscard]] std::string fixedWidthRwoObjectName(std::string name)
		{
			name = sanitizePathPart(std::move(name));
			std::erase(name, ' ');
			if (name.size() > 9)
				name.resize(9);
			return name;
		}

		[[nodiscard]] double positiveOr(double value, double fallback) noexcept
		{
			return std::isfinite(value) && value > 0.0 ? value : fallback;
		}

		[[nodiscard]] int decimalDigitsForAccuracy(double accuracy) noexcept
		{
			if (!(std::isfinite(accuracy) && accuracy > 0.0))
				return 6;
			const double digits = -std::log10(accuracy);
			const int rounded = static_cast<int>(std::lround(digits));
			return std::clamp(rounded, 0, 10);
		}

		[[nodiscard]] std::string dateFieldsForRwo(double jd, double accuracyDays)
		{
			int year = 0;
			int month = 0;
			double dayFrac = 0.0;
			Convertors::jdToYMD(jd, year, month, dayFrac);

			const int decimals = decimalDigitsForAccuracy(accuracyDays);
			const int dayWidth = decimals > 0 ? 3 + decimals : 2;

			std::ostringstream day;
			day << std::setfill('0') << std::fixed << std::setprecision(decimals) << std::setw(dayWidth) << dayFrac;

			std::ostringstream os;
			os << std::setw(4) << year << ' ' << std::setfill('0') << std::setw(2) << month << ' ' << std::setfill(' ')
			   << std::left << std::setw(13) << day.str() << std::right << std::setfill(' ');
			return os.str();
		}

		[[nodiscard]] std::string raFields(double radians)
		{
			int h = 0;
			int m = 0;
			double s = 0.0;
			Convertors::radiansToHMS(Helpers::wrapTwoPi(radians), h, m, s);

			std::ostringstream os;
			os << std::setfill('0') << std::setw(2) << h << ' ' << std::setw(2) << m << ' ' << std::fixed
			   << std::setprecision(3) << std::setw(6) << s << std::setfill(' ');
			return os.str();
		}

		[[nodiscard]] std::string decFields(double radians)
		{
			int d = 0;
			int m = 0;
			double s = 0.0;
			Convertors::radiansToDMS(radians, d, m, s);

			std::ostringstream os;
			os << (radians < 0.0 ? '-' : '+') << std::setfill('0') << std::setw(2) << std::abs(d) << ' ' << std::setw(2)
			   << m << ' ' << std::fixed << std::setprecision(2) << std::setw(5) << s << std::setfill(' ');
			return os.str();
		}

		[[nodiscard]] const Residual* residualForObservation(const FitResult& fit, const ObsNode& obs,
		                                                     std::size_t index)
		{
			const double obsJdTdb = obs.epochJdTdb();
			const auto& residuals = fit.residuals.data();
			if (!fit.keptObsIndices.empty()) {
				for (std::size_t row = 0; row < fit.keptObsIndices.size() && row < residuals.size(); ++row) {
					if (fit.keptObsIndices[row] == index)
						return &residuals[row];
				}
			} else if (index < residuals.size()) {
				return &residuals[index];
			}

			const auto& allResiduals = fit.allResiduals.data();
			if (index < allResiduals.size())
				return &allResiduals[index];

			if (index < residuals.size()) {
				const Residual& candidate = residuals[index];
				if (std::abs(candidate.jd - obsJdTdb) < 1.0e-7 && candidate.obsCode == obs.obsCode) {
					return &candidate;
				}
			}

			for (const Residual& r : residuals) {
				if (std::abs(r.jd - obsJdTdb) < 1.0e-7 && r.obsCode == obs.obsCode) {
					return &r;
				}
			}
			return nullptr;
		}

		[[nodiscard]] const char* errmodName(AstrometricWeightProfile profile) noexcept
		{
			switch (profile) {
				case AstrometricWeightProfile::Vfcc17:
					return "vfcc17";
				case AstrometricWeightProfile::Uniform:
					return "uniform";
			}
			return "unknown";
		}

		struct WeightedResidualRow {
			const ObsNode* obs{nullptr};
			const Residual* residual{nullptr};
			double draArcsec{0.0};
			double ddecArcsec{0.0};
			double chi{0.0};
			double sigmaRaArcsec{1.0};
			double sigmaDecArcsec{1.0};
			bool forceWeightRa{false};
			bool forceWeightDec{false};
			std::string magString;
			double magSigma{-1.0};
			double magResidual{0.0};
			int magSelected{0};
			bool magResidualDefined{false};
			bool selected{true};
		};

		[[nodiscard]] std::string formatFixedField(double value, int width, int precision)
		{
			if (!std::isfinite(value))
				throw std::invalid_argument("Non-finite RWO field");
			std::ostringstream os;
			os.imbue(std::locale::classic());
			os << std::fixed << std::setprecision(precision) << std::setw(width) << value;
			std::string text = os.str();
			if (static_cast<int>(text.size()) <= width &&
			    (value == 0 || std::abs(value) >= 0.5 * std::pow(10.0, -precision)))
				return text;
			for (int digits = std::max(0, width - 6); digits >= 0; --digits) {
				os.str({});
				os.clear();
				os << std::uppercase << std::scientific << std::setprecision(digits) << std::setw(width) << value;
				if (static_cast<int>(os.str().size()) <= width)
					return os.str();
			}
			throw std::invalid_argument("RWO numeric field cannot fit its fixed width");
		}

		[[nodiscard]] double computePhotometricResiduals(std::vector<WeightedResidualRow>& rows, const FitResult& fit,
		                                                 bool comet, const PhotometryOptions& options)
		{
			const auto solution = photometry::fit(photometry::observations(fit).data(), options.asteroidG);
			const auto cometSolution = comet ? photometry::fitComet(photometry::observations(fit).data(), options)
			                                 : photometry::CometSolution{};
			double squaredResiduals = 0.0;
			double squaredWeights = 0.0;
			for (WeightedResidualRow& row : rows) {
				row.magString = photometry::formatMagnitude(row.obs->MAG, row.obs->magBand, row.obs->magString);
				const auto evaluation = comet ? photometry::evaluateComet(*row.residual, cometSolution,
				                                                          fit.allResiduals.count() != 0 || row.selected)
				                              : photometry::evaluate(*row.residual, solution,
				                                                     fit.allResiduals.count() != 0 || row.selected);
				row.magSigma = evaluation.sigma > 0.0 ? evaluation.sigma : -1.0;
				row.magSelected = evaluation.selected ? 1 : 0;
				row.magResidualDefined = evaluation.residual.has_value();
				row.magResidual = evaluation.residual.value_or(0.0);
				if (comet && evaluation.selected && evaluation.residual && evaluation.sigma > 0.0) {
					const double weight = 1.0 / evaluation.sigma;
					squaredResiduals += std::pow(*evaluation.residual * weight, 2);
					squaredWeights += weight * weight;
				}
			}
			return comet ? (squaredWeights > 0.0 ? std::sqrt(squaredResiduals / squaredWeights) : -1.0) : solution.rms;
		}

		[[nodiscard]] std::string magnitudeBlockForRwo(const WeightedResidualRow& row)
		{
			if (row.magString.find_first_not_of(' ') == std::string::npos)
				return std::string(19, ' ');

			std::string text(19, ' ');
			text.replace(0, 6, row.magString.substr(0, 6));
			if (row.magSigma >= 0.0)
				text.replace(7, 5, formatFixedField(row.magSigma, 5, 2));
			if (row.magResidualDefined)
				text.replace(13, 6, formatFixedField(row.magResidual, 6, 2));
			return text;
		}

		void writeWeightedResidualsRwo(const ExportContext& context, const std::filesystem::path& dir,
		                               ExportReport& report)
		{
			const auto file = dir / "weighted_residuals_rwo.txt";
			const FitResult& fit = *context.fitResult;
			const std::string name = fixedWidthRwoObjectName(Helpers::unpackMpcPackedDesignation(context.objectName));

			std::unordered_set<std::size_t> kept;
			kept.insert(fit.keptObsIndices.begin(), fit.keptObsIndices.end());
			const bool allKept = kept.empty() && fit.residuals.count() == context.observations.size();
			for (auto index : kept)
				if (index >= context.observations.size())
					throw std::invalid_argument("RWO selection index out of range");

			std::vector<WeightedResidualRow> rows;
			rows.reserve(context.observations.size());
			double astSumSq = 0.0;
			std::size_t astSelectedCount = 0;

			for (std::size_t i = 0; i < context.observations.size(); ++i) {
				const ObsNode& obs = context.observations[i];
				if (obs.radar)
					continue;
				const Residual* residual = residualForObservation(fit, obs, i);
				if (!residual || !residual->astrometricWeight)
					throw std::runtime_error("RWO export requires final-fit residual weights");
				const auto& weight = *residual->astrometricWeight;
				if (!(std::isfinite(weight.sigmaRaArcsec) && weight.sigmaRaArcsec > 0.0 &&
				      std::isfinite(weight.sigmaDecArcsec) && weight.sigmaDecArcsec > 0.0 &&
				      std::isfinite(weight.chiSquared) && weight.chiSquared >= 0.0) ||
				    (weight.postFitChiSquared &&
				     (!std::isfinite(*weight.postFitChiSquared) || *weight.postFitChiSquared < 0.0)))
					throw std::runtime_error("RWO export received invalid final-fit residual weights");

				WeightedResidualRow row;
				row.obs = &obs;
				row.residual = residual;
				row.draArcsec = residual ? residual->dRA_arcsec : 0.0;
				row.ddecArcsec = residual ? residual->dDEC_arcsec : 0.0;
				row.sigmaRaArcsec = weight.sigmaRaArcsec;
				row.sigmaDecArcsec = weight.sigmaDecArcsec;
				row.forceWeightRa = weight.forceWeightRa;
				row.forceWeightDec = weight.forceWeightDec;
				row.selected = allKept || kept.contains(i);

				row.chi = std::sqrt(weight.postFitChiSquared.value_or(weight.chiSquared));

				if (row.selected && residual) {
					astSumSq += weight.chiSquared;
					astSelectedCount += 2;
				}

				rows.push_back(row);
			}

			const double rmsMag = computePhotometricResiduals(rows, fit, context.cometPhotometry, context.photometry);

			for (std::size_t i = 0; i < context.observations.size(); ++i) {
				const auto& obs = context.observations[i];
				if (!obs.radar || (!allKept && !kept.contains(i)))
					continue;
				const auto* residual = residualForObservation(fit, obs, i);
				if (!residual || !residual->radar)
					throw std::runtime_error("RWO export requires radar residuals");
				astSumSq += residual->radar->chiSquared;
				astSelectedCount += obs.radar->dimensions();
			}

			const double rmsAst = astSelectedCount > 0 ? std::sqrt(astSumSq / astSelectedCount) : 0.0;

			auto os = openTextFile(file);
			os << "version =   2\n";
			if (context.cometPhotometry)
				os << "! TRAJECTA_PHOTOMETRY comet K1=" << context.photometry.cometK1
				   << " K2=" << context.photometry.cometK2
				   << " PHCOF=" << context.photometry.cometNuclearPhaseCoefficient << '\n';
			os << "! OC_OBJECT " << nlohmann::json(context.objectName).dump() << '\n';
			os << "errmod  = '" << errmodName(context.weightProfile) << "'\n";
			os << std::uppercase << std::scientific << std::setprecision(5);
			os << "RMSast  = " << std::setw(13) << rmsAst << '\n';
			if (rmsMag >= 0.0)
				os << "RMSmag  = " << std::setw(13) << rmsMag << '\n';
			os << std::defaultfloat << std::setprecision(17);
			os << "END_OF_HEADER\n";
			os << "! Object   Obser ============= Date ============= ================== Right Ascension "
			      "=================  "
			      "================= Declination ===================== ==== Magnitude ==== Ast Obs  Residual SEL\n";
			os << "! Design   K T N YYYY MM DD.dddddddddd   Accuracy HH MM SS.sss  Accuracy      RMS  F     Bias    "
			      "Resid "
			      "sDD MM "
			      "SS.ss  Accuracy      RMS  F     Bias    Resid Val  B   RMS  Resid Cat Cod       Chi A M\n";

			for (const WeightedResidualRow& row : rows) {
				const ObsNode& obs = *row.obs;
				os << "! OC_INPUT " << obs::rwoInputMetadata(obs, context.gaiaPhotocenterRadiusKm).dump() << '\n';
				os << ' ' << std::left << std::setw(9) << name << std::right;
				if (obs.rovingObserver && (obs.isSatellite || !obs.rovingObserver->valid()))
					throw std::invalid_argument("Invalid roving observer geometry for RWO");
				if (!obs.rovingObserver && (obs.note2 == 'V' || obs.obsCode == "247" || obs.obsCode == "270"))
					throw std::invalid_argument("RWO roving observer requires original terrestrial coordinates");
				const char tech = obs.rovingObserver ? 'V' : obs.isSatellite ? 'S' : obs.note2 == ' ' ? 'C' : obs.note2;
				os << (obs.isSatellite ? " S " : " O ") << tech << ' ';
				if (obs.note1 != ' ')
					os << obs.note1;
				else
					os << ' ';
				const double timeAccuracyDays = positiveOr(obs.timeAccuracyDays, 1.0e-6);
				const double raAccuracyArcsec = positiveOr(obs.raAccuracyArcsec, 1.5e-2);
				const double decAccuracyArcsec = positiveOr(obs.decAccuracyArcsec, 1.0e-2);

				const double rawRa = obs.submittedRaRad();
				const double rawDec = obs.submittedDecRad();

				const double rwoJd = obs::rwoEpochJd(obs.JD, obs.timeScale);
				os << ' ' << dateFieldsForRwo(rwoJd, timeAccuracyDays);
				os << ' ' << std::uppercase << std::scientific << std::setprecision(3) << std::setw(10)
				   << timeAccuracyDays << std::defaultfloat << ' ';
				os << raFields(rawRa);
				os << ' ' << std::uppercase << std::scientific << std::setprecision(3) << std::setw(10)
				   << raAccuracyArcsec << std::defaultfloat << ' ';
				os << formatFixedField(row.sigmaRaArcsec, 8, 3) << ' ' << (row.forceWeightRa ? 'T' : 'F') << ' '
				   << formatFixedField(obs.biasRaArcsec * std::cos(rawDec), 8, 3)
				   << formatFixedField(row.draArcsec, 9, 3) << ' ';
				os << decFields(rawDec);
				os << ' ' << std::uppercase << std::scientific << std::setprecision(3) << std::setw(10)
				   << decAccuracyArcsec << std::defaultfloat << ' ';
				os << formatFixedField(row.sigmaDecArcsec, 8, 3) << ' ' << (row.forceWeightDec ? 'T' : 'F') << ' '
				   << formatFixedField(obs.biasDecArcsec, 8, 3) << formatFixedField(row.ddecArcsec, 9, 3) << ' ';
				os << magnitudeBlockForRwo(row) << "   ";
				os << (obs.catalogCode == ' ' ? ' ' : obs.catalogCode) << ' ' << std::setw(3) << obs.obsCode << "    "
				   << formatFixedField(row.chi, 6, 2) << ' ' << (row.selected ? 1 : 0) << ' ' << row.magSelected
				   << '\n';
				os << std::defaultfloat << std::setprecision(17);
				if (obs.isSatellite) {
					const auto eq = Helpers::eclToEcv(obs.satelliteOffsetEclAu);
					std::ostringstream prefix;
					prefix << ' ' << std::left << std::setw(9) << name << std::right << " S s " << obs.note1 << ' '
					       << dateFieldsForRwo(rwoJd, timeAccuracyDays);
					os << prefix.str().substr(0, 33) << " 2" << std::fixed << std::setprecision(16) << std::setw(24)
					   << eq.x << std::setw(24) << eq.y << std::setw(24) << eq.z << ' ' << obs.obsCode << '\n'
					   << std::defaultfloat << std::setprecision(17);
				} else if (obs.rovingObserver) {
					const auto& position = *obs.rovingObserver;
					const double longitude =
					    position.longitudeDeg < 0.0 ? position.longitudeDeg + 360.0 : position.longitudeDeg;
					std::ostringstream prefix;
					prefix << ' ' << std::left << std::setw(9) << name << std::right << " O v " << obs.note1 << ' '
					       << dateFieldsForRwo(rwoJd, timeAccuracyDays);
					os << prefix.str().substr(0, 33) << ' ' << formatFixedField(longitude, 10, 6) << ' '
					   << formatFixedField(position.latitudeDeg, 10, 6) << ' '
					   << formatFixedField(position.heightMeters, 8, 1) << ' ' << obs.obsCode << '\n';
				}
			}

			report.files.push_back(file);
			if (std::ranges::any_of(context.observations, [](const auto& obs) { return obs.radar.has_value(); })) {
				os << "! Object   Obser ====== Date =======  ============ Radar range/range rate (km or km/d) "
				      "============= "
				      "Station    Residual\n";
				os << "! Design   K T N YYYY MM DD hh:mm:ss        Measure     Accuracy    rms    F      Bias       "
				      "Resid TRX "
				      "RCX  "
				      "   Chi   S\n";
			}
			for (std::size_t i = 0; i < context.observations.size(); ++i) {
				const auto& obs = context.observations[i];
				if (!obs.radar)
					continue;
				const auto* residual = residualForObservation(fit, obs, i);
				if (!residual || !residual->radar)
					throw std::runtime_error("Radar export requires final-fit residuals");
				int year, month;
				double day;
				Convertors::jdToYMD(obs::rwoEpochJd(obs.JD, obs.timeScale) + 0.5 / Constants::SECONDS_PER_DAY, year,
				                    month, day);
				const int seconds = static_cast<int>(std::floor((day - std::floor(day)) * 86400.0));
				const auto write = [&](const auto& channel, char type, double factor) {
					if (!channel)
						return;
					os << "! OC_RADAR " << std::setprecision(17) << obs.radar->frequencyHz << ' '
					   << obs.radar->frequencyIsReference << ' ' << obs.radar->surfaceRadiusKm.value_or(0.0) << ' '
					   << channel->observed << ' ' << channel->sigma << ' ' << obs.JD << ' ' << obs.hasSplitUtc << ' '
					   << obs.utcMjdDay << ' ' << obs.utcSecondsOfDay;
					if (obs.timeScale != time::legacyObservationScale(obs.JD))
						os << ' ' << static_cast<int>(obs.timeScale);
					os << '\n';
					os << ' ' << std::left << std::setw(9) << name << std::right << ' ' << type << ' '
					   << (obs.radar->returnType == obs::RadarReturn::Surface ? 's' : 'c') << "   " << std::setfill('0')
					   << std::setw(4) << year << ' ' << std::setw(2) << month << ' ' << std::setw(2)
					   << static_cast<int>(day) << ' ' << std::setw(2) << seconds / 3600 << ':' << std::setw(2)
					   << seconds / 60 % 60 << ':' << std::setw(2) << seconds % 60 << std::setfill(' ')
					   << formatFixedField(channel->observed * factor, 16, 5) << ' '
					   << formatFixedField(channel->sigma * std::abs(factor), 9, 5) << ' '
					   << formatFixedField(channel->sigma * std::abs(factor), 9, 5) << " F "
					   << formatFixedField(0.0, 11, 5) << ' ' << formatFixedField(channel->residual * factor, 11, 5)
					   << ' ' << obs.radar->transmitter << ' ' << obs.obsCode << ' '
					   << formatFixedField(std::abs(channel->residual / channel->sigma), 9, 2) << ' '
					   << (allKept || kept.contains(i)) << '\n'
					   << std::defaultfloat << std::setprecision(17);
				};
				constexpr double cKmSec = Constants::AU_PER_DAY * Constants::AU_IN_METERS / 1000.0 / 86400.0;
				write(residual->radar->delaySeconds, 'R', cKmSec / 2.0);
				write(residual->radar->dopplerHz, 'V', -cKmSec * 86400.0 / (2.0 * obs.radar->frequencyHz));
			}
		}

		template <class Matrix>
		void writeMatrix(std::ostream& os, const char* label, const Matrix& matrix)
		{
			os << '[' << label << "]\n";
			os << "rows " << matrix.rows() << " cols " << matrix.cols() << '\n';
			for (Eigen::Index r = 0; r < matrix.rows(); ++r) {
				for (Eigen::Index c = 0; c < matrix.cols(); ++c) {
					if (c > 0)
						os << ' ';
					os << matrix(r, c);
				}
				os << '\n';
			}
			os << '\n';
		}

		void writeNonGrav(std::ostream& os, const NonGravModel& ng)
		{
			os << "[non_grav]\n";
			os << "enabled " << (ng.enabled() ? 1 : 0) << '\n';
			os << "model " << static_cast<int>(ng.model()) << '\n';
			os << "fit_mask " << static_cast<int>(static_cast<std::uint8_t>(ng.fitParameterMask())) << '\n';
			os << "A1_au_day2 " << ng.A1() << '\n';
			os << "A2_au_day2 " << ng.A2() << '\n';
			os << "A3_au_day2 " << ng.A3() << '\n';
			os << "DT_day " << ng.DT() << "\n\n";
		}

		void writeSummary(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			const auto file = dir / "fit_summary.txt";
			auto os = openTextFile(file);
			const FitResult& fit = *context.fitResult;

			os << "object " << sanitizePathPart(context.objectName) << '\n';
			os << "generated_local_date " << localDateStamp() << '\n';
			os << "converged " << (fit.converged ? 1 : 0) << '\n';
			os << "cancelled " << (fit.cancelled ? 1 : 0) << '\n';
			os << "iterations " << fit.iterations << '\n';
			os << "cost " << fit.cost << '\n';
			os << "rms_arcsec " << fit.residuals.rms() << '\n';
			os << "residual_rows " << fit.residuals.count() << '\n';
			os << "input_observations " << context.observations.size() << '\n';
			os << "kept_observations " << fit.keptObsIndices.size() << '\n';
			os << "sigma_clip_passes " << fit.sigmaClipPassesRun << '\n';
			os << "sigma_clip_rejected " << fit.sigmaClipRejected << '\n';
			os << "covariance_objective_scale " << fit.uncertainty.covarianceObjectiveScale << "\n\n";
			if (context.cometPhotometry) {
				const auto photometry = orbital_core::photometry::fitComet(
				    orbital_core::photometry::observations(fit).data(), context.photometry);
				os << "[photometry_comet]\n";
				if (photometry.total.absoluteMagnitude)
					os << "M1 " << *photometry.total.absoluteMagnitude << '\n';
				if (photometry.nuclear.absoluteMagnitude)
					os << "M2 " << *photometry.nuclear.absoluteMagnitude << '\n';
				os << "K1_fixed " << photometry.total.slope << '\n'
				   << "K2_fixed " << photometry.nuclear.slope << '\n'
				   << "PHCOF_fixed_mag_per_degree " << photometry.nuclear.phaseCoefficient << '\n'
				   << "total_observations " << photometry.total.count << '\n'
				   << "nuclear_observations " << photometry.nuclear.count << '\n'
				   << "unknown_type_observations " << photometry.unknownTypeCount << "\n\n";
			} else {
				const auto photometry = orbital_core::photometry::fit(
				    orbital_core::photometry::observations(fit).data(), context.photometry.asteroidG);
				os << "[photometry_asteroid]\n";
				if (photometry.absoluteMagnitude)
					os << "H " << *photometry.absoluteMagnitude << '\n';
				os << "G_fixed " << context.photometry.asteroidG << '\n'
				   << "observations " << photometry.count << "\n\n";
			}

			writeNonGrav(os, fit.finalNonGrav);
			report.files.push_back(file);
		}

		void writeOrbitElements(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			const Orbit_p& orbit = exportStateOrbit(*context.fitResult);
			if (orbit && cometary::usesCometaryElements(orbit->getExcentricity())) {
				writeCartesianOef(context, dir, report);
				writeCometaryOef(context, dir, report);
				return;
			}
			writeEquinoctialOef(context, dir, report);
			writeKeplerianOef(context, dir, report);
		}

		void writeStateVectors(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			writeCartesianOef(context, dir, report);
		}

		void writeCovariance(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			const auto file = dir / "covariance.txt";
			auto os = openTextFile(file);
			const auto& unc = context.fitResult->uncertainty;

			os << "has_element_covariance " << (unc.hasElementCovariance ? 1 : 0) << '\n';
			os << "covariance_objective_scale " << unc.covarianceObjectiveScale << "\n\n";

			writeMatrix(os, "element_covariance_a_e_i_node_peri_nu", unc.elementCovariance);
			writeMatrix(os, "cartesian_state_covariance_fit_epoch_R_AU_V_AU_day", unc.cartesianStateCovariance);
			writeMatrix(os, "cartesian_state_covariance_formal_fit_epoch_R_AU_V_AU_day",
			            unc.cartesianStateCovarianceFormal);
			writeMatrix(os, "fit_parameter_covariance", unc.fitParameterCovariance);
			writeMatrix(os, "fit_parameter_covariance_formal", unc.fitParameterCovarianceFormal);

			if (unc.cometaryElementCovariance) {
				const auto& cot = *unc.cometaryElementCovariance;
				os << "cometary_element_covariance_epoch_jd_tdb " << cot.epochJdTdb << '\n';
				os << "cometary_element_covariance_units AU 1 rad rad rad rad\n\n";
				writeMatrix(os, "cometary_element_covariance_q_e_i_node_peri_nu", cot.covariance);
				writeMatrix(os, "cometary_element_covariance_formal_q_e_i_node_peri_nu", cot.formalCovariance);
			}

			report.files.push_back(file);
		}

		void writeResiduals(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			writeWeightedResidualsRwo(context, dir, report);
		}

		void writeObservations(const ExportContext& context, const std::filesystem::path& dir, ExportReport& report)
		{
			writeWeightedResidualsRwo(context, dir, report);
		}

		[[nodiscard]] bool includesAll(std::span<const ExportProduct> products)
		{
			return std::find(products.begin(), products.end(), ExportProduct::All) != products.end();
		}

		[[nodiscard]] std::vector<ExportProduct> expandProducts(std::span<const ExportProduct> products)
		{
			if (includesAll(products)) {
				return {ExportProduct::OrbitElements, ExportProduct::StateVectors, ExportProduct::Residuals};
			}
			return {products.begin(), products.end()};
		}

	}  // namespace

	ExportReport Exporter::save(const ExportContext& context, std::initializer_list<ExportProduct> products)
	{
		return save(context, std::span<const ExportProduct>(products.begin(), products.size()));
	}

	ExportReport Exporter::save(const ExportContext& context, std::span<const ExportProduct> products)
	{
		requireFit(context);

		ExportReport report;
		report.directory = resolveDirectory(context);

		const std::vector<ExportProduct> expanded = expandProducts(products);
		for (const ExportProduct product : expanded) {
			switch (product) {
				case ExportProduct::FitSummary:
					writeSummary(context, report.directory, report);
					break;
				case ExportProduct::OrbitElements:
					writeOrbitElements(context, report.directory, report);
					break;
				case ExportProduct::StateVectors:
					writeStateVectors(context, report.directory, report);
					break;
				case ExportProduct::Covariance:
					writeCovariance(context, report.directory, report);
					break;
				case ExportProduct::Residuals:
					writeResiduals(context, report.directory, report);
					break;
				case ExportProduct::Observations:
					writeObservations(context, report.directory, report);
					break;
				case ExportProduct::All:
					break;
			}
		}

		return report;
	}

}  // namespace orbital_core
