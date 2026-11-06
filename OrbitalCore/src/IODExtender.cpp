

#include "IODExtender.h"
#include "fit/OppositionCoverage.h"
#include <utility>

#include "IodScoring.h"
#include "NBodyPropagator.h"
#include "PerturberSelector.h"
#include "logger/Logger.h"
#include "fit/FitObservations.h"
#include "obs/AstrometricResidual.h"
#include "obs/LightDeflection.h"
#include "obs/LightTime.h"
#include "obs/ObserverGeometry.h"
#include "obs/PreparedObservation.h"
#include "orbital_core/Constants.h"
#include "orbital_core/FitConfig.h"
#include "orbital_core/Types.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace Calculations {

	namespace {


		constexpr double kResidualLimitSigma = 200.0;
		constexpr double kWindowMaxSigma = 2.0;

		constexpr double kArcLimitDays = 7300.0;


		constexpr int kMaxOuterIter = 30;
		constexpr int kMaxForceModelFits = 5;
		constexpr int kMaxAutoRejectRounds = 4;
		constexpr double kAutoRejectSigma = 3.0;
		constexpr double kAutoRejectOuterScale = 1.1;
		constexpr double kAutoRejectRetryScale = 1.5;
		constexpr double kAutoRejectMaxSigma = 60.0;
		constexpr double kLovDelta = 0.0001;
		constexpr int kLovDivisions = 10000;

		struct ResidualEval {
			std::vector<double> totalArcsec;
			std::vector<double> totalSigma;
			double rmsAll{1.0e30};
			uint16_t autoPerturberMask{0};
		};

		struct ProjectedRaDec {
			double ra{0.0};
			double dec{0.0};
			double rho{0.0};
			bool valid{false};
		};

		[[nodiscard]] ProjectedRaDec projectRaDec(const Vec3D& helioPosition,
		                                          const orbital_core::obs::PreparedObservation& po, double cosObl,
		                                          double sinObl) noexcept
		{
			const Vec3D R_topo = orbital_core::obs::opticalObservationDirection(helioPosition, po);
			Vec3D R_eq(R_topo.x, R_topo.y * cosObl - R_topo.z * sinObl, R_topo.y * sinObl + R_topo.z * cosObl);

			const double rho = R_eq.magnitude();
			if (rho <= 0.0 || !std::isfinite(rho))
				return {};

			double ra = std::atan2(R_eq.y, R_eq.x);
			if (ra < 0.0)
				ra += 2.0 * Constants::PI;

			return ProjectedRaDec{ra, std::asin(std::clamp(R_eq.z / rho, -1.0, 1.0)), rho, true};
		}

		[[nodiscard]] ResidualEval evaluateResidualsArcsec(
		    const Orbit_p& orbit, std::span<const orbital_core::obs::PreparedObservation> prepared,
		    NBodyPropagator& nbody, const ArcScheduledPerturberSelector* perturberSelector = nullptr)
		{
			const Vec3D R0 = orbit->getR();
			const Vec3D V0_real = orbit->getV() * Constants::GAUSS_K;
			const double epochTdb = orbit->getEpoch();
			const double cosObl = std::cos(Constants::EARTH_INCL);
			const double sinObl = std::sin(Constants::EARTH_INCL);
			const auto& ng = nbody.nonGrav();

			const std::size_t N = prepared.size();
			ResidualEval out;
			out.totalArcsec.assign(N, 1.0e30);
			out.totalSigma.assign(N, 1.0e30);

			std::vector<double> jdsTdb(N);
			for (std::size_t i = 0; i < N; ++i)
				jdsTdb[i] = prepared[i].jdTdb;

			const auto receptionStates = nbody.propagateAllWithVelocity(R0, V0_real, epochTdb, jdsTdb, ng);
			if (receptionStates.size() != N)
				return out;

			std::vector<NBodyPropagator::StateVel> states(N);
			std::vector<Vec3D> positions(N);
			for (std::size_t i = 0; i < N; ++i) {
				const auto apparent =
				    orbital_core::obs::solveLightTimeTwoPass(jdsTdb[i], receptionStates[i].R, receptionStates[i].V,
				                                             prepared[i].earthR + prepared[i].observerOffsetEclAu);
				if (!apparent.valid)
					return out;
				states[i] = {apparent.position, apparent.velocityAuPerDay};
				positions[i] = apparent.position;
			}

			bool hasTimingSigma = false;
			for (const auto& po : prepared) {
				if (po.timeSigmaDays > 0.0) {
					hasTimingSigma = true;
					break;
				}
			}

			if (perturberSelector)
				out.autoPerturberMask = perturberSelector->detectAutomaticMask(positions, prepared);

			double sumSq = 0.0;
			std::size_t nFinite = 0;
			for (std::size_t i = 0; i < N; ++i) {
				const auto& po = prepared[i];
				const ProjectedRaDec predicted = projectRaDec(positions[i], po, cosObl, sinObl);
				if (!predicted.valid) {
					out.totalArcsec[i] = 1.0e30;
					out.totalSigma[i] = 1.0e30;
					continue;
				}

				const double dRA = Helpers::wrapPi(po.raObs - predicted.ra);
				const double dDEC = Helpers::wrapPi(po.decObs - predicted.dec);
				const double cosDecObs = po.raMetricCosDec;
				const double dRA_arc = (cosDecObs * dRA) / Constants::ARCSEC_RAD;
				const double dDEC_arc = dDEC / Constants::ARCSEC_RAD;
				const double total = std::sqrt(dRA_arc * dRA_arc + dDEC_arc * dDEC_arc);

				const Vec3D* velocity = hasTimingSigma ? &states[i].V : nullptr;
				const auto transform = orbital_core::obs::astrometricResidualTransform(po, positions[i], velocity);
				const auto normalized =
				    orbital_core::obs::applyAstrometricResidualTransform(cosDecObs * dRA, dDEC, transform);

				if (std::isfinite(total)) {
					out.totalArcsec[i] = total;
					out.totalSigma[i] = normalized.valid ? normalized.norm : 1.0e30;
					sumSq += total * total;
					++nFinite;
				} else {
					out.totalArcsec[i] = 1.0e30;
					out.totalSigma[i] = 1.0e30;
				}
			}

			out.rmsAll = (nFinite > 0) ? std::sqrt(sumSq / static_cast<double>(nFinite)) : 1.0e30;
			return out;
		}

		[[nodiscard]] std::size_t countIncluded(const std::vector<std::uint8_t>& included) noexcept
		{
			std::size_t n = 0;
			for (const auto v : included)
				if (v)
					++n;
			return n;
		}

		[[nodiscard]] bool firstLastIncluded(const std::vector<std::uint8_t>& included, std::size_t& firstIdx,
		                                     std::size_t& lastIdx) noexcept
		{
			firstIdx = included.size();
			lastIdx = 0;
			for (std::size_t i = 0; i < included.size(); ++i) {
				if (!included[i])
					continue;
				if (firstIdx == included.size())
					firstIdx = i;
				lastIdx = i;
			}
			return firstIdx != included.size();
		}

		[[nodiscard]] std::vector<ObsNode> buildKeptFromMask(const std::vector<ObsNode>& allObs,
		                                                     const std::vector<std::uint8_t>& included)
		{
			std::vector<ObsNode> kept;
			kept.reserve(countIncluded(included));
			for (std::size_t i = 0; i < allObs.size(); ++i)
				if (included[i])
					kept.push_back(allObs[i]);
			return kept;
		}

		[[nodiscard]] std::vector<std::size_t> buildKeptIndices(const std::vector<std::uint8_t>& included)
		{
			std::vector<std::size_t> idx;
			idx.reserve(countIncluded(included));
			for (std::size_t i = 0; i < included.size(); ++i)
				if (included[i])
					idx.push_back(i);
			return idx;
		}

		[[nodiscard]] double rmsArcsecFromMask(const std::vector<double>& totalArcsec,
		                                       const std::vector<std::uint8_t>& included) noexcept
		{
			double sumSq = 0.0;
			std::size_t n = 0;
			const std::size_t count = std::min(totalArcsec.size(), included.size());
			for (std::size_t i = 0; i < count; ++i) {
				if (!included[i])
					continue;
				const double r = totalArcsec[i];
				if (std::isfinite(r) && r < 1.0e29) {
					sumSq += r * r;
					++n;
				}
			}

			return n > 0 ? std::sqrt(sumSq / static_cast<double>(n)) : 1.0e30;
		}

		[[nodiscard]] double residualValueAt(const std::vector<double>& values, std::size_t idx) noexcept
		{
			if (idx >= values.size())
				return std::numeric_limits<double>::infinity();
			return values[idx];
		}

		using LovVector6 = std::array<double, 6>;

		[[nodiscard]] std::array<double, 3> rotateObsVectorFO(const ObsNode& obs, const Vec3D& vector) noexcept
		{
			const double r = std::sqrt(obs.vector.x * obs.vector.x + obs.vector.y * obs.vector.y);
			if (r <= 0.0 || !std::isfinite(r))
				return {0.0, 0.0, 1.0};

			const Vec3D sideways(obs.vector.y / r, -obs.vector.x / r, 0.0);
			const Vec3D cross = sideways.cross(obs.vector);

			return {sideways.dot(vector), cross.dot(vector), obs.vector.dot(vector)};
		}

		[[nodiscard]] double pseudoResidFO(double y, double x) noexcept
		{
			if (y < x && y > -x)
				return y / x;
			return 2.0 - x / std::abs(y);
		}

		[[nodiscard]] double inverseErfFO(double y) noexcept
		{
			double x = 0.0;
			double diff = 0.0;
			constexpr double kErfLimit = 0.915;
			constexpr double kSqrtPi = 1.7724538509055160272981674833411451827975494561223871282138077898529113;

			if (y < -kErfLimit)
				return -inverseErfFO(-y);
			if (y > kErfLimit)
				x = std::sqrt(-std::log(1.0 - y)) - 0.34;
			else
				x = y * (0.8963 + y * y * (0.0889 + 0.494 * y * y));

			int nIterations = 0;
			do {
				const double dy = std::erf(x) - y;
				const double slope = (2.0 / kSqrtPi) * std::exp(-x * x);
				diff = -dy / (slope + x * dy);
				x += diff;
			} while (nIterations++ < 10 && std::abs(diff) > 1.0e-14);

			return x;
		}

		[[nodiscard]] double parabolicMinimumPointFO(const double x[3], const double y[3]) noexcept
		{
			const double x1 = x[1] - x[0];
			const double x2 = x[2] - x[0];
			const double y1 = y[1] - y[0];
			const double y2 = y[2] - y[0];
			const double bottom = y1 * x2 - y2 * x1;
			const double top = y1 * x2 * x2 - y2 * x1 * x1;
			return x[0] + top / (2.0 * bottom);
		}

		[[nodiscard]] double lovSearchScoreFO(const std::vector<std::array<double, 3>>& xyz,
		                                      const std::vector<std::array<double, 3>>& slopes, double sigmas) noexcept
		{
			double score = 0.0;
			const std::size_t n = std::min(xyz.size(), slopes.size());
			for (std::size_t i = 0; i < n; ++i) {
				const double x = xyz[i][0] + sigmas * slopes[i][0];
				const double y = xyz[i][1] + sigmas * slopes[i][1];
				const double z = xyz[i][2] + sigmas * slopes[i][2];
				const double dx = pseudoResidFO(x, z);
				const double dy = pseudoResidFO(y, z);
				score += dx * dx + dy * dy;
			}
			return score;
		}

		[[nodiscard]] double lovScoreMinimumFO(const std::vector<std::array<double, 3>>& xyz,
		                                       const std::vector<std::array<double, 3>>& slopes, const double startX[3],
		                                       const double startY[3], double& xmin) noexcept
		{
			double x[3] = {startX[0], startX[1], startX[2]};
			double y[3] = {startY[0], startY[1], startY[2]};

			for (int iter = 0; iter < 5; ++iter) {
				const double newX = parabolicMinimumPointFO(x, y);
				x[0] = x[1];
				y[0] = y[1];
				x[1] = x[2];
				y[1] = y[2];
				x[2] = newX;
				y[2] = lovSearchScoreFO(xyz, slopes, newX);
			}

			xmin = x[2];
			return y[2];
		}

		[[nodiscard]] std::vector<double> stateVectorFromOrbitFO(const Orbit& orbit)
		{
			const Vec3D R = orbit.getR();
			const Vec3D V = orbit.getV();
			return {R.x, R.y, R.z, V.x, V.y, V.z};
		}

		[[nodiscard]] Orbit_p orbitFromStateStepFO(const Orbit& orbit, const LovVector6& lov, double sigmas)
		{
			std::vector<double> state = stateVectorFromOrbitFO(orbit);
			for (std::size_t i = 0; i < lov.size(); ++i)
				state[i] += sigmas * lov[i];

			return make_orbit_from_states_vectors(
			    StateVectors(Vec3D(state[0], state[1], state[2]), Vec3D(state[3], state[4], state[5])),
			    orbit.getEpoch());
		}

		[[nodiscard]] std::vector<Vec3D> topocentricVectorsForLovFO(
		    const Orbit& orbit, std::span<const orbital_core::obs::PreparedObservation> prepared,
		    NBodyPropagator& nbody)
		{
			const std::size_t n = prepared.size();
			std::vector<double> jdsTdb(n);
			for (std::size_t i = 0; i < n; ++i)
				jdsTdb[i] = prepared[i].jdTdb;

			const Vec3D R0 = orbit.getR();
			const Vec3D V0Real = orbit.getV() * Constants::GAUSS_K;
			const double epochTdb = orbit.getEpoch();

			const auto receptionStates = nbody.propagateAllWithVelocity(R0, V0Real, epochTdb, jdsTdb, nbody.nonGrav());
			if (receptionStates.size() != n)
				return {};

			std::vector<Vec3D> topo(n);
			for (std::size_t i = 0; i < n; ++i) {
				const Vec3D observer = prepared[i].earthR + prepared[i].observerOffsetEclAu;
				const auto apparent = orbital_core::obs::solveLightTimeTwoPass(jdsTdb[i], receptionStates[i].R,
				                                                               receptionStates[i].V, observer);
				if (!apparent.valid)
					return {};
				topo[i] = apparent.position - observer;
			}

			return topo;
		}

		[[nodiscard]] bool oneSigmaLovFromCovarianceFO(const orbital_core::StateCov6x6& stateCov,
		                                               LovVector6& lov) noexcept
		{
			Eigen::Matrix<double, 6, 6> cov = Eigen::Matrix<double, 6, 6>::Zero();
			Eigen::Matrix<double, 6, 6> toCanonical = Eigen::Matrix<double, 6, 6>::Identity();
			toCanonical(3, 3) = 1.0 / Constants::GAUSS_K;
			toCanonical(4, 4) = 1.0 / Constants::GAUSS_K;
			toCanonical(5, 5) = 1.0 / Constants::GAUSS_K;

			cov = toCanonical * stateCov * toCanonical.transpose();
			cov = 0.5 * (cov + cov.transpose()).eval();

			Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(cov);
			if (solver.info() != Eigen::Success)
				return false;

			constexpr int kLargestCovarianceAxis = 5;
			const double eigenvalue = solver.eigenvalues()(kLargestCovarianceAxis);
			if (!(eigenvalue > 0.0) || !std::isfinite(eigenvalue))
				return false;

			const double sigma = std::sqrt(eigenvalue);
			const auto axis = solver.eigenvectors().col(kLargestCovarianceAxis);
			for (int i = 0; i < 6; ++i)
				lov[static_cast<std::size_t>(i)] = axis(i) * sigma;

			return std::all_of(lov.begin(), lov.end(), [](double v) { return std::isfinite(v); });
		}

		struct LovImproveResult {
			Orbit_p orbit;
			double sigmas{0.0};
			bool applied{false};
		};

		[[nodiscard]] LovImproveResult improveAlongLovFO(
		    const Orbit_p& orbit, const LovVector6& lov, const std::vector<ObsNode>& allObs,
		    std::span<const orbital_core::obs::PreparedObservation> prepared, const std::vector<std::uint8_t>& included,
		    NBodyPropagator& nbody)
		{
			LovImproveResult result;
			result.orbit = orbit;
			if (!orbit || allObs.empty() || prepared.size() != allObs.size() || included.size() != allObs.size()) {
				return result;
			}

			std::size_t first = allObs.size();
			std::size_t last = 0;
			if (!firstLastIncluded(included, first, last))
				return result;

			Orbit_p deltaOrbit = orbitFromStateStepFO(*orbit, lov, kLovDelta);
			if (!deltaOrbit)
				return result;

			const std::vector<Vec3D> nominalTopo = topocentricVectorsForLovFO(*orbit, prepared, nbody);
			const std::vector<Vec3D> deltaTopo = topocentricVectorsForLovFO(*deltaOrbit, prepared, nbody);
			if (nominalTopo.size() != allObs.size() || deltaTopo.size() != allObs.size())
				return result;

			std::vector<std::array<double, 3>> xyz;
			std::vector<std::array<double, 3>> slopes;
			xyz.reserve(last - first + 1);
			slopes.reserve(last - first + 1);

			for (std::size_t i = first; i <= last; ++i) {
				if (included[i] && allObs[i].note2 != 'R') {
					std::array<double, 3> loc = rotateObsVectorFO(allObs[i], nominalTopo[i]);
					std::array<double, 3> slopeVec = {(deltaTopo[i].x - nominalTopo[i].x) / kLovDelta,
					                                  (deltaTopo[i].y - nominalTopo[i].y) / kLovDelta,
					                                  (deltaTopo[i].z - nominalTopo[i].z) / kLovDelta};
					std::array<double, 3> slope =
					    rotateObsVectorFO(allObs[i], Vec3D(slopeVec[0], slopeVec[1], slopeVec[2]));

					double sigmaArcsec = std::hypot(prepared[i].sigmaDec, prepared[i].sigmaRa) / Constants::ARCSEC_RAD;
					if (!(sigmaArcsec > 0.0) || !std::isfinite(sigmaArcsec))
						sigmaArcsec = std::sqrt(2.0);

					loc[0] /= sigmaArcsec;
					loc[1] /= sigmaArcsec;
					slope[0] /= sigmaArcsec;
					slope[1] /= sigmaArcsec;

					xyz.push_back(loc);
					slopes.push_back(slope);
				} else {
					xyz.push_back({0.0, 0.0, 1.0});
					slopes.push_back({0.0, 0.0, 0.0});
				}
			}

			std::vector<double> x(static_cast<std::size_t>(kLovDivisions));
			std::vector<double> score(static_cast<std::size_t>(kLovDivisions));
			for (int i = 0; i < kLovDivisions; ++i) {
				x[static_cast<std::size_t>(i)] =
				    inverseErfFO(2.0 * (static_cast<double>(i) + 0.5) / static_cast<double>(kLovDivisions) - 1.0) * 3.0;
				score[static_cast<std::size_t>(i)] = lovSearchScoreFO(xyz, slopes, x[static_cast<std::size_t>(i)]);
			}

			int edgeIter = 0;
			while (score[0] < score[1] && score[0] < score[2] && edgeIter++ < 20) {
				x[0] += x[0] - x[1];
				score[0] = lovSearchScoreFO(xyz, slopes, x[0]);
			}

			edgeIter = 0;
			const std::size_t lastDiv = static_cast<std::size_t>(kLovDivisions - 1);
			while (score[lastDiv] < score[lastDiv - 1] && score[lastDiv] < score[lastDiv - 2] && edgeIter++ < 20) {
				x[lastDiv] += x[lastDiv] - x[lastDiv - 1];
				score[lastDiv] = lovSearchScoreFO(xyz, slopes, x[lastDiv]);
			}

			double bestSigma = 0.0;
			double lowestScore = 1.0e200;
			for (int i = 0; i < kLovDivisions - 2; ++i) {
				const std::size_t idx = static_cast<std::size_t>(i);
				if (score[idx + 1] < score[idx] && score[idx + 1] < score[idx + 2]) {
					double localMin = 0.0;
					const double sx[3] = {x[idx], x[idx + 1], x[idx + 2]};
					const double sy[3] = {score[idx], score[idx + 1], score[idx + 2]};
					const double newScore = lovScoreMinimumFO(xyz, slopes, sx, sy, localMin);
					if (lowestScore > newScore) {
						lowestScore = newScore;
						bestSigma = localMin;
					}
				}
			}

			Orbit_p improved = orbitFromStateStepFO(*orbit, lov, bestSigma);
			if (!improved)
				return result;

			result.orbit = improved;
			result.sigmas = bestSigma;
			result.applied = true;
			return result;
		}

		struct AutoRejectWithinArcResult {
			bool changed{false};
			bool safeguardRollback{false};
			int flippedIn{0};
			int flippedOut{0};
			double appliedSigma{0.0};
			std::size_t keptInArc{0};
		};

		[[nodiscard]] AutoRejectWithinArcResult autoRejectWithinArcFO(std::vector<std::uint8_t>& included,
		                                                              const std::vector<ObsNode>& allObs,
		                                                              const std::vector<double>& residualsSigma,
		                                                              std::size_t arcFirst, std::size_t arcLast)
		{
			AutoRejectWithinArcResult out;
			if (allObs.empty() || included.size() != allObs.size() || arcFirst >= allObs.size() ||
			    arcLast >= allObs.size() || arcFirst > arcLast) {
				return out;
			}

			const std::size_t segmentCount = arcLast - arcFirst + 1;
			const std::size_t minKept = (segmentCount * 2) / 3;
			const double originalSpan = allObs[arcLast].JD - allObs[arcFirst].JD;

			double rejectLimit = kAutoRejectSigma;
			while (rejectLimit < kAutoRejectMaxSigma) {
				const double inner2 = rejectLimit * rejectLimit;
				const double outer = rejectLimit * kAutoRejectOuterScale;
				const double outer2 = outer * outer;

				std::vector<std::uint8_t> originalIncluded = included;
				orbital_core::fit::OppositionCoverage coverage(allObs, included);
				int flippedIn = 0;
				int flippedOut = 0;
				std::size_t nFound = 0;
				std::size_t start = arcLast;
				std::size_t end = arcFirst;

				for (std::size_t i = arcFirst; i <= arcLast; ++i) {
					const double rho = (i < residualsSigma.size() && std::isfinite(residualsSigma[i]))
					                       ? residualsSigma[i]
					                       : std::numeric_limits<double>::infinity();
					const double r2 = rho * rho;

					const bool wasIncluded = (included[i] != 0);
					if (!wasIncluded && r2 < inner2) {
						included[i] = 1;
						coverage.include(i);
						++flippedIn;
					} else if (wasIncluded && r2 > outer2 && !coverage.isLast(i)) {
						included[i] = 0;
						coverage.exclude(i);
						++flippedOut;
					}

					if (included[i]) {
						if (nFound == 0)
							start = i;
						end = i;
						++nFound;
					}
				}

				const double finalSpan = (nFound > 0) ? (allObs[end].JD - allObs[start].JD) : 0.0;

				if (finalSpan >= originalSpan * 0.5 && nFound >= minKept) {
					out.changed = (flippedIn + flippedOut) > 0;
					out.flippedIn = flippedIn;
					out.flippedOut = flippedOut;
					out.appliedSigma = rejectLimit;
					out.keptInArc = nFound;
					return out;
				}

				included = std::move(originalIncluded);
				rejectLimit *= kAutoRejectRetryScale;
			}

			out.safeguardRollback = true;
			out.appliedSigma = rejectLimit;
			return out;
		}

		[[nodiscard]] int extendInclusionMaskFO(std::vector<std::uint8_t>& included, const std::vector<ObsNode>& allObs,
		                                        const std::vector<double>& residualsSigma)
		{
			std::size_t firstIdx = 0;
			std::size_t lastIdx = 0;
			if (allObs.empty() || !firstLastIncluded(included, firstIdx, lastIdx))
				return 0;

			const std::size_t initialCount = countIncluded(included);
			int nAdded = 0;

			auto markIncluded = [&](std::size_t idx) {
				if (!included[idx]) {
					included[idx] = 1;
					++nAdded;
				}
			};

			const double jdLow = allObs[firstIdx].JD * 2.0 - allObs[lastIdx].JD;
			const double jdHigh = allObs[lastIdx].JD * 2.0 - allObs[firstIdx].JD;

			std::size_t scanFirst = firstIdx;
			std::size_t scanLast = lastIdx;

			while (scanFirst > 0) {
				const std::size_t cand = scanFirst - 1;
				if (allObs[scanLast].JD - allObs[cand].JD >= kArcLimitDays)
					break;
				if (residualValueAt(residualsSigma, cand) >= kResidualLimitSigma)
					break;
				markIncluded(cand);
				scanFirst = cand;
			}

			while (scanFirst > 0) {
				const std::size_t cand = scanFirst - 1;
				if (allObs[cand].JD <= jdLow)
					break;
				if (residualValueAt(residualsSigma, cand) < kWindowMaxSigma)
					markIncluded(cand);
				scanFirst = cand;
			}

			while (scanLast + 1 < allObs.size()) {
				const std::size_t cand = scanLast + 1;
				if (allObs[cand].JD - allObs[scanFirst].JD >= kArcLimitDays)
					break;
				if (residualValueAt(residualsSigma, cand) >= kResidualLimitSigma)
					break;
				markIncluded(cand);
				scanLast = cand;
			}

			while (scanLast + 1 < allObs.size()) {
				const std::size_t cand = scanLast + 1;
				if (allObs[cand].JD >= jdHigh)
					break;
				if (residualValueAt(residualsSigma, cand) < kWindowMaxSigma)
					markIncluded(cand);
				scanLast = cand;
			}

			if (nAdded == 0) {
				int direction = 0;
				if (scanFirst == 0) {
					if (scanLast + 1 < allObs.size())
						direction = +1;
				} else if (scanLast + 1 >= allObs.size()) {
					direction = -1;
				} else {
					const double timeBefore = allObs[scanFirst].JD - allObs[scanFirst - 1].JD;
					const double timeAfter = allObs[scanLast + 1].JD - allObs[scanLast].JD;
					direction = (timeBefore > timeAfter) ? +1 : -1;
				}

				if (direction == -1) {
					const std::size_t cand = scanFirst - 1;
					if (allObs[scanLast].JD - allObs[cand].JD < kArcLimitDays)
						markIncluded(cand);
				} else if (direction == +1) {
					const std::size_t cand = scanLast + 1;
					if (allObs[cand].JD - allObs[scanFirst].JD < kArcLimitDays)
						markIncluded(cand);
				}
			}

			const std::size_t finalCount = countIncluded(included);
			return (finalCount == initialCount) ? 0 : static_cast<int>(finalCount - initialCount);
		}

		struct FitSubsetResult {
			Orbit_p orbit;
			Orbit_p fitEpochOrbit;
			orbital_core::StateCov6x6 stateCovariance{orbital_core::StateCov6x6::Zero()};
			bool hasCovariance{false};
			bool converged{false};
		};

		[[nodiscard]] FitSubsetResult fitOnSubset(const Orbit_p& seed,
		                                          std::span<const orbital_core::obs::PreparedObservation> kept,
		                                          const OrbitOptimizerConfig& fitCfg)
		{
			FitSubsetResult out;
			try {
				OrbitOptimizer fitOpt(
				    seed, std::vector<orbital_core::obs::PreparedObservation>(kept.begin(), kept.end()), fitCfg);
				const auto& fit = fitOpt.optimize();
				out.converged = fit.converged;
				out.orbit = fit.orbit;
				out.fitEpochOrbit = fit.fitEpochOrbit ? fit.fitEpochOrbit : out.orbit;
				if (fitCfg.computeUncertainty && fit.uncertainty.hasStateCovariance()) {
					out.stateCovariance = fit.uncertainty.cartesianStateCovariance;
					out.hasCovariance = !out.stateCovariance.isZero(0.0);
				}
				return out;
			} catch (const std::exception& e) {
				LOG_ERROR << "[ATT-EXT] OrbitOptimizer threw: " << e.what();
				return out;
			}
		}

		void runBestArcAutoRejectFO(const std::vector<ObsNode>& allObs,
		                            std::span<const orbital_core::obs::PreparedObservation> prepared,
		                            NBodyPropagator& nbody, const OrbitOptimizerConfig& baseCfg,
		                            uint16_t bestPerturberMask, Orbit_p& bestOrbit, std::vector<std::uint8_t>& bestMask,
		                            std::size_t arcFirst, std::size_t arcLast, double& bestRMS)
		{
			if (!bestOrbit || allObs.empty() || bestMask.size() != allObs.size() || arcFirst >= allObs.size() ||
			    arcLast >= allObs.size() || arcFirst > arcLast) {
				return;
			}

			for (std::size_t i = 0; i < arcFirst; ++i)
				bestMask[i] = 0;
			for (std::size_t i = arcLast + 1; i < allObs.size(); ++i)
				bestMask[i] = 0;

			const double arcDays = allObs[arcLast].JD - allObs[arcFirst].JD;
			const uint16_t mask = bestPerturberMask;

			OrbitOptimizerConfig fitCfg = baseCfg;
			fitCfg.nbody.enabled = (mask != 0);
			fitCfg.nbody.planetMask = mask;
			fitCfg.nonGrav.enabled = false;
			fitCfg.outputEpochTdb = 0.0;
			fitCfg.computeUncertainty = false;
			fitCfg.populateResiduals = false;
			fitCfg.logDiagnostics = false;
			fitCfg.logCacheStats = false;

			nbody.setPlanetMask(mask);

			LOG_INFO << "[ATT-EXT] auto-reject: start arc [" << arcFirst << ".." << arcLast << "] ("
			         << (arcLast - arcFirst + 1) << " obs over " << arcDays << "d), planetMask=0x" << std::hex << mask
			         << std::dec;

			ResidualEval res = evaluateResidualsArcsec(bestOrbit, prepared, nbody);
			for (int round = 0; round < kMaxAutoRejectRounds; ++round) {
				const auto maskBeforeReject = bestMask;
				const Orbit_p orbitBeforeReject = bestOrbit;

				bestRMS = rmsArcsecFromMask(res.totalArcsec, bestMask);

				AutoRejectWithinArcResult reject =
				    autoRejectWithinArcFO(bestMask, allObs, res.totalSigma, arcFirst, arcLast);

				if (reject.safeguardRollback) {
					bestMask = maskBeforeReject;
					LOG_INFO << "[ATT-EXT] auto-reject round " << round
					         << ": every sigma escalation violated safeguards"
					         << " (last sigma=" << reject.appliedSigma << ")";
					break;
				}

				if (!reject.changed) {
					LOG_INFO << "[ATT-EXT] auto-reject round " << round
					         << ": no mask changes at sigma=" << reject.appliedSigma << ", kept-RMS=" << bestRMS
					         << "\"";
					break;
				}

				auto kept = orbital_core::fit::selectObservations(prepared, bestMask);
				if (kept.size() < 3) {
					bestMask = maskBeforeReject;
					LOG_INFO << "[ATT-EXT] auto-reject round " << round
					         << ": kept set fell below 3 obs after safeguards"
					         << " -- restoring previous mask";
					break;
				}

				FitSubsetResult refined = fitOnSubset(bestOrbit, kept, fitCfg);
				if (!refined.orbit || !refined.converged) {
					bestMask = maskBeforeReject;
					bestOrbit = orbitBeforeReject;
					LOG_INFO << "[ATT-EXT] auto-reject round " << round
					         << ": refit did not converge -- restoring previous fitted mask";
					break;
				}

				bestOrbit = refined.orbit;
				res = evaluateResidualsArcsec(bestOrbit, prepared, nbody);
				bestRMS = rmsArcsecFromMask(res.totalArcsec, bestMask);

				LOG_INFO << "[ATT-EXT] auto-reject round " << round << ": sigma=" << reject.appliedSigma
				         << " in=" << reject.flippedIn << " out=" << reject.flippedOut
				         << " keptInArc=" << reject.keptInArc << "/" << (arcLast - arcFirst + 1)
				         << " totalKept=" << kept.size() << " refit-RMS=" << bestRMS << "\"";
			}
		}

	}  // namespace

	ExtensionResult attemptExtensions(const Orbit_p& seed, const std::vector<ObsNode>& allObs, std::size_t subarcStart,
	                                  std::size_t subarcEnd, const OrbitOptimizerConfig& baseCfg,
	                                  std::span<const orbital_core::obs::PreparedObservation> prepared)
	{
		ExtensionResult result;
		result.orbit = seed;

		if (!seed || allObs.empty() || subarcStart >= allObs.size() || subarcEnd >= allObs.size() ||
		    subarcStart > subarcEnd) {
			LOG_INFO << "[ATT-EXT] Invalid inputs — returning seed unchanged";
			result.keptObs = allObs;
			result.firstIncluded = 0;
			result.lastIncluded = allObs.empty() ? 0 : allObs.size() - 1;
			result.arcFirst = result.firstIncluded;
			result.arcLast = result.lastIncluded;
			return result;
		}

		for (std::size_t i = 1; i < allObs.size(); ++i) {
			if (allObs[i].JD < allObs[i - 1].JD) {
				LOG_ERROR << "[ATT-EXT] allObs not JD-sorted at index " << i << " — caller must sort upstream";
				result.keptObs = allObs;
				result.firstIncluded = 0;
				result.lastIncluded = allObs.size() - 1;
				result.arcFirst = result.firstIncluded;
				result.arcLast = result.lastIncluded;
				return result;
			}
		}

		const bool interstellarObject =
		    std::any_of(allObs.begin(), allObs.end(), [](const ObsNode& o) { return o.isInterstellar; });
		if (interstellarObject) {
			LOG_INFO << "[ATT-EXT] interstellar object flag enabled: "
			         << "hyperbolic eccentricity penalty disabled";
		}

		if (prepared.size() != allObs.size()) {
			LOG_ERROR << "[ATT-EXT] prepared observation count " << prepared.size() << " obs, expected "
			          << allObs.size() << " — abort";
			result.keptObs = allObs;
			result.firstIncluded = 0;
			result.lastIncluded = allObs.size() - 1;
			result.arcFirst = result.firstIncluded;
			result.arcLast = result.lastIncluded;
			return result;
		}

		NBodyPropagator nbody;
		nbody.setStepSize(baseCfg.nbody.stepSize);
		nbody.setUseSb441(baseCfg.nbody.useSb441);
		nbody.setTargetNaifId(baseCfg.nbody.targetNaifId);
		nbody.setInactiveInteriorMassPolicy(baseCfg.nbody.inactiveInteriorMassPolicy);
		if (baseCfg.nbody.enabled) {
			if (auto* ephemeris = getDefaultEphemerisSource())
				nbody = nbody.forEncounterMonitoring(*ephemeris);
		}
		const ArcScheduledPerturberSelector perturberSelector({baseCfg.nbody.enabled, baseCfg.nbody.planetMask,
		                                                       baseCfg.alwaysIncludedPerturberMask,
		                                                       baseCfg.perturberSelectionMode});

		std::vector<std::uint8_t> included(allObs.size(), 0);
		for (std::size_t i = subarcStart; i <= subarcEnd; ++i)
			included[i] = 1;

		Orbit_p currentOrbit = seed;
		Orbit_p bestOrbit = seed;
		orbital_core::StateCov6x6 currentCovariance{orbital_core::StateCov6x6::Zero()};
		bool currentHasCovariance{false};
		std::vector<std::uint8_t> bestMask = included;
		std::size_t bestFirstIdx = subarcStart;
		std::size_t bestLastIdx = subarcEnd;
		double bestRMS = 1.0e30;
		uint16_t carriedAutomaticMask = 0;
		const double initialArcDays = allObs[subarcEnd].JD - allObs[subarcStart].JD;
		double carriedArcDays = initialArcDays;
		PerturberSelection bestPerturberSelection = perturberSelector.select(carriedArcDays, carriedAutomaticMask);
		bool currentArcNeedsActiveForceFit = (baseCfg.nbody.enabled && bestPerturberSelection.effectiveMask != 0);
		std::optional<ResidualEval> nextResiduals;
		uint16_t nextResidualsMask = 0;
		LOG_INFO << "[ATT-EXT] Start: subarc [" << subarcStart << ".." << subarcEnd << "] ("
		         << (subarcEnd - subarcStart + 1) << " obs over " << initialArcDays << "d)";

		for (result.iterations = 0; result.iterations < kMaxOuterIter; ++result.iterations) {
			std::size_t firstIdx = 0;
			std::size_t lastIdx = 0;
			if (!firstLastIncluded(included, firstIdx, lastIdx))
				break;

			const double currentArcDays = allObs[lastIdx].JD - allObs[firstIdx].JD;
			PerturberSelection gateSelection = perturberSelector.select(currentArcDays, carriedAutomaticMask);
			nbody.setPlanetMask(gateSelection.effectiveMask);

			ResidualEval extensionRes =
			    nextResiduals && nextResidualsMask == gateSelection.effectiveMask
			        ? std::move(*nextResiduals)
			        : evaluateResidualsArcsec(currentOrbit, prepared, nbody, &perturberSelector);
			nextResiduals.reset();
			const uint16_t detectedAutoMask =
			    static_cast<uint16_t>(extensionRes.autoPerturberMask & ~carriedAutomaticMask);
			bool refitCurrentArcForPerturbers = false;
			if (detectedAutoMask != 0) {
				const uint16_t previousEffectiveMask = gateSelection.effectiveMask;
				carriedAutomaticMask = static_cast<uint16_t>(carriedAutomaticMask | detectedAutoMask);
				gateSelection = perturberSelector.select(currentArcDays, carriedAutomaticMask);
				nbody.setPlanetMask(gateSelection.effectiveMask);
				refitCurrentArcForPerturbers = (gateSelection.effectiveMask != previousEffectiveMask);
				if (baseCfg.logDiagnostics) {
					LOG_INFO << "[ATT-EXT] iter " << result.iterations << " automatic perturbers added mask=0x"
					         << std::hex << detectedAutoMask << " effective=0x" << gateSelection.effectiveMask
					         << " blocked=0x" << gateSelection.blockedMask << std::dec;
				}
				extensionRes = evaluateResidualsArcsec(currentOrbit, prepared, nbody, &perturberSelector);
			}
			const int nAdded = extendInclusionMaskFO(included, allObs, extensionRes.totalSigma);
			const bool refitCurrentArcForActiveForce =
			    (currentArcNeedsActiveForceFit && gateSelection.effectiveMask != 0);
			if (!nAdded && !refitCurrentArcForPerturbers && !refitCurrentArcForActiveForce) {
				LOG_INFO << "[ATT-EXT] iter " << result.iterations << ": no obs within extension gates (edge < "
				         << kResidualLimitSigma << " sigma, window < " << kWindowMaxSigma
				         << " sigma) -- extension stalled, converged";
				result.converged = true;
				break;
			}
			if (!nAdded) {
				if (refitCurrentArcForPerturbers) {
					LOG_INFO << "[ATT-EXT] iter " << result.iterations
					         << ": no obs added, but perturber mask changed to 0x" << std::hex
					         << gateSelection.effectiveMask << std::dec << "; refitting current arc";
				} else {
					LOG_INFO << "[ATT-EXT] iter " << result.iterations << ": no obs added before active-force fit; "
					         << "refitting current arc";
				}
			}

			if (!firstLastIncluded(included, firstIdx, lastIdx))
				break;

			const double arcDays = allObs[lastIdx].JD - allObs[firstIdx].JD;
			PerturberSelection fitSelection = perturberSelector.select(arcDays, carriedAutomaticMask);
			uint16_t mask = fitSelection.effectiveMask;
			OrbitOptimizerConfig fitCfg = baseCfg;
			fitCfg.nbody.enabled = (mask != 0);
			fitCfg.nbody.planetMask = mask;
			fitCfg.nonGrav.enabled = false;
			fitCfg.outputEpochTdb = 0.0;
			fitCfg.computeUncertainty = true;
			fitCfg.populateResiduals = false;
			fitCfg.logDiagnostics = false;
			fitCfg.logCacheStats = false;
			const auto kept = orbital_core::fit::selectObservations(prepared, included);
			LOG_INFO << "[ATT-EXT] iter " << result.iterations << ": extension added " << nAdded << " obs; trying=["
			         << firstIdx << ".." << lastIdx << "] (" << kept.size() << " obs, " << arcDays
			         << "d), planetMask=0x" << std::hex << mask << std::dec << " blocked=0x" << std::hex
			         << fitSelection.blockedMask << std::dec << ", useNBody=" << fitCfg.nbody.enabled;

			Orbit_p latestOrbit;
			ResidualEval latestRes;
			double latestScore = iod_scoring::kBadScore;
			double latestWeightedRMS = iod_scoring::kBadScore;
			double latestSubRMS = 1.0e30;
			bool fitFailed = false;
			bool catastrophic = false;
			bool convergedThisArc = false;

			for (int pass = 0; pass < kMaxForceModelFits; ++pass) {
				nbody.setPlanetMask(mask);

				if (currentHasCovariance && currentOrbit) {
					LovVector6 lov{};
					if (oneSigmaLovFromCovarianceFO(currentCovariance, lov)) {
						ResidualEval lovRes = pass == 0 && mask == gateSelection.effectiveMask
						                          ? std::move(extensionRes)
						                          : evaluateResidualsArcsec(currentOrbit, prepared, nbody);
						double lovRms = iod_scoring::weightedRmsFromTotalSigma(lovRes.totalSigma, included);
						double prevLovRms = 1.0e70;
						for (int lovIter = 0;
						     lovIter < 5 && lovRms > 2.0 && (lovRms > 2000.0 || prevLovRms > lovRms + 0.1); ++lovIter) {
							prevLovRms = lovRms;
							const LovImproveResult lovStep =
							    improveAlongLovFO(currentOrbit, lov, allObs, prepared, included, nbody);
							if (!lovStep.applied || !lovStep.orbit)
								break;

							currentOrbit = lovStep.orbit;
							lovRes = evaluateResidualsArcsec(currentOrbit, prepared, nbody);
							lovRms = iod_scoring::weightedRmsFromTotalSigma(lovRes.totalSigma, included);
							LOG_INFO << "[ATT-EXT] iter " << result.iterations << " pass " << pass << " LOV " << lovIter
							         << ": minimum at " << lovStep.sigmas << " sigmas; weightedRMS=" << lovRms;
						}
					}
				}

				FitSubsetResult refined = fitOnSubset(currentOrbit, kept, fitCfg);
				if (!refined.orbit || !refined.converged) {
					LOG_INFO << "[ATT-EXT] iter " << result.iterations << " pass " << pass
					         << ": fit did not converge -- aborting extension";
					fitFailed = true;
					break;
				}

				currentOrbit = refined.orbit;
				currentCovariance = refined.stateCovariance;
				currentHasCovariance = refined.hasCovariance;
				latestOrbit = refined.orbit;

				nbody.setPlanetMask(mask);
				latestRes = evaluateResidualsArcsec(refined.orbit, prepared, nbody, &perturberSelector);
				const uint16_t fittedMask = mask;
				const uint16_t passAutoMask =
				    static_cast<uint16_t>(latestRes.autoPerturberMask & ~carriedAutomaticMask);
				if (passAutoMask != 0) {
					carriedAutomaticMask = static_cast<uint16_t>(carriedAutomaticMask | passAutoMask);
					fitSelection = perturberSelector.select(arcDays, carriedAutomaticMask);
					mask = fitSelection.effectiveMask;
					fitCfg.nbody.planetMask = mask;
					fitCfg.nbody.enabled = (mask != 0);
					nbody.setPlanetMask(mask);
					if (baseCfg.logDiagnostics) {
						LOG_INFO << "[ATT-EXT] iter " << result.iterations << " pass " << pass
						         << " automatic perturbers added mask=0x" << std::hex << passAutoMask << " fitMask=0x"
						         << mask << " blocked=0x" << fitSelection.blockedMask << std::dec;
					}
					latestRes = evaluateResidualsArcsec(refined.orbit, prepared, nbody, &perturberSelector);
				}

				double subSumSq = 0.0;
				std::size_t subCount = 0;
				for (std::size_t i = 0; i < allObs.size(); ++i) {
					if (!included[i])
						continue;
					if (std::isfinite(latestRes.totalArcsec[i]) && latestRes.totalArcsec[i] < 1.0e29) {
						subSumSq += latestRes.totalArcsec[i] * latestRes.totalArcsec[i];
						++subCount;
					}
				}
				latestSubRMS = (subCount > 0) ? std::sqrt(subSumSq / static_cast<double>(subCount)) : 1.0e30;

				const std::size_t includedCount = countIncluded(included);
				latestWeightedRMS = iod_scoring::weightedRmsFromTotalSigma(latestRes.totalSigma, included);
				latestScore =
				    iod_scoring::scoreOrbit(refined.orbit, latestWeightedRMS, includedCount, interstellarObject);

				convergedThisArc = (mask == fittedMask);

				if (latestScore > iod_scoring::kFailureScore) {
					LOG_INFO << "[ATT-EXT] iter " << result.iterations << " pass " << pass << ": score=" << latestScore
					         << " > " << iod_scoring::kFailureScore << " — LM catastrophically diverged, reverting";
					catastrophic = true;
					break;
				}

				LOG_INFO << "[ATT-EXT] iter " << result.iterations << " pass " << pass << ": score=" << latestScore
				         << " weightedRMS=" << latestWeightedRMS << " kept-RMS=" << latestSubRMS
				         << "\" full-RMS=" << latestRes.rmsAll << "\" converged=" << convergedThisArc
				         << " refitRequired=" << !convergedThisArc;

				if (convergedThisArc)
					break;
			}

			if (fitFailed || catastrophic)
				break;

			if (!latestOrbit) {
				LOG_INFO << "[ATT-EXT] iter " << result.iterations << ": no usable fit result — aborting extension";
				break;
			}

			currentArcNeedsActiveForceFit = false;

			if (convergedThisArc) {
				if (latestScore > iod_scoring::kAcceptableExtensionScore) {
					LOG_INFO << "[ATT-EXT] iter " << result.iterations << ": converged score=" << latestScore << " > "
					         << iod_scoring::kAcceptableExtensionScore << " (weightedRMS=" << latestWeightedRMS
					         << ", kept-RMS=" << latestSubRMS << "\") -- not promoting this extension";
					break;
				}

				bestOrbit = latestOrbit;
				bestMask = included;
				bestRMS = latestSubRMS;
				bestFirstIdx = firstIdx;
				bestLastIdx = lastIdx;
				bestPerturberSelection = fitSelection;
				carriedArcDays = arcDays;
				nextResiduals = std::move(latestRes);
				nextResidualsMask = mask;
			} else {
				LOG_INFO << "[ATT-EXT] iter " << result.iterations << ": force model did not stabilize after "
				         << kMaxForceModelFits << " fits; restoring previous fitted arc";
				break;
			}
		}

		runBestArcAutoRejectFO(allObs, prepared, nbody, baseCfg, bestPerturberSelection.effectiveMask, bestOrbit,
		                       bestMask, bestFirstIdx, bestLastIdx, bestRMS);

		std::size_t finalFirstIdx = bestFirstIdx;
		std::size_t finalLastIdx = bestLastIdx;
		if (!firstLastIncluded(bestMask, finalFirstIdx, finalLastIdx)) {
			finalFirstIdx = bestFirstIdx;
			finalLastIdx = bestLastIdx;
		}

		result.orbit = bestOrbit;
		result.firstIncluded = finalFirstIdx;
		result.lastIncluded = finalLastIdx;
		result.arcFirst = bestFirstIdx;
		result.arcLast = bestLastIdx;
		result.finalRMS = bestRMS;
		result.includedMask = bestMask;
		result.keptObs = buildKeptFromMask(allObs, bestMask);
		result.keptObsIndices = buildKeptIndices(bestMask);
		result.configuredPerturberMask = perturberSelector.userMask();
		result.effectivePerturberMask = bestPerturberSelection.effectiveMask;
		result.automaticPerturberMask = bestPerturberSelection.automaticMask;
		result.blockedPerturberMask = bestPerturberSelection.blockedMask;

		LOG_INFO << "[ATT-EXT] Done: [" << finalFirstIdx << ".." << finalLastIdx << "] (" << result.keptObs.size()
		         << " obs over " << (allObs[finalLastIdx].JD - allObs[finalFirstIdx].JD) << "d), "
		         << "planetMask=0x" << std::hex << result.effectivePerturberMask << " auto=0x"
		         << result.automaticPerturberMask << " blocked=0x" << result.blockedPerturberMask << std::dec << ", "
		         << "iters=" << result.iterations << ", final kept-RMS=" << result.finalRMS << "\""
		         << ", converged=" << result.converged;

		return result;
	}

}  // namespace Calculations
