
#include "InitialOrbit.h"
#include <utility>
#include <vector>


#include "Gauss.h"
#include "Herget.h"
#include "TwoBodyPropagation.h"
#include "IodScoring.h"
#include "Vaisala.h"
#include "logger/Logger.h"
#include "obs/AstrometricSigma.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/Types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <span>

namespace Calculations {
	namespace {
		constexpr double kBadCost = 1.0e30;
		constexpr double kBadScore = 1.0e30;

		constexpr double kAcceptableScore = 5.0;

		constexpr double kMaxInitialArclenDays = 730.0;

		constexpr double kMaxArclenForScansDays = 230.0;

		constexpr double kVaisalaFlippedOutScore = 5.0e4;

		constexpr double kArclenShrink = 0.7;

		constexpr double kMinArclenDays = 0.001;

		constexpr int kImprovementMaxIter = 5;
		constexpr double kImprovementTriggerScore = 50.0;
		constexpr double kImprovementAcceptSlack = 0.01;
		constexpr double kFullStepAdvantage = 0.3;
		constexpr double kImprovementFailureScore = 10000.0;

		constexpr double kVaisalaScanRStart = 0.1;
		constexpr double kVaisalaScanRMax = 200.0;

		constexpr double kHergetScanRMax = 5.0;
		constexpr double kHergetScanRBase = 0.004;
		constexpr double kHergetScanRExpo = 0.6666;

		constexpr double kScanStepFactor = 1.2;

		[[nodiscard]] double evaluateInitialOrbit(const Vec3D& R, const Vec3D& V, double epochJD,
		                                          const std::vector<ObsNode>& obs, bool interstellar,
		                                          orbital_core::AstrometricWeightProfile weightProfile) noexcept
		{
			return iod_scoring::scoreStateKepler(R, V, epochJD, obs, interstellar, weightProfile);
		}

		[[nodiscard]] bool hasInterstellarObject(const std::vector<ObsNode>& obs) noexcept
		{
			return std::any_of(obs.begin(), obs.end(), [](const ObsNode& o) { return o.isInterstellar; });
		}

		[[nodiscard]] double scoreOrbitArc(const ObsNode* obs, std::size_t n_obs) noexcept
		{
			if (n_obs <= 1)
				return -1.0;

			constexpr double kFiveDegrees = 5.0 * Constants::PI / 180.0;
			const Vec3D xprod = obs[0].vector.cross(obs[n_obs - 1].vector);

			double minimum = 0.0;
			double maximum = 0.0;
			double max_dist = (obs[0].vector - obs[n_obs - 1].vector).magnitude();

			for (std::size_t i = 1; i + 1 < n_obs; ++i) {
				const double dot = xprod.dot(obs[i].vector);
				if (dot > maximum)
					maximum = dot;
				else if (dot < minimum)
					minimum = dot;

				if (max_dist < kFiveDegrees) {
					const double d1 = (obs[i].vector - obs[0].vector).magnitude();
					const double d2 = (obs[i].vector - obs[n_obs - 1].vector).magnitude();
					if (d1 > max_dist)
						max_dist = d1;
					if (d2 > max_dist)
						max_dist = d2;
				}
			}
			if (max_dist > kFiveDegrees)
				max_dist = kFiveDegrees;

			return (maximum - minimum) * max_dist;
		}

		void lookForBestSubarc(const std::vector<ObsNode>& obs, double max_arc_len, std::size_t& start,
		                       std::size_t& end) noexcept
		{
			const double cos_45_deg = 1.414213 / 2.0;
			double best_score = -999.0;
			start = 0;
			end = 0;

			const std::size_t n = obs.size();
			if (n < 3)
				return;

			std::size_t scanStart = 0;

			std::size_t j = scanStart;
			for (std::size_t i = scanStart; i + 2 < n; ++i) {
				if (j < i)
					j = i;

				while (j + 1 < n && (obs[j + 1].epochJdTdb() - obs[i].epochJdTdb()) <= max_arc_len &&
				       obs[i].vector.dot(obs[j + 1].vector) > cos_45_deg) {
					++j;
				}

				std::size_t ts = i;
				std::size_t te = j;
				if (obs[ts].epochJdTdb() == obs[te].epochJdTdb())
					te = ts;
				if (te <= ts + 1)
					continue;

				const double s = scoreOrbitArc(obs.data() + ts, te - ts + 1);
				if (s > best_score) {
					best_score = s;
					start = ts;
					end = te;
				}
			}
		}

		[[nodiscard]] std::size_t pickGaussMiddleIdx(const std::vector<ObsNode>& subarc) noexcept
		{
			const std::size_t n = subarc.size();
			if (n <= 3)
				return n / 2;

			const Vec3D normal = subarc.front().vector.cross(subarc.back().vector);

			std::size_t bestIdx = n / 2;
			double bestAbs = -1.0;
			for (std::size_t k = 1; k + 1 < n; ++k) {
				const double d = std::abs(normal.dot(subarc[k].vector));
				if (d > bestAbs) {
					bestAbs = d;
					bestIdx = k;
				}
			}
			return bestIdx;
		}

		struct IodImprovementSeed {
			StateVectors state;
			double epoch{0.0};
			double rawScore{kBadScore};
			double selectionScore{kBadScore};
			bool changed{false};
		};

		[[nodiscard]] double scoreIodState(const StateVectors& sv, double epoch, const std::vector<ObsNode>& obs,
		                                   orbital_core::AstrometricWeightProfile weightProfile) noexcept
		{
			if (!sv.isValid())
				return kBadScore;
			const bool interstellar = hasInterstellarObject(obs);
			const double score = evaluateInitialOrbit(sv.first, sv.second, epoch, obs, interstellar, weightProfile);
			return std::isfinite(score) ? score : kBadScore;
		}

		[[nodiscard]] double topocentricRhoFromState(const StateVectors& sv, double epoch, const ObsNode& obs) noexcept
		{
			const StateVectors atObs = propagateKeplerState(sv, obs.epochJdTdb() - epoch);
			if (!atObs.isValid())
				return -1.0;

			const double rho = (atObs.first + obs.obsVector).magnitude();
			return std::isfinite(rho) ? rho : -1.0;
		}

        class IodResidualModel final {
		   public:
            IodResidualModel(std::span<const ObsNode> obs,
			                   orbital_core::AstrometricWeightProfile weightProfile) noexcept
			    : m_obs(obs), m_weightProfile(weightProfile)
			{
			}

			[[nodiscard]] bool compute(const StateVectors& sv, double epoch, Eigen::VectorXd& out) const
			{
				if (!sv.isValid() || m_obs.empty())
					return false;

				out.resize(static_cast<Eigen::Index>(2 * m_obs.size()));
				for (std::size_t i = 0; i < m_obs.size(); ++i) {
					const auto [raModel, decModel] = modelRaDec(sv, epoch, m_obs[i]);
					if (!std::isfinite(raModel) || !std::isfinite(decModel) || std::abs(raModel) > Constants::TWOPI ||
					    std::abs(decModel) > Constants::PI) {
						return false;
					}

					const auto [sigmaRa, sigmaDec] = sigmasFor(m_obs[i]);
					const double dRA = Helpers::wrapPi(m_obs[i].RA - raModel);
					const double dDec = m_obs[i].DEC - decModel;

					out[static_cast<Eigen::Index>(2 * i)] = std::cos(m_obs[i].DEC) * dRA / sigmaRa;
					out[static_cast<Eigen::Index>(2 * i + 1)] = dDec / sigmaDec;

					if (!std::isfinite(out[static_cast<Eigen::Index>(2 * i)]) ||
					    !std::isfinite(out[static_cast<Eigen::Index>(2 * i + 1)])) {
						return false;
					}
				}
				return true;
			}

		   private:
			std::span<const ObsNode> m_obs;
			orbital_core::AstrometricWeightProfile m_weightProfile;

			[[nodiscard]] static std::pair<double, double> modelRaDec(const StateVectors& sv, double epoch,
			                                                          const ObsNode& obs)
			{
				const StateVectors atObs = propagateKeplerState(sv, obs.epochJdTdb() - epoch);
				if (!atObs.isValid())
					return {kBadScore, kBadScore};

				const Vec3D losEcl = atObs.first + obs.obsVector;
				const Vec3D losEq(
				    losEcl.x, losEcl.y * std::cos(Constants::EARTH_INCL) - losEcl.z * std::sin(Constants::EARTH_INCL),
				    losEcl.y * std::sin(Constants::EARTH_INCL) + losEcl.z * std::cos(Constants::EARTH_INCL));
				const double rho = losEq.magnitude();
				if (rho <= 0.0 || !std::isfinite(rho))
					return {kBadScore, kBadScore};

				double ra = std::atan2(losEq.y, losEq.x);
				if (ra < 0.0)
					ra += Constants::TWOPI;
				const double dec = std::asin(std::clamp(losEq.z / rho, -1.0, 1.0));
				return {ra, dec};
			}

			[[nodiscard]] std::pair<double, double> sigmasFor(const ObsNode& obs) const
			{
				const orbital_core::ObservationUncertainty uncertainty =
				    orbital_core::obs::lookupUncertainty(obs, obs.JD, m_weightProfile);
				double sigmaRa = uncertainty.sigmaRaRad;
				double sigmaDec = uncertainty.sigmaDecRad;
				return {sigmaRa > 0.0 ? sigmaRa : Constants::ARCSEC_RAD,
				        sigmaDec > 0.0 ? sigmaDec : Constants::ARCSEC_RAD};
			}
		};

        class FixedDistanceEndpointFitter final {
		   public:
            explicit FixedDistanceEndpointFitter(orbital_core::AstrometricWeightProfile weightProfile)
			    : m_weightProfile(weightProfile)
			{
				orbital_core::obs::ensureSigmaTable(m_weightProfile);
			}

			[[nodiscard]] IodImprovementSeed adjust(const StateVectors& seed, double epoch,
			                                        const std::vector<ObsNode>& obs) const
			{
				IodImprovementSeed out;
				if (!seed.isValid() || obs.size() < 3)
					return out;

				const double rho1 = topocentricRhoFromState(seed, epoch, obs.front());
				const double rhoN = topocentricRhoFromState(seed, epoch, obs.back());
				if (rho1 <= 0.0 || rhoN <= 0.0)
					return out;

                const IodResidualModel residualModel(obs, m_weightProfile);
				Eigen::VectorXd base;
				if (!residualModel.compute(seed, epoch, base))
					return out;

				Eigen::MatrixXd slopes(base.size(), kFitParams);
				for (int p = 0; p < kFitParams; ++p) {
					Params params{};
					params[p] = -kDeltaRadians;

					const StateVectors tweaked = transferWithEndpointOffsets(obs, rho1, rhoN, params);
					Eigen::VectorXd shifted;
					if (!residualModel.compute(tweaked, obs.front().epochJdTdb(), shifted) ||
					    shifted.size() != base.size()) {
						return out;
					}

					slopes.col(p) = (shifted - base) / kDeltaRadians;
				}

				const Eigen::VectorXd rhs = slopes.transpose() * base;
				const Eigen::MatrixXd normal = slopes.transpose() * slopes;
				const Eigen::VectorXd delta = normal.completeOrthogonalDecomposition().solve(rhs);
				if (delta.size() != kFitParams || !delta.allFinite())
					return out;

				Params params{};
				for (int i = 0; i < kFitParams; ++i)
					params[i] = delta[i];

				const StateVectors adjusted = transferWithEndpointOffsets(obs, rho1, rhoN, params);
				const double score = scoreIodState(adjusted, obs.front().epochJdTdb(), obs, m_weightProfile);
				if (!adjusted.isValid() || !std::isfinite(score) || score >= kBadScore)
					return out;

				out.state = adjusted;
				out.epoch = obs.front().epochJdTdb();
				out.rawScore = score;
				out.selectionScore = score;
				out.changed = true;
				return out;
			}

		   private:
			orbital_core::AstrometricWeightProfile m_weightProfile;
			static constexpr int kFitParams = 4;
			static constexpr double kDeltaRadians = 1.0e-6;
			using Params = std::array<double, kFitParams>;

			[[nodiscard]] static ObsNode offsetEndpoint(const ObsNode& obs, double dRa, double dDec)
			{
				ObsNode out = obs;
				out.RA += dRa;
				out.DEC = std::clamp(out.DEC + dDec, -Constants::PI / 2.0, Constants::PI / 2.0);

				const Vec3D equatorial(std::cos(out.DEC) * std::cos(out.RA), std::cos(out.DEC) * std::sin(out.RA),
				                       std::sin(out.DEC));
				out.vector = Convertors::ecvToEcl(equatorial, Constants::EARTH_INCL);
				return out;
			}

			[[nodiscard]] static StateVectors transferWithEndpointOffsets(const std::vector<ObsNode>& obs, double rho1,
			                                                              double rhoN, const Params& params)
			{
				if (obs.size() < 2)
					return {};

				const ObsNode first = offsetEndpoint(obs.front(), params[0], params[1]);
				const ObsNode last = offsetEndpoint(obs.back(), params[2], params[3]);
				return transferOrbitFromEndpointDistances(first, last, rho1, rhoN);
			}
		};

		class FoFullImprovementFitter final {
		   public:
			explicit FoFullImprovementFitter(orbital_core::AstrometricWeightProfile weightProfile) noexcept
			    : m_weightProfile(weightProfile)
			{
			}

			[[nodiscard]] IodImprovementSeed improve(const IodImprovementSeed& seed,
			                                         const std::vector<ObsNode>& subarc) const
			{
				IodImprovementSeed out;
				if (!seed.state.isValid() || subarc.size() < 3)
					return out;

                const IodResidualModel residualModel(subarc, m_weightProfile);
				Eigen::VectorXd base;
				if (!residualModel.compute(seed.state, seed.epoch, base))
					return out;

				Eigen::MatrixXd jac(base.size(), kParams);
				if (!computeJacobian(seed, subarc, residualModel, base, jac))
					return out;

				const Eigen::Matrix<double, kParams, kParams> normal = jac.transpose() * jac;
				const Eigen::VectorXd rhs = jac.transpose() * base;
				const Eigen::VectorXd delta = normal.completeOrthogonalDecomposition().solve(rhs);
				if (!delta.allFinite())
					return out;

				const double scale = stepScale(seed, subarc, delta);
				Params params = paramsFromState(seed.state);
				params -= delta / scale;

				const StateVectors candidate = stateFromParams(params);
				const double score = scoreIodState(candidate, seed.epoch, subarc, m_weightProfile);
				if (!candidate.isValid() || !std::isfinite(score) || score >= kBadScore)
					return out;

				out.state = candidate;
				out.epoch = seed.epoch;
				out.rawScore = score;
				out.selectionScore = score - kFullStepAdvantage;
				out.changed = true;
				return out;
			}

		   private:
			orbital_core::AstrometricWeightProfile m_weightProfile;
			static constexpr int kParams = 6;
			static constexpr int kMaxDeltaIterations = 100;
			static constexpr int kMaxDeltaTweaks = 30;
			static constexpr double kMinChangeSigmas = 0.03;
			static constexpr double kMaxChangeSigmas = 3.0;
			static constexpr double kOptimalChangeSigmas = 1.0;
			static constexpr double kShortArcDays = 7.0;
			static constexpr double kStateStepLimitFraction = 0.7;
			inline static constexpr std::array<double, kParams> kDefaultDeltaVals{1.0e-4, 1.0e-4, 1.0e-5,
			                                                                      1.0e-5, 1.0e-3, 1.0e-3};

			using Params = Eigen::Matrix<double, kParams, 1>;
			mutable std::array<double, kParams> m_deltaVals{kDefaultDeltaVals};

			[[nodiscard]] static Params paramsFromState(const StateVectors& state) noexcept
			{
				Params p;
				p << state.first.x, state.first.y, state.first.z, state.second.x, state.second.y, state.second.z;
				return p;
			}

			[[nodiscard]] static StateVectors stateFromParams(const Params& p) noexcept
			{
				return StateVectors(Vec3D(p[0], p[1], p[2]), Vec3D(p[3], p[4], p[5]));
			}

			[[nodiscard]] static bool residualsForParams(const Params& params, double epoch,
                                                         const IodResidualModel& residualModel, Eigen::VectorXd& out)
			{
				return residualModel.compute(stateFromParams(params), epoch, out);
			}

			[[nodiscard]] bool computeJacobian(const IodImprovementSeed& seed, const std::vector<ObsNode>& subarc,
                                               const IodResidualModel& residualModel, const Eigen::VectorXd& base,
			                                   Eigen::MatrixXd& jac) const
			{
				const Params params = paramsFromState(seed.state);
				const double integrationLength =
				    std::max(1.0e-6, std::max(std::abs(seed.epoch - subarc.front().epochJdTdb()),
				                              std::abs(seed.epoch - subarc.back().epochJdTdb())));
				const bool preferSymmetric = (subarc.back().epochJdTdb() - subarc.front().epochJdTdb()) < kShortArcDays;

				for (int col = 0; col < kParams; ++col) {
					bool acceptedColumn = false;
					double lowDelta = 0.0;
					double highDelta = 0.0;
					double lowChange = 0.0;
					double highChange = 0.0;

					for (int iter = 0; iter < kMaxDeltaIterations && !acceptedColumn; ++iter) {
						double deltaVal = m_deltaVals[col] / (integrationLength * integrationLength);
						if (deltaVal <= 0.0 || !std::isfinite(deltaVal))
							return false;

						Eigen::VectorXd minus;
						bool minusOk = false;
						for (int tweak = 0; tweak < kMaxDeltaTweaks; ++tweak) {
							Params tweaked = params;
							tweaked[col] -= deltaVal;
							if (residualsForParams(tweaked, seed.epoch, residualModel, minus)) {
								minusOk = true;
								break;
							}

							deltaVal *= 0.5;
							m_deltaVals[col] *= 0.5;
						}
						if (!minusOk)
							return false;

						bool useSymmetric = preferSymmetric;
						Eigen::VectorXd plus;
						if (useSymmetric) {
							Params tweaked = params;
							tweaked[col] += deltaVal;
							useSymmetric = residualsForParams(tweaked, seed.epoch, residualModel, plus);
						}

						const Eigen::VectorXd change = useSymmetric ? (plus - minus) : (minus - base);
						double worstChangeSquared = 0.0;
						for (Eigen::Index row = 0; row + 1 < change.size(); row += 2) {
							const double change2 = change[row] * change[row] + change[row + 1] * change[row + 1];
							worstChangeSquared = std::max(worstChangeSquared, change2);
						}

						double worstChange = std::sqrt(worstChangeSquared) + 1.0e-10;
						if (worstChange > kMinChangeSigmas && worstChange < kMaxChangeSigmas) {
							if (useSymmetric)
								jac.col(col) = (plus - minus) / (2.0 * deltaVal);
							else
								jac.col(col) = (base - minus) / deltaVal;
							acceptedColumn = true;
							break;
						}

						if (worstChange <= kOptimalChangeSigmas) {
							lowDelta = m_deltaVals[col];
							lowChange = worstChange;
						} else {
							highDelta = m_deltaVals[col];
							highChange = worstChange;
						}

						double rescale = kOptimalChangeSigmas / worstChange;
						if (lowChange > 0.0 && highChange > 0.0 && (iter & 1)) {
							const double slope = std::log(highDelta / lowDelta) / std::log(highChange / lowChange);
							rescale = std::exp(std::log(rescale) * slope);
						}
						m_deltaVals[col] *= rescale;
					}

					if (!acceptedColumn)
						return false;
				}

				return jac.allFinite();
			}

			[[nodiscard]] static double stepScale(const IodImprovementSeed& seed, const std::vector<ObsNode>& subarc,
			                                      const Params& delta) noexcept
			{
				double scale = 1.0;
				double rho = topocentricRhoFromState(seed.state, seed.epoch, subarc.front());
				if (rho <= 0.0 || !std::isfinite(rho))
					rho = seed.state.first.magnitude();

				const double span = std::max(1.0e-9, subarc.back().epochJdTdb() - subarc.front().epochJdTdb());
				for (int i = 0; i < kParams; ++i) {
					double maxDifference = kStateStepLimitFraction * rho;
					if (i >= 3)
						maxDifference = (maxDifference / (span * 0.5)) / Constants::GAUSS_K;
					if (maxDifference > 0.0 && std::isfinite(maxDifference))
						scale = std::max(scale, std::abs(delta[i]) / maxDifference);
				}

				return scale;
			}
		};

		[[nodiscard]] IodImprovementSeed attemptIodImprovements(
		    const StateVectors& seedState, double seedEpoch, const std::vector<ObsNode>& subarc, double seedScore,
		    orbital_core::AstrometricWeightProfile weightProfile,
            const FixedDistanceEndpointFitter& fixedDistanceFitter,
		    const FoFullImprovementFitter& fullImprovementFitter)
		{
			if (!seedState.isValid() || subarc.size() < 3 || !std::isfinite(seedScore) || seedScore >= kBadScore) {
				return IodImprovementSeed{seedState, seedEpoch, seedScore, seedScore, false};
			}

			StateVectors stateAtFirst = seedState;
			double epochAtFirst = seedEpoch;
			if (std::abs(seedEpoch - subarc.front().epochJdTdb()) > 1.0e-12) {
				const StateVectors propagated =
				    propagateKeplerState(seedState, subarc.front().epochJdTdb() - seedEpoch);
				if (propagated.isValid()) {
					stateAtFirst = propagated;
					epochAtFirst = subarc.front().epochJdTdb();
				}
			}

			const double currentScore = scoreIodState(stateAtFirst, epochAtFirst, subarc, weightProfile);
			IodImprovementSeed current{stateAtFirst, epochAtFirst, currentScore, currentScore, false};
			if (!std::isfinite(currentScore) || currentScore >= kBadScore)
				return current;

			IodImprovementSeed best = current;

			for (int method = 0; method < 2; ++method) {
				bool errorOccurred = false;

				for (int iter = 0; !errorOccurred && iter < kImprovementMaxIter; ++iter) {
					IodImprovementSeed trial;

					if (method == 0) {
						const double rho1 = topocentricRhoFromState(current.state, current.epoch, subarc.front());
						const double rhoN = topocentricRhoFromState(current.state, current.epoch, subarc.back());
						if (rho1 <= 0.0 || rhoN <= 0.0) {
							errorOccurred = true;
							break;
						}

						trial.epoch = subarc.front().epochJdTdb();
						const StateVectors rawHerget = herget(subarc, rho1, rhoN);
						trial = fixedDistanceFitter.adjust(rawHerget, trial.epoch, subarc);
						if (!trial.changed) {
							errorOccurred = true;
							break;
						}
					} else {
						trial = fullImprovementFitter.improve(current, subarc);
					}

					if (!trial.state.isValid() || !std::isfinite(trial.rawScore) || trial.rawScore >= kBadScore) {
						errorOccurred = true;
						break;
					}

					if (trial.rawScore > kImprovementFailureScore)
						errorOccurred = true;

					if (trial.selectionScore < current.selectionScore + kImprovementAcceptSlack) {
						current = trial;
						best = trial;
					}
				}
			}

			return best;
		}

	}  // namespace

	InitialOrbitResult bestInitialOrbit(const std::vector<ObsNode>& groundObs, const std::vector<ObsNode>& allObs,
	                                    orbital_core::AstrometricWeightProfile weightProfile)
	{
		InitialOrbitResult best;
		best.costRad2 = kBadCost;
		best.score = kBadScore;

		if (groundObs.size() < 3 || allObs.empty())
			return best;

		const bool interstellarObject = hasInterstellarObject(allObs) || hasInterstellarObject(groundObs);
		if (interstellarObject) {
			LOG_INFO << "[IOD] interstellar object flag enabled: "
			         << "hyperbolic eccentricity penalty disabled";
		}

        const FixedDistanceEndpointFitter fixedDistanceFitter(weightProfile);
		const FoFullImprovementFitter fullImprovementFitter(weightProfile);

		const double totalSpanDays = groundObs.back().epochJdTdb() - groundObs.front().epochJdTdb();
		double arclen = std::min(totalSpanDays, kMaxInitialArclenDays);

		bool hasGeocentricObs = false;
		for (const auto& o : groundObs) {
			if (o.obsVector.magnitude() < 2.0 * Constants::EARTH_RADIUS_IN_AU) {
				hasGeocentricObs = true;
				break;
			}
		}

		const std::vector<ObsNode>* subarcRef = nullptr;

		std::size_t winningSubarcStart = static_cast<std::size_t>(-1);
		std::size_t winningSubarcEnd = static_cast<std::size_t>(-1);

		auto tryCandidateWithSelection = [&](const StateVectors& sv, double epoch, const char* name, bool verbose,
		                                     double selectionScoreOverride) -> double {
			if (!sv.isValid())
				return kBadScore;

			const auto& scoringObs = (subarcRef != nullptr) ? *subarcRef : allObs;
			const double score =
			    evaluateInitialOrbit(sv.first, sv.second, epoch, scoringObs, interstellarObject, weightProfile);
			if (!std::isfinite(score) || score >= kBadScore)
				return kBadScore;

			const double costAllObs = keplerResidualCost(sv.first, sv.second, epoch, allObs);

			if (verbose) {
				const double cost = keplerResidualCost(sv.first, sv.second, epoch, scoringObs);
				const double rms = std::sqrt(cost / static_cast<double>(scoringObs.size())) / Constants::ARCSEC_RAD;
				LOG_INFO << "[IOD] " << name << " RMS=" << rms << "\"  score=" << score;
			}

			const double selectionScore = (std::isfinite(selectionScoreOverride) && selectionScoreOverride < kBadScore)
			                                  ? selectionScoreOverride
			                                  : score;

			if (selectionScore < best.score) {
				best.state = sv;
				best.epoch = epoch;
				best.method = name;
				best.score = selectionScore;
				best.costRad2 = costAllObs;
				best.subarcStart = winningSubarcStart;
				best.subarcEnd = winningSubarcEnd;
			}
			return score;
		};

		auto tryCandidate = [&](const StateVectors& sv, double epoch, const char* name, bool verbose) -> double {
			return tryCandidateWithSelection(sv, epoch, name, verbose, kBadScore);
		};

		LOG_INFO << "[IOD] deterministic search: initialArclen=" << arclen << "d minArclen=" << kMinArclenDays
		         << "d shrink=" << kArclenShrink;

		int pass = 0;
		std::size_t lastStart = ~std::size_t{0};
		std::size_t lastEnd = ~std::size_t{0};

		while (best.score > kAcceptableScore && arclen > kMinArclenDays) {
			++pass;

			std::size_t start = 0;
			std::size_t end = 0;
			lookForBestSubarc(groundObs, arclen, start, end);

			if (end < start + 2) {
				LOG_INFO << "[IOD] pass " << pass << " arclen=" << arclen << "d: no valid subarc, shrinking";
				arclen *= kArclenShrink;
				continue;
			}

			const std::size_t n_sub = end - start + 1;
			const double subSpan = groundObs[end].epochJdTdb() - groundObs[start].epochJdTdb();

			if (start == lastStart && end == lastEnd) {
				LOG_INFO << "[IOD] pass " << pass << " arclen=" << arclen
				         << "d: same subarc as previous pass, skipping";
				arclen *= kArclenShrink;
				continue;
			}
			lastStart = start;
			lastEnd = end;

			winningSubarcStart = start;
			winningSubarcEnd = end;

			const std::vector<ObsNode> subarc(groundObs.begin() + start, groundObs.begin() + end + 1);

			subarcRef = &subarc;

			LOG_INFO << "[IOD] pass " << pass << " arclen-budget=" << arclen << "d  →  subarc [" << start << ".." << end
			         << "] (" << n_sub << " obs over " << subSpan << "d)";

			const std::size_t midIdx = pickGaussMiddleIdx(subarc);
			const std::vector<ObsNode> triplet = {subarc.front(), subarc[midIdx], subarc.back()};

			const auto gaussRoots = gaussAllRoots(triplet[0], triplet[1], triplet[2]);
			LOG_INFO << "[IOD]   Gauss surviving roots: " << gaussRoots.size();

			static const char* kGaussNames[6] = {
			    "Gauss[0]", "Gauss[0]+Herget-refined", "Gauss[1]", "Gauss[1]+Herget-refined",
			    "Gauss[2]", "Gauss[2]+Herget-refined",
			};
			static const char* kGaussImprovedNames[3] = {
			    "Gauss[0]+attempt_improvements",
			    "Gauss[1]+attempt_improvements",
			    "Gauss[2]+attempt_improvements",
			};

			for (std::size_t r = 0; r < gaussRoots.size(); ++r) {
				const StateVectors& sv = gaussRoots[r];
				const char* rawName = (r < 3) ? kGaussNames[2 * r] : "Gauss[?]";
				const double rawScore = tryCandidate(sv, triplet[1].epochJdTdb(), rawName, true);

				if (!std::isfinite(rawScore) || rawScore >= 1000.0)
					continue;

				const auto improved = attemptIodImprovements(sv, triplet[1].epochJdTdb(), subarc, rawScore,
				                                             weightProfile, fixedDistanceFitter, fullImprovementFitter);
				if (improved.changed) {
					const char* improvedName = (r < 3) ? kGaussImprovedNames[r] : "Gauss[?]+attempt_improvements";
					tryCandidateWithSelection(improved.state, improved.epoch, improvedName, true,
					                          improved.selectionScore);
				}

				const double epochTdb = triplet[1].epochJdTdb();
				auto orb = make_orbit_from_states_vectors(sv, epochTdb);
				if (!orb)
					continue;

				auto propagateRho = [&](const ObsNode& o) -> double {
					const Vec3D R = orb->getType() == OrbitType::Parabolic
					                    ? Calculations::propagateKeplerState(sv, o.epochJdTdb() - orb->getEpoch()).first
					                    : Calculations::positionAt(
					                          orb->getSemiMajorAxis(), orb->getExcentricity(), orb->getInclination(),
					                          orb->getLongitudeOfAscendingNode(), orb->getArgumentOfPericentre(),
					                          orb->getMeanAnomaly(), orb->getEpoch(), o.epochJdTdb(), Constants::GM_SUN,
					                          orb->getType());
					return (R + o.obsVector).magnitude();
				};

				const double rho_f = propagateRho(subarc.front());
				const double rho_b = propagateRho(subarc.back());
				if (!std::isfinite(rho_f) || !std::isfinite(rho_b) || rho_f <= 0.0 || rho_b <= 0.0)
					continue;

				const StateVectors refinedRaw = herget(subarc, rho_f, rho_b);
				const auto refined = fixedDistanceFitter.adjust(refinedRaw, subarc.front().epochJdTdb(), subarc);
				if (!refined.changed)
					continue;
				const char* refName = (r < 3) ? kGaussNames[2 * r + 1] : "Gauss[?]+Herget-refined";
				tryCandidate(refined.state, refined.epoch, refName, true);
			}

			{
				const auto adjusted = fixedDistanceFitter.adjust(herget(triplet), triplet.front().epochJdTdb(), subarc);
				if (adjusted.changed)
					tryCandidate(adjusted.state, adjusted.epoch, "Herget", true);
			}
			{
				const auto adjusted = fixedDistanceFitter.adjust(vaisala(triplet.front(), triplet.back()),
				                                                 triplet.front().epochJdTdb(), subarc);
				if (adjusted.changed)
					tryCandidate(adjusted.state, adjusted.epoch, "Vaisala", true);
			}

			if (subSpan < kMaxArclenForScansDays && best.score > kAcceptableScore) {
				int vaisalaCalls = 0;
				int hergetCalls = 0;
				double vaisalaBest = kBadScore;
				double hergetBest = kBadScore;

				{
					ScopedLogSilence quietScan(LogLevel::Error);

					{
						double pseudoR = kVaisalaScanRStart;
						bool reasonable = true;
						while (pseudoR < kVaisalaScanRMax && reasonable) {
							const auto raw = vaisala(subarc.front(), subarc.back(), pseudoR, pseudoR);
							const auto adjusted = fixedDistanceFitter.adjust(raw, subarc.front().epochJdTdb(), subarc);
							const double s = adjusted.changed
							                     ? tryCandidate(adjusted.state, adjusted.epoch, "Vaisala-scan", false)
							                     : kBadScore;
							++vaisalaCalls;
							if (s < vaisalaBest)
								vaisalaBest = s;
							if (s > kVaisalaFlippedOutScore && vaisalaBest < kBadScore)
								reasonable = false;
							pseudoR *= kScanStepFactor;
						}
					}

					{
						double pseudoR = kHergetScanRBase * std::pow(subSpan, kHergetScanRExpo);
						if (hasGeocentricObs)
							pseudoR += Constants::EARTH_RADIUS_IN_AU;

						while (pseudoR < kHergetScanRMax) {
							const auto raw = herget(subarc, pseudoR, pseudoR);
							const auto adjusted = fixedDistanceFitter.adjust(raw, subarc.front().epochJdTdb(), subarc);
							const double s = adjusted.changed
							                     ? tryCandidate(adjusted.state, adjusted.epoch, "Herget-scan", false)
							                     : kBadScore;
							++hergetCalls;
							if (s < hergetBest)
								hergetBest = s;
							pseudoR *= kScanStepFactor;
						}
					}
				}

				LOG_INFO << "[IOD]   Vaisala-scan " << vaisalaCalls << " calls (best " << vaisalaBest
				         << "), Herget-scan " << hergetCalls << " calls (best " << hergetBest << ")";
			}

			if (best.score < kImprovementTriggerScore && groundObs.size() > 2) {
				const double preScore = scoreIodState(best.state, best.epoch, subarc, weightProfile);
				const auto improved = attemptIodImprovements(best.state, best.epoch, subarc, preScore, weightProfile,
				                                             fixedDistanceFitter, fullImprovementFitter);
				if (improved.changed) {
					const double rawScore = tryCandidateWithSelection(
					    improved.state, improved.epoch, "attempt_improvements", true, improved.selectionScore);
					LOG_INFO << "[IOD]   attempt_improvements: score " << preScore << " -> " << rawScore
					         << " (selection " << improved.selectionScore << ")";
				}
			}

			arclen *= kArclenShrink;
		}

		if (best.score > kAcceptableScore && arclen <= kMinArclenDays) {
			LOG_INFO << "[IOD] subarc search collapsed below " << kMinArclenDays
			         << "d; using best-so-far score=" << best.score;
		}

		if (best.method != nullptr) {
			const double rms = std::sqrt(best.costRad2 / static_cast<double>(allObs.size())) / Constants::ARCSEC_RAD;

			const Vec3D& R = best.state.first;
			const Vec3D& V = best.state.second;
			const double ecc = eccentricityVec(R, V).magnitude();
			const Vec3D H = angMomentum(R, V);
			const double inclDeg = inclination(H) * Constants::RAD_TO_DEG;
			const double energy = specificEnergy(R, V, Constants::MU_IN_UNIT);
			const double a = (std::abs(energy) > 1e-20) ? semiMajorAxis(Constants::MU_IN_UNIT, energy) : 0.0;
			const double q = (ecc < 1.0) ? a * (1.0 - ecc) : std::abs(a) * (ecc - 1.0);

			LOG_INFO << "[IOD] Winner: " << best.method << "  RMS=" << rms << "\""
			         << "  score=" << best.score << "  (passes=" << pass << ")";
			LOG_INFO << "[IOD]   seed orbit: a=" << a << " AU"
			         << "  e=" << ecc << "  i=" << inclDeg << "°"
			         << "  q=" << q << " AU";
		} else {
			LOG_ERROR << "[IOD] All methods failed after " << pass << " passes";
		}

		return best;
	}
}  // namespace Calculations
