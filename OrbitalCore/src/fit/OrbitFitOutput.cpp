#include "fit/OrbitFitOutput.h"
#include "fit/FitCancellation.h"
#include "fit/PostFitResidual.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/Photometry.h"
#include "logger/Logger.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numbers>

namespace orbital_core::fit {
	namespace {
		void setPhotometricGeometry(Residual& row, const Vec3D& position, const Vec3D& observerHelio)
		{
			const Vec3D topo = position - observerHelio;
			row.helioDistAu = position.magnitude();
			row.geoDistAu = topo.magnitude();
			row.phaseAngleRad = std::numeric_limits<double>::quiet_NaN();
			if (row.helioDistAu > 1.0e-12 && row.geoDistAu > 1.0e-12) {
				const double observerDistance = observerHelio.magnitude();
				const double denom = 2.0 * row.helioDistAu * row.geoDistAu;
				const double cosine = (row.helioDistAu * row.helioDistAu + row.geoDistAu * row.geoDistAu -
				                       observerDistance * observerDistance) /
				                      denom;
				row.phaseAngleRad = std::acos(std::clamp(cosine, -1.0, 1.0));
			}
		}

		std::vector<double> postFitChiSquares(const OrbitFitModel& model,
		                                      const WeightedObservationLinearization& linear,
		                                      const Eigen::MatrixXd& covariance, std::span<const std::uint8_t> included)
		{
			std::vector<double> chi(model.observations().size(), std::numeric_limits<double>::quiet_NaN());
			for (std::size_t i = 0; i < chi.size(); ++i) {
				const int row = model.residualOffset(i);
				const int dimensions = model.observations()[i].measurementDimensions();
				const Eigen::MatrixXd jac = linear.jacobian.middleRows(row, dimensions);
				const Eigen::MatrixXd projected = jac * covariance * jac.transpose();
				const bool active = included.empty() || included[i] != 0;
				if (dimensions == 1) {
					const double variance = 1.0 + (active ? -projected(0, 0) : projected(0, 0));
					if (variance > 0.0 && std::isfinite(variance))
						chi[i] = linear.residuals[row] * linear.residuals[row] / variance;
				} else {
					const double scale = model.usesWeights() ? 1.0 : 1.0 / Constants::ARCSEC_RAD;
					const auto statistic = postFitResidualStatistic(
					    scale * Eigen::Vector2d(linear.residuals[row], linear.residuals[row + 1]), projected, active);
					if (statistic.valid)
						chi[i] = statistic.chiSquared;
				}
			}
			return chi;
		}
	}  // namespace

	std::vector<double> OrbitFitOutput::selectionChiSquares(const std::vector<double>& params) const
	{
		std::vector<double> chi(model_.observations().size(), std::numeric_limits<double>::quiet_NaN());
		const auto linear = problem_.weightedLinearization(params);
		if (!linear.valid)
			return chi;
		const auto& snapshot = problem_.snapshot();
		if (!snapshot.matches(params))
			return chi;
		Eigen::MatrixXd covariance =
		    uncertainty_.computeFitCovariance(params, snapshot.rawResiduals, false).formalFullCovariance;
		for (int i = 3; i < 6; ++i) {
			covariance.row(i) *= Constants::GAUSS_K;
			covariance.col(i) *= Constants::GAUSS_K;
		}
		if (!covariance.allFinite())
			return chi;
		return postFitChiSquares(model_, linear, covariance, problem_.included());
	}

	double OrbitFitOutput::reportRms(std::span<const double> rawResiduals) const
	{
		if (rawResiduals.empty())
			return 0.0;
		double sumSq = 0.0;
		int count = 0;
		for (std::size_t observation = 0; observation < model_.observations().size(); ++observation) {
			if (!problem_.isIncluded(observation) || model_.observations()[observation].radar)
				continue;
			const auto i = model_.residualOffset(observation);
			const double ra = Convertors::radiansToDegrees(rawResiduals[i]) * 3600.0;
			const double dec = Convertors::radiansToDegrees(rawResiduals[i + 1]) * 3600.0;
			sumSq += ra * ra + dec * dec;
			count += 2;
		}
		return count ? std::sqrt(sumSq / static_cast<double>(count)) : 0.0;
	}

	Residuals OrbitFitOutput::populateResidualsFromRaw(const std::vector<double>& rawRes,
	                                                   const std::vector<double>& params) const
	{
		model_.nonGravProfile().updateFromParams(params, 6);
		const std::vector<Vec3D> positions = model_.computeLightTimePositions(params, model_.propagator().nonGrav());

		Residuals residuals;
		for (std::size_t i = 0; i < model_.observations().size(); ++i) {
			const auto& obs = model_.observations()[i];
			int offset = model_.residualOffset(i);
			if (obs.radar) {
				Residual row{};
				row.jd = obs.jdTdb;
				row.obsCode = obs.observer.obsCode;
				row.raObs = row.decObs = row.dRA_arcsec = row.dDEC_arcsec = row.total_arcsec =
				    std::numeric_limits<double>::quiet_NaN();
				row.radar.emplace();
				row.radar->transmitter = obs.radar->transmitter;
				row.radar->frequencyHz = obs.radar->frequencyHz;
				const auto channel = [&](const auto& measurement) -> std::optional<Residual::RadarChannel> {
					if (!measurement)
						return std::nullopt;
					const double residual = rawRes[offset++];
					const double normalized = residual / measurement->sigma;
					row.radar->chiSquared += normalized * normalized;
					return Residual::RadarChannel{measurement->value, measurement->value - residual, residual,
					                              measurement->sigma};
				};
				row.radar->delaySeconds = channel(obs.radar->delaySeconds);
				row.radar->dopplerHz = channel(obs.radar->dopplerHz);
				residuals.add(row);
				continue;
			}
			const double dRA_as = Convertors::radiansToDegrees(rawRes[offset]) * 3600.0;
			const double dDEC_as = Convertors::radiansToDegrees(rawRes[offset + 1]) * 3600.0;
			const double total = std::sqrt(dRA_as * dRA_as + dDEC_as * dDEC_as);

			Residual row;
			row.jd = obs.jdTdb;
			row.raObs = obs.raObs;
			row.decObs = obs.decObs;
			row.dRA_arcsec = dRA_as;
			row.dDEC_arcsec = dDEC_as;
			row.total_arcsec = total;
			row.obsCode = obs.observer.obsCode;
			row.mag = obs.mag;
			row.magBand = obs.magBand;
			row.nuclearMagnitude = obs.nuclearMagnitude;
			row.magSigma = obs.magSigma;
			const bool hasFinalTransform = problem_.snapshot().matches(params) &&
			                               problem_.snapshot().transforms.size() == model_.observations().size();
			row.astrometricWeight = orbital_core::obs::residualWeightForReport(
			    obs, rawRes[offset], rawRes[offset + 1], model_.usesWeights(),
			    hasFinalTransform ? problem_.snapshot().transforms[i]
			                      : orbital_core::obs::AstrometricResidualTransform{});

			if (i < positions.size()) {
				const Vec3D observerHelio = obs.earthR + obs.observerOffsetEclAu;
				setPhotometricGeometry(row, positions[i], observerHelio);
			}

			residuals.add(row);
		}
		return residuals;
	}

	Orbit_p OrbitFitOutput::propagateToTargetEpoch(const std::vector<double>& finalParams, const Orbit& fitOrbit) const
	{
		if (config_.outputEpochTdb <= 0.0)
			return nullptr;

		const double outputEpochTdb = config_.outputEpochTdb;

		if (std::abs(outputEpochTdb - fitOrbit.getEpoch()) <= 1e-6)
			return nullptr;

		const Vec3D R0(finalParams[0], finalParams[1], finalParams[2]);
		const Vec3D V0_canonical(finalParams[3], finalParams[4], finalParams[5]);
		const Vec3D V0_real = V0_canonical * Constants::GAUSS_K;

		const auto stv =
		    model_.propagator().propagateAllWithVelocity(R0, V0_real, fitOrbit.getEpoch(), {outputEpochTdb});
		if (stv.empty())
			return nullptr;

		const auto& sv = stv.front();
		const Vec3D V_canonical = sv.V * (1.0 / Constants::GAUSS_K);
		auto orbAt = make_orbit_from_states_vectors(StateVectors(sv.R, V_canonical), outputEpochTdb);

		if (config_.logDiagnostics)
			LOG_DEBUG << "[FitOutput] propagated orbit to JD " << std::fixed << std::setprecision(6) << outputEpochTdb
			          << " TDB";

		return orbAt;
	}

	FitResult OrbitFitOutput::build(const std::vector<double>& params, const FitSolveSummary& summary,
	                                std::stop_token cancel) const
	{
		checkCancellation(cancel);
		FitResult result;
		result.iterations = summary.iterations;
		result.converged = summary.converged;
		result.cancelled = summary.cancelled;
		auto fitOrbit = model_.orbitFromStateVector(params);
		std::vector<double> rawRes;
		const auto& snapshot = problem_.snapshot();
		const bool hasFinalResiduals =
		    snapshot.matches(params) && snapshot.rawResiduals.size() == model_.numResiduals();
		if (hasFinalResiduals) {
			rawRes = snapshot.rawResiduals;
		} else {
			rawRes = model_.computeResiduals(params);
		}

		checkCancellation(cancel);
		if (model_.usesWeights() && !rawRes.empty()) {
			auto weightedResiduals = rawRes;
			problem_.weightResiduals(weightedResiduals,
			                         hasFinalResiduals
			                             ? std::span<const obs::AstrometricResidualTransform>(snapshot.transforms)
			                             : std::span<const obs::AstrometricResidualTransform>{});
			problem_.maskResiduals(weightedResiduals);
			if (problem_.activeResidualCount() > 0)
				result.normalizedRms = std::sqrt(problem_.sumSquares(weightedResiduals) /
				                                 static_cast<double>(problem_.activeResidualCount()));
		}
		if (config_.populateResiduals)
			result.residuals = populateResidualsFromRaw(rawRes, params);

		checkCancellation(cancel);
		result.fitEpochOrbit = fitOrbit;

		finalize(result, params, rawRes, cancel);
		result.cost = config_.populateResiduals ? result.residuals.rms() : reportRms(rawRes);
		result.converged = result.orbit && summary.converged;
		result.finalNonGrav = model_.propagator().nonGrav();
		return result;
	}

	void OrbitFitOutput::finalize(FitResult& result, const std::vector<double>& params,
	                              const std::vector<double>& rawRes, std::stop_token cancel) const
	{
		checkCancellation(cancel);
		if (config_.computeUncertainty) {
			populateUncertainty(result, params, rawRes, cancel);
		} else {
			result.orbit = propagateToTargetEpoch(params, *result.fitEpochOrbit);
			if (!result.orbit)
				result.orbit = result.fitEpochOrbit;
		}

		checkCancellation(cancel);
	}

	PhotometricOrbitSigmas OrbitFitOutput::photometricOrbitSigmas(const Residuals& nominal,
	                                                              const std::vector<double>& params,
	                                                              const Eigen::MatrixXd& covariance,
	                                                              std::stop_token cancel) const
	{
		PhotometricOrbitSigmas result;
		const auto n = static_cast<Eigen::Index>(params.size());
		const auto& observations = model_.observations();
		if (n < 6 || covariance.rows() != n || covariance.cols() != n || !covariance.allFinite() ||
		    nominal.count() != observations.size())
			return result;

		const auto asteroid = photometry::fit(nominal.data(), config_.photometry.asteroidG);
		const auto comet = photometry::fitComet(nominal.data(), config_.photometry);
		if (!asteroid.absoluteMagnitude && !comet.total.absoluteMagnitude && !comet.nuclear.absoluteMagnitude)
			return result;

		try {
			checkCancellation(cancel);
			const auto derivatives = model_.computePhotometricPositionDerivatives(params, cancel);
			if (derivatives.positions.size() != observations.size() ||
			    derivatives.jacobians.size() != observations.size())
				return result;

			struct Gradient {
				Eigen::VectorXd value;
				double weightSum{0.0};
				std::size_t count{0};
				bool valid{true};
			};

			Gradient h{Eigen::VectorXd::Zero(n)}, m1{Eigen::VectorXd::Zero(n)}, m2{Eigen::VectorXd::Zero(n)};
			const auto& rows = nominal.data();
			constexpr double invLog10 = 1.0 / std::numbers::ln10;
			for (std::size_t i = 0; i < rows.size(); ++i) {
				checkCancellation(cancel);
				const auto hp = asteroid.absoluteMagnitude ? photometry::evaluate(rows[i], asteroid, true)
				                                           : photometry::Evaluation{};
				const auto cp = rows[i].nuclearMagnitude ? photometry::evaluateComet(rows[i], comet, true)
				                                         : photometry::Evaluation{};
				Gradient* cg = cp.residual ? (*rows[i].nuclearMagnitude ? &m2 : &m1) : nullptr;
				if (!hp.residual && !cg)
					continue;
				const auto fail = [&] {
					if (hp.residual)
						h.valid = false;
					if (cg)
						cg->valid = false;
				};
				if (observations[i].radar || derivatives.jacobians[i].rows() != 3 ||
				    derivatives.jacobians[i].cols() != n) {
					fail();
					continue;
				}

				const auto& p = derivatives.positions[i];
				const Vec3D observer = observations[i].earthR + observations[i].observerOffsetEclAu;
				const Eigen::Vector3d r(p.x, p.y, p.z);
				const Eigen::Vector3d d(p.x - observer.x, p.y - observer.y, p.z - observer.z);
				const double rnorm = r.norm(), dnorm = d.norm();
				if (!(rnorm > 0.0) || !(dnorm > 0.0) || !std::isfinite(rnorm) || !std::isfinite(dnorm)) {
					fail();
					continue;
				}
				const Eigen::Vector3d rGradient = r / (rnorm * rnorm);
				const Eigen::Vector3d dGradient = d / (dnorm * dnorm);
				const double alpha = rows[i].phaseAngleRad;
				const double sine = std::sin(alpha);
				const auto cosineGradient = [&]() -> Eigen::Vector3d {
					return (r + d) / (rnorm * dnorm) - std::cos(alpha) * (rGradient + dGradient);
				};
				if (hp.residual) {
					if (!std::isfinite(alpha) || !(sine > 1e-12)) {
						h.valid = false;
					} else {
						const auto [phase, phaseRate] =
						    Calculations::hgPhaseFunctionWithDerivative(alpha, config_.photometry.asteroidG);
						if (!(phase > 0.0) || !std::isfinite(phaseRate)) {
							h.valid = false;
						} else {
							const Eigen::Vector3d gradient =
							    -5.0 * invLog10 * (rGradient + dGradient) -
							    2.5 * invLog10 * phaseRate / (phase * sine) * cosineGradient();
							h.value += derivatives.jacobians[i].transpose() * gradient / hp.sigma;
							h.weightSum += 1.0 / hp.sigma;
							++h.count;
						}
					}
				}
				if (cg) {
					const auto& component = *rows[i].nuclearMagnitude ? comet.nuclear : comet.total;
					Eigen::Vector3d gradient = -invLog10 * (component.slope * rGradient + 5.0 * dGradient);
					if (component.phaseCoefficient != 0.0) {
						if (!std::isfinite(alpha) || !(sine > 1e-12)) {
							cg->valid = false;
							continue;
						}
						gradient += component.phaseCoefficient * Constants::RAD_TO_DEG / sine * cosineGradient();
					}
					cg->value += derivatives.jacobians[i].transpose() * gradient / cp.sigma;
					cg->weightSum += 1.0 / cp.sigma;
					++cg->count;
				}
			}

			const auto sigma = [&](const Gradient& gradient, std::size_t expected) -> std::optional<double> {
				if (!gradient.valid || expected == 0 || gradient.count != expected || !(gradient.weightSum > 0.0))
					return std::nullopt;
				const Eigen::VectorXd mean = gradient.value / gradient.weightSum;
				const double variance = mean.dot(covariance * mean);
				return std::isfinite(variance) && variance >= 0.0 ? std::optional<double>{std::sqrt(variance)}
				                                                  : std::nullopt;
			};
			result.asteroid = sigma(h, asteroid.count);
			result.cometTotal = sigma(m1, comet.total.count);
			result.cometNuclear = sigma(m2, comet.nuclear.count);
			return result;
		} catch (const FitCancelled&) {
			throw;
		} catch (const EphemerisError&) {
			throw;
		} catch (const std::exception& e) {
			checkCancellation(cancel);
			LOG_WARN << "[Photometry] Orbit uncertainty unavailable: " << e.what();
			return result;
		}
	}
	void OrbitFitOutput::populateUncertainty(FitResult& result, const std::vector<double>& params,
	                                         const std::vector<double>& rawRes, std::stop_token cancel) const
	{
		const FitCovariance fitCov = uncertainty_.computeFitCovariance(params, rawRes);
		checkCancellation(cancel);
		if (!fitCov.fullCovariance.allFinite() || !fitCov.formalFullCovariance.allFinite()) {
			result.orbit = propagateToTargetEpoch(params, *result.fitEpochOrbit);
			if (!result.orbit)
				result.orbit = result.fitEpochOrbit;
			return;
		}

		const int nParams = static_cast<int>(fitCov.fullCovariance.rows());
		Eigen::MatrixXd paramScale = Eigen::MatrixXd::Identity(nParams, nParams);
		if (nParams >= 6) {
			paramScale(3, 3) = Constants::GAUSS_K;
			paramScale(4, 4) = Constants::GAUSS_K;
			paramScale(5, 5) = Constants::GAUSS_K;
		}
		result.uncertainty.fitParameterCovariance = paramScale * fitCov.fullCovariance * paramScale.transpose();
		result.uncertainty.fitParameterCovarianceFormal =
		    paramScale * fitCov.formalFullCovariance * paramScale.transpose();
		result.uncertainty.cartesianStateCovariance = result.uncertainty.fitParameterCovariance.topLeftCorner(6, 6);
		result.uncertainty.cartesianStateCovarianceFormal =
		    result.uncertainty.fitParameterCovarianceFormal.topLeftCorner(6, 6);
		result.uncertainty.covarianceObjectiveScale =
		    (std::isfinite(fitCov.varianceFactor) && fitCov.varianceFactor > 0.0) ? fitCov.varianceFactor : 1.0;

		const auto photometrySigmas = photometricOrbitSigmas(
		    result.allResiduals.count() == model_.observations().size() ? result.allResiduals : result.residuals,
		    params, fitCov.fullCovariance, cancel);
		result.uncertainty.absoluteMagnitudeOrbitSigma = photometrySigmas.asteroid;
		result.uncertainty.cometTotalMagnitudeOrbitSigma = photometrySigmas.cometTotal;
		result.uncertainty.cometNuclearMagnitudeOrbitSigma = photometrySigmas.cometNuclear;

		auto orbit = propagateToTargetEpoch(params, *result.fitEpochOrbit);
		checkCancellation(cancel);
		if (!orbit)
			orbit = result.fitEpochOrbit;

		auto saveOutputStateCovariance = [&](const Eigen::Matrix<double, 6, Eigen::Dynamic>& stateJ) {
			if (stateJ.cols() != fitCov.fullCovariance.rows() || !stateJ.allFinite())
				return;
			Eigen::Matrix<double, 6, Eigen::Dynamic> physicalJ = stateJ;
			physicalJ.bottomRows<3>() *= Constants::GAUSS_K;
			const StateCov6x6 covariance = physicalJ * fitCov.fullCovariance * physicalJ.transpose();
			if (covariance.allFinite()) {
				result.uncertainty.outputStateCovariance = covariance;
				result.uncertainty.outputStateCovarianceEpochJdTdb = orbit->getEpoch();
			}
		};
		const bool useCometaryElements = cometary::usesCometaryElements(orbit->getExcentricity());
		if (useCometaryElements) {
			const auto derivatives =
			    uncertainty_.computeCometaryDerivatives(params, *result.fitEpochOrbit, *orbit, cancel);
			checkCancellation(cancel);
			saveOutputStateCovariance(derivatives.stateJ);
			uncertainty_.populateCometaryUncertainty(result.uncertainty, derivatives, fitCov, *orbit);
		} else {
			const auto derivs = uncertainty_.computeElementDerivatives(params, *result.fitEpochOrbit, cancel);
			checkCancellation(cancel);
			saveOutputStateCovariance(derivs.stateJ);
			const Eigen::Matrix<double, 6, 6> covOrb = derivs.elemJ * fitCov.fullCovariance * derivs.elemJ.transpose();
			const Eigen::MatrixXd qualityElementCov =
			    derivs.qualityElemJ * fitCov.fullCovariance * derivs.qualityElemJ.transpose();
			result.uncertainty.elementCovariance = covOrb;
			result.uncertainty.hasElementCovariance = true;
			result.uncertainty.derivedSigmas = uncertainty_.computeDerivedSigmasFromQualityCovariance(
			    *orbit, result.uncertainty.elementCovariance, qualityElementCov);
		}
		const auto& sig = result.uncertainty.derivedSigmas;
		result.uncertainty.orbitQualityU = uncertainty_.orbitQualityU(*orbit, sig);

		uncertainty_.populateOrbitElementSigmas(*orbit, sig);

		if (config_.logDiagnostics)
			uncertainty_.logOrbitQualityU(*orbit, sig);

		result.orbit = orbit;
	}

	PostFitDiagnostics OrbitFitOutput::diagnostics(const FitResult& result,
	                                               const orbital_core::FitParameterCovariance& formalCovariance,
	                                               std::span<const std::uint8_t> included) const
	{
		PostFitDiagnostics out;
		out.residuals = result.residuals;
		const auto& snapshot = problem_.snapshot();
		if (!snapshot.valid || snapshot.params.size() < 6 ||
		    (!included.empty() && included.size() != model_.observations().size()))
			return out;
		const auto linear = problem_.weightedLinearization(snapshot.params);
		if (!linear.valid || out.residuals.count() != model_.observations().size())
			return out;

		if (formalCovariance.size() != 0) {
			out.formalCovariance = formalCovariance;
		} else if (result.uncertainty.fitParameterCovarianceFormal.size() != 0) {
			out.formalCovariance = result.uncertainty.fitParameterCovarianceFormal;
		} else {
			const auto covariance = uncertainty_.computeFitCovariance(snapshot.params, snapshot.rawResiduals);
			out.formalCovariance = covariance.formalFullCovariance;
			for (int i = 3; i < 6; ++i) {
				out.formalCovariance.row(i) *= Constants::GAUSS_K;
				out.formalCovariance.col(i) *= Constants::GAUSS_K;
			}
		}
		if (out.formalCovariance.rows() != linear.jacobian.cols() ||
		    out.formalCovariance.cols() != linear.jacobian.cols() || !out.formalCovariance.allFinite())
			return out;

		Residuals corrected;
		std::size_t invalid = 0;
		const auto chi = postFitChiSquares(model_, linear, out.formalCovariance, included);
		for (std::size_t i = 0; i < model_.observations().size(); ++i) {
			Residual row = out.residuals.data()[i];
			const bool valid = std::isfinite(chi[i]);
			if (row.astrometricWeight) {
				row.astrometricWeight->postFitChiSquared.reset();
				if (valid)
					row.astrometricWeight->postFitChiSquared = chi[i];
			}
			if (!valid)
				++invalid;
			if (row.radar) {
				row.radar->postFitChiSquared.reset();
				if (valid)
					row.radar->postFitChiSquared = chi[i];
			}
			corrected.add(row);
		}
		if (invalid)
			LOG_WARN << "[PostFitResiduals] " << invalid
			         << " singular/invalid residual covariances; no rejection for these rows";
		out.residuals = std::move(corrected);
		out.valid = true;
		return out;
	}

}  // namespace orbital_core::fit
