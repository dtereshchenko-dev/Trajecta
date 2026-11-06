#include "fit/OrbitFitModel.h"
#include "TwoBodyPropagation.h"
#include "orbital_core/Types.h"
#include "obs/LightDeflection.h"
#include "obs/LightTime.h"
#include "obs/RadarPrediction.h"
#include "logger/Logger.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace orbital_core::fit {
	using obs::PreparedObservation;

	OrbitFitModel::OrbitFitModel(Orbit_p initial, std::vector<PreparedObservation> observations,
	                             const OrbitOptimizerConfig& config)
	    : initialOrbit_(initial ? std::make_shared<Orbit>(*initial) : nullptr),
	      observations_(std::move(observations)),
	      cosObliquity_(std::cos(Constants::EARTH_INCL)),
	      sinObliquity_(std::sin(Constants::EARTH_INCL)),
	      useWeights_(config.useWeights),
	      useNBody_(config.nbody.enabled || config.nonGrav.enabled)
	{
		if (!initialOrbit_)
			throw std::invalid_argument("OrbitFitModel requires an initial orbit");

		configureDynamics(config);
		for (const auto& observation : observations_) {
			residualOffsets_.push_back(residualCount_);
			residualCount_ += observation.measurementDimensions();
			if (observation.radar && !useWeights_)
				throw std::invalid_argument("Mixed radar/optical fitting requires uncertainty weighting");
		}
	}

	void OrbitFitModel::configureDynamics(const OrbitOptimizerConfig& config)
	{
		nbody_.setStepSize(config.nbody.stepSize);
		nbody_.setPlanetMask(config.nbody.enabled ? config.nbody.planetMask : 0);
		nbody_.setUseSb441(config.nbody.useSb441);
		nbody_.setTargetNaifId(config.nbody.targetNaifId);
		nbody_.setInactiveInteriorMassPolicy(config.nbody.inactiveInteriorMassPolicy);
		if (useNBody_) {
			if (auto* ephemeris = getDefaultEphemerisSource())
				nbody_ = nbody_.forEncounterMonitoring(*ephemeris);
		}

		nbody_.nonGrav().setEnabled(config.nonGrav.enabled);
		if (config.nonGrav.enabled) {
			nbody_.nonGrav().setModel(static_cast<NgModel>(config.nonGrav.model));
			const auto requestedMask = static_cast<NgParameterMask>(config.nonGrav.solveMask);
			nbody_.nonGrav().setFitParameterMask(
			    requestedMask == NgParameterMask::None
			        ? NonGravModel::defaultFitParameterMaskFor(static_cast<NgModel>(config.nonGrav.model))
			        : requestedMask);
			if (config.nonGrav.useInitialValues) {
				nbody_.nonGrav().setParams(config.nonGrav.initialA1, config.nonGrav.initialA2,
				                           config.nonGrav.initialA3);
				nbody_.nonGrav().setDT(config.nonGrav.initialDT);

				if (config.logDiagnostics) {
					LOG_INFO << "Non-grav initial values: A1=" << std::scientific << std::setprecision(6)
					         << config.nonGrav.initialA1 << " A2=" << config.nonGrav.initialA2
					         << " A3=" << config.nonGrav.initialA3 << " DT=" << config.nonGrav.initialDT
					         << " skipBootstrap=" << (config.skipNonGravBootstrap ? "true" : "false");
				}
			}
		}
	}

	Vec3D OrbitFitModel::toTopocentric(const Vec3D& R_helio, const PreparedObservation& obs) const
	{
		const Vec3D R_topo = orbital_core::obs::opticalObservationDirection(R_helio, obs);

		Vec3D R_eq(R_topo.x, R_topo.y * cosObliquity_ - R_topo.z * sinObliquity_,
		           R_topo.y * sinObliquity_ + R_topo.z * cosObliquity_);

		return R_eq;
	}

	std::pair<double, double> OrbitFitModel::toRADEC(const Vec3D& R_eq) noexcept
	{
		const double rho = R_eq.magnitude();
		if (rho <= 0.0 || !std::isfinite(rho))
			return {std::nan(""), std::nan("")};

		double ra = std::atan2(R_eq.y, R_eq.x);
		if (ra < 0.0)
			ra += 2.0 * Constants::PI;

		double dec = std::asin(std::clamp(R_eq.z / rho, -1.0, 1.0));
		return {ra, dec};
	}

	std::vector<double> OrbitFitModel::computeResiduals(const std::vector<double>& params) const
	{
		nonGravProfile_.updateFromParams(params, 6);
		return computeResiduals(params, nbody_.nonGrav());
	}

	OrbitFitModel::ResidualsAndStateJac OrbitFitModel::computeResidualsAndStateJac(
	    const std::vector<double>& params) const
	{
		ResidualsAndStateJac out;
		const int N = static_cast<int>(observations_.size());

		if (!useNBody_) {
			out.residuals = computeResiduals(params);
			return out;
		}

		nonGravProfile_.updateFromParams(params, 6);
		const NonGravModel& ng = nbody_.nonGrav();

		const int nModelParams = ng.numFitParams();
		const int nAcols = std::clamp(static_cast<int>(params.size()) - 6, 0, nModelParams);
		const int jacCols = 6 + nAcols;
		constexpr int kMaxJacCols = 10;
		if (jacCols > kMaxJacCols) {
			out.residuals = computeResiduals(params);
			return out;
		}

		auto orbit = orbitFromStateVector(params);
		const Vec3D R0 = orbit->getR();
		const Vec3D V0_real = orbit->getV() * Constants::GAUSS_K;
		const double epochMjd = orbit->epochMjd();

		std::vector<double> receptionMjds(N);
		for (int i = 0; i < N; ++i)
			receptionMjds[i] = observations_[i].mjdTdb;

		std::vector<NBodyPropagator::PropagationResultFull> fullResults;
		std::vector<NBodyPropagator::PropagationResult> stmResults;
		if (nAcols > 0)
			fullResults = nbody_.propagateAllWithStmAndSensMjd(R0, V0_real, epochMjd, receptionMjds, ng);
		else
			stmResults = nbody_.propagateAllWithStmMjd(R0, V0_real, epochMjd, receptionMjds, ng);

		if ((nAcols > 0 && fullResults.size() != static_cast<std::size_t>(N)) ||
		    (nAcols == 0 && stmResults.size() != static_cast<std::size_t>(N))) {
			out.residuals = computeResiduals(params);
			return out;
		}

		out.residuals.resize(residualCount_);
		out.stateJ.resize(residualCount_, jacCols);
		out.transforms.resize(static_cast<std::size_t>(N));
		std::vector<Vec3D> emissionPositions(static_cast<std::size_t>(N));
		std::vector<Vec3D> emissionVelocities(static_cast<std::size_t>(N));
		std::vector<Eigen::Matrix<double, 3, 6>> emissionJacobians(static_cast<std::size_t>(N));

		for (int i = 0; i < N; ++i) {
			const Vec3D& rReception =
			    nAcols > 0 ? fullResults[static_cast<std::size_t>(i)].R : stmResults[static_cast<std::size_t>(i)].R;
			const Vec3D& vReception =
			    nAcols > 0 ? fullResults[static_cast<std::size_t>(i)].V : stmResults[static_cast<std::size_t>(i)].V;
			const auto& observation = observations_[i];
			if (observation.radar) {
				auto* eph = getDefaultEphemerisSource();
				if (!eph)
					throw std::runtime_error("Radar requires planetary ephemerides");
				const auto prediction = obs::predictRadar(observation, rReception, vReception, *eph);
				int row = residualOffset(i);
				if (observation.radar->delaySeconds)
					out.residuals[row++] = observation.radar->delaySeconds->value - prediction.delaySeconds;
				if (observation.radar->dopplerHz)
					out.residuals[row] = observation.radar->dopplerHz->value - prediction.dopplerHz;
				Eigen::MatrixXd transport = Eigen::MatrixXd::Zero(6, jacCols);
				transport.leftCols(6) = nAcols > 0 ? fullResults[i].Phi : stmResults[i].Phi;
				transport.middleCols(3, 3) *= Constants::GAUSS_K;
				for (int col = 0; col < nAcols; ++col)
					transport.col(6 + col) = fullResults[i].SensA.col(ng.sensitivityColumnForFitParameter(col));
				out.stateJ.middleRows(residualOffset(i), observation.measurementDimensions()) =
				    obs::radarResidualJacobian(observation, rReception, vReception, *eph) * transport;
				continue;
			}
			const auto apparent = orbital_core::obs::solveLightTimeTwoPass(
			    observations_[static_cast<std::size_t>(i)].jdTdb, rReception, vReception,
			    observations_[static_cast<std::size_t>(i)].earthR +
			        observations_[static_cast<std::size_t>(i)].observerOffsetEclAu,
			    &emissionJacobians[static_cast<std::size_t>(i)]);
			if (!apparent.valid) {
				out.residuals = computeResiduals(params);
				out.stateJ.resize(0, 0);
				out.transforms.clear();
				emissionPositions.clear();
				return out;
			}
			emissionPositions[static_cast<std::size_t>(i)] = apparent.position;
			emissionVelocities[static_cast<std::size_t>(i)] = apparent.velocityAuPerDay;
		}

#pragma omp parallel for
		for (int i = 0; i < N; ++i) {
			const auto& obs = observations_[i];
			if (obs.radar)
				continue;
			const int row = residualOffset(i);
			const Vec3D& r_em = emissionPositions[static_cast<std::size_t>(i)];
			const Vec3D& v_em = emissionVelocities[static_cast<std::size_t>(i)];

			Vec3D R_eq = toTopocentric(r_em, obs);
			auto [raModel, decModel] = toRADEC(R_eq);
			const double dRA = Helpers::wrapPi(obs.raObs - raModel);
			const double dDEC = Helpers::wrapPi(obs.decObs - decModel);
			const double cosDecObs = obs.raMetricCosDec;
			out.residuals[row] = cosDecObs * dRA;
			out.residuals[row + 1] = dDEC;
			out.transforms[static_cast<std::size_t>(i)] =
			    orbital_core::obs::astrometricResidualTransform(obs, r_em, &v_em);

			const Eigen::Matrix<double, 2, 3> dRes_dRem =
			    orbital_core::obs::astrometricAngularJacobianEcl(R_eq, cosDecObs) *
			    orbital_core::obs::opticalObservationDirectionJacobian(r_em, obs);

			Eigen::Matrix<double, 6, kMaxJacCols> transport = Eigen::Matrix<double, 6, kMaxJacCols>::Zero();
			transport.leftCols<6>() = nAcols > 0 ? fullResults[i].Phi : stmResults[i].Phi;
			transport.middleCols<3>(3) *= Constants::GAUSS_K;
			if (nAcols > 0) {
				for (int col = 0; col < nAcols; ++col)
					transport.col(6 + col) = fullResults[i].SensA.col(ng.sensitivityColumnForFitParameter(col));
			}

			out.stateJ.block(row, 0, 2, jacCols) =
			    dRes_dRem * emissionJacobians[static_cast<std::size_t>(i)] * transport.leftCols(jacCols);
		}

		if (!out.stateJ.allFinite())
			out.stateJ.resize(0, 0);
		return out;
	}

	std::vector<double> OrbitFitModel::computeResiduals(const std::vector<double>& params, const NonGravModel& ng) const
	{
		const auto states = receptionStates(params, ng);

		const int N = static_cast<int>(observations_.size());
		std::vector<double> res(residualCount_, 1.0e30);
		if (states.size() != observations_.size())
			return res;

		for (int i = 0; i < N; ++i) {
			const auto& obs = observations_[i];
			int row = residualOffset(i);
			if (obs.radar) {
				auto* eph = getDefaultEphemerisSource();
				if (!eph)
					throw std::runtime_error("Radar requires planetary ephemerides");
				const auto prediction = orbital_core::obs::predictRadar(obs, states[i].R, states[i].V, *eph);
				if (obs.radar->delaySeconds)
					res[row++] = obs.radar->delaySeconds->value - prediction.delaySeconds;
				if (obs.radar->dopplerHz)
					res[row] = obs.radar->dopplerHz->value - prediction.dopplerHz;
				continue;
			}
			const auto apparent = orbital_core::obs::solveLightTimeTwoPass(obs.jdTdb, states[i].R, states[i].V,
			                                                               obs.earthR + obs.observerOffsetEclAu);
			if (!apparent.valid)
				continue;
			Vec3D R_eq = toTopocentric(apparent.position, obs);
			auto [raModel, decModel] = toRADEC(R_eq);

			const double dRA = Helpers::wrapPi(obs.raObs - raModel);
			const double dDEC = Helpers::wrapPi(obs.decObs - decModel);

			res[row] = obs.raMetricCosDec * dRA;
			res[row + 1] = dDEC;
		}
		return res;
	}

	std::vector<NBodyPropagator::StateVel> OrbitFitModel::receptionStates(const std::vector<double>& params,
	                                                                      const NonGravModel& ng) const
	{
		auto orbit = orbitFromStateVector(params);

		std::vector<double> receptionMjds(observations_.size());
		for (size_t i = 0; i < observations_.size(); ++i)
			receptionMjds[i] = observations_[i].mjdTdb;

		std::vector<NBodyPropagator::StateVel> receptionStates;
		if (useNBody_) {
			const Vec3D velocityAuPerDay = orbit->getV() * Constants::GAUSS_K;
			receptionStates = nbody_.propagateAllWithVelocityMjd(orbit->getR(), velocityAuPerDay, orbit->epochMjd(),
			                                                     receptionMjds, ng);
		}

		else {
			receptionStates.resize(observations_.size());
			const StateVectors initial(orbit->getR(), orbit->getV());
			for (size_t i = 0; i < observations_.size(); ++i) {
				const StateVectors state =
				    Calculations::propagateKeplerState(initial, receptionMjds[i] - orbit->epochMjd());
				if (!state.isValid())
					return {};
				receptionStates[i] = {state.first, state.second * Constants::GAUSS_K};
			}
		}

		if (receptionStates.size() != observations_.size())
			return {};
		return receptionStates;
	}

	std::vector<Vec3D> OrbitFitModel::computeLightTimePositions(const std::vector<double>& params,
	                                                            const NonGravModel& ng) const
	{
		const auto states = receptionStates(params, ng);
		if (states.size() != observations_.size())
			return {};
		std::vector<Vec3D> positions(observations_.size());
		for (size_t i = 0; i < observations_.size(); ++i) {
			if (observations_[i].radar) {
				positions[i] = states[i].R;
				continue;
			}
			const auto apparent = orbital_core::obs::solveLightTimeTwoPass(
			    observations_[i].jdTdb, states[i].R, states[i].V,
			    observations_[i].earthR + observations_[i].observerOffsetEclAu);
			if (!apparent.valid)
				return {};
			positions[i] = apparent.position;
		}
		return positions;
	}

	double OrbitFitModel::computeCost(const std::vector<double>& params) const
	{
		return opticalRms(computeResiduals(params));
	}

	double OrbitFitModel::opticalRms(const std::vector<double>& res) const
	{
		double sumSq = 0.0;
		int count = 0;
		for (std::size_t i = 0; i < observations_.size(); ++i) {
			if (observations_[i].radar)
				continue;
			const int row = residualOffset(i);
			sumSq += res[row] * res[row];
			sumSq += res[row + 1] * res[row + 1];
			count += 2;
		}

		double rms = count ? std::sqrt(sumSq / static_cast<double>(count)) : 0.0;
		return Convertors::radiansToDegrees(rms) * 3600.0;
	}

	double OrbitFitModel::computeCostWeighted(const std::vector<double>& params) const
	{
		const auto res = computeResiduals(params);
		const std::size_t N = observations_.size();
		const bool useW = useWeights_;
		double sumSq = 0.0;
		for (std::size_t i = 0; i < N; ++i) {
			int row = residualOffset(i);
			if (observations_[i].radar) {
				const auto& radar = *observations_[i].radar;
				for (const auto& measurement : {radar.delaySeconds, radar.dopplerHz}) {
					if (!measurement)
						continue;
					const double value = res[row++] / measurement->sigma;
					sumSq += value * value;
				}
				continue;
			}
			const double wr = useW ? 1.0 / observations_[i].sigmaRa : 1.0;
			const double wd = useW ? 1.0 / observations_[i].sigmaDec : 1.0;
			const double rRA = res[row] * wr;
			const double rDec = res[row + 1] * wd;
			sumSq += rRA * rRA + rDec * rDec;
		}
		return sumSq;
	}

	Eigen::MatrixXd OrbitFitModel::computeJacobian(const std::vector<double>& x) const
	{
		const int n = static_cast<int>(x.size());
		const int m = residualCount_;
		Eigen::MatrixXd J(m, n);

		const int stateCols = std::min(6, n);

		for (int j = 0; j < stateCols; ++j) {
			double dv = (j < static_cast<int>(deltaVals_.size())) ? deltaVals_[j] : 1e-8;
			double h = std::max(std::abs(x[j]) * 1e-8, dv);
			std::vector<double> xp = x, xm = x;
			xp[j] += h;
			xm[j] -= h;
			const auto fp = computeResiduals(xp);
			const auto fm = computeResiduals(xm);
			for (int i = 0; i < m; ++i)
				J(i, j) = (fp[i] - fm[i]) / (2.0 * h);
		}

		for (int j = stateCols; j < n; ++j) {
			double dv = (j < static_cast<int>(deltaVals_.size())) ? deltaVals_[j] : 1e-8;
			double h = std::max(std::abs(x[j]) * 1e-8, dv);
			std::vector<double> xp = x, xm = x;
			xp[j] += h;
			xm[j] -= h;

			NonGravModel localNg = nbody_.nonGrav();
			NonGravProfile::updateFromParams(localNg, xp, 6);

			const auto fp = computeResiduals(xp, localNg);
			NonGravProfile::updateFromParams(localNg, xm, 6);

			const auto fm = computeResiduals(xm, localNg);
			for (int i = 0; i < m; ++i)
				J(i, j) = (fp[i] - fm[i]) / (2.0 * h);
		}
		return J;
	}

	Orbit_p OrbitFitModel::orbitFromStateVector(const std::vector<double>& params) const
	{
		if (params.size() < 6)
			throw std::runtime_error("Expected at least 6 params: Rx Ry Rz Vx Vy Vz");

		Vec3D R(params[0], params[1], params[2]);
		Vec3D V(params[3], params[4], params[5]);
		auto orbit = make_orbit_from_states_vectors(StateVectors(R, V), initialOrbit_->getEpoch());

		if (!orbit)
			throw std::runtime_error("Failed to build orbit from state vectors");

		return orbit;
	}

	std::vector<double> OrbitFitModel::stateVectorFromOrbit(const Orbit& orbit)
	{
		Vec3D R = orbit.getR();
		Vec3D V = orbit.getV();
		return {R.x, R.y, R.z, V.x, V.y, V.z};
	}

}  // namespace orbital_core::fit
