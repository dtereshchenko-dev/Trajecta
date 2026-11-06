#include "fit/OrbitFitUncertainty.h"
#include "fit/FitCancellation.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/Types.h"
#include "fit/APosterioriCovarianceScaling.h"
#include "logger/Logger.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>

namespace orbital_core::fit {

	FitCovariance OrbitFitUncertainty::computeFitCovariance(const std::vector<double>& params,
	                                                        const std::vector<double>& rawRes, bool log) const
	{
		std::vector<double> residualsForSigma = rawRes;
		std::vector<orbital_core::obs::AstrometricResidualTransform> transformsForSigma;
		Eigen::MatrixXd Jfull;
		if (problem_.snapshot().matches(params) &&
		    problem_.snapshot().rawJacobian.cols() == static_cast<Eigen::Index>(params.size()) &&
		    problem_.snapshot().rawResiduals.size() ==
		        static_cast<std::size_t>(problem_.snapshot().rawJacobian.rows())) {
			residualsForSigma = problem_.snapshot().rawResiduals;
			Jfull = problem_.snapshot().rawJacobian;
			transformsForSigma = problem_.snapshot().transforms;
		}

		else if (model_.usesNBody()) {
			auto sj = model_.computeResidualsAndStateJac(params);
			if (sj.stateJ.size() != 0) {
				residualsForSigma = std::move(sj.residuals);
				Jfull = std::move(sj.stateJ);
				transformsForSigma = std::move(sj.transforms);
			}
		}

		if (Jfull.size() == 0) {
			Jfull = model_.computeJacobian(params);
		}

		if (residualsForSigma.empty())
			residualsForSigma = model_.computeResiduals(params);

		Eigen::MatrixXd Jw = Jfull;
		problem_.weightJacobian(Jw, transformsForSigma);
		problem_.maskJacobian(Jw);

		problem_.weightResiduals(residualsForSigma, transformsForSigma);
		problem_.maskResiduals(residualsForSigma);
		const double chi2 = problem_.sumSquares(residualsForSigma);

		const int nParams = static_cast<int>(params.size());
		const int activeCount = problem_.activeResidualCount();
		const double dof = (activeCount > nParams) ? static_cast<double>(activeCount - nParams) : 1.0;

		const double sigma2 = orbital_core::fit::aPosterioriVarianceFactor(chi2, activeCount, nParams);

		if (log)
			LOG_INFO << "[Covariance] chi2=" << chi2 << "  n=" << activeCount << "  params=" << nParams
			         << "  dof=" << dof << "  chi2/dof=" << (dof > 0.0 ? chi2 / dof : 0.0)
			         << "  n/dof=" << (dof > 0.0 ? static_cast<double>(activeCount) / dof : 0.0)
			         << "  sigma2=" << sigma2;

		// Normalize parameter columns before decomposition: state and force
		// parameters can differ by many orders of magnitude in their units.
		Eigen::VectorXd columnScale(nParams);
		double smallestSingular = std::numeric_limits<double>::quiet_NaN();
		double rankThreshold = std::numeric_limits<double>::quiet_NaN();
		int zeroColumn = -1;
		bool usable = Jw.rows() >= nParams && Jw.cols() == nParams && Jw.allFinite();
		for (int i = 0; i < nParams && usable; ++i) {
			const double norm = Jw.col(i).norm();
			usable = norm > 0.0 && std::isfinite(norm);
			if (!usable)
				zeroColumn = i;
			if (usable)
				columnScale(i) = 1.0 / norm;
		}
		Eigen::MatrixXd formal = Eigen::MatrixXd::Constant(nParams, nParams, std::numeric_limits<double>::quiet_NaN());
		if (usable) {
			const Eigen::MatrixXd scaled = Jw * columnScale.asDiagonal();
			const Eigen::JacobiSVD<Eigen::MatrixXd> svd(scaled, Eigen::ComputeThinV);
			const auto& singular = svd.singularValues();
			const double threshold = std::numeric_limits<double>::epsilon() *
			                         static_cast<double>(std::max(scaled.rows(), scaled.cols())) * singular(0);
			smallestSingular = singular(nParams - 1);
			rankThreshold = threshold;
			usable = svd.info() == Eigen::Success && singular.allFinite() && singular(nParams - 1) > threshold;
			if (usable) {
				const Eigen::MatrixXd transform = columnScale.asDiagonal() * svd.matrixV();
				formal = transform * singular.array().square().inverse().matrix().asDiagonal() * transform.transpose();
				usable = formal.allFinite();
			}
		}
		if (!usable && log)
			LOG_WARN << "[Covariance] Jacobian is non-finite or rank deficient; covariance unavailable"
			         << " (column=" << zeroColumn << ", smallest normalized singular=" << smallestSingular
			         << ", threshold=" << rankThreshold << ")";

		FitCovariance out;
		out.formalFullCovariance = formal;
		out.fullCovariance = formal * sigma2;
		out.varianceFactor = sigma2;

		return out;
	}

	double OrbitFitUncertainty::timeSincePerihelionOnNominalBranch(const cometary::Vector6& sampleCot,
	                                                               const cometary::Vector6& nominalCot)
	{
		const double nan = std::numeric_limits<double>::quiet_NaN();
		const double nominalNu = std::remainder(nominalCot[5], Constants::TWOPI);
		const double sampleNu = std::remainder(sampleCot[5], Constants::TWOPI);
		const double localNu = nominalNu + Helpers::wrapPi(sampleNu - nominalNu);
		const double turns = std::round((localNu - sampleNu) / Constants::TWOPI);
		double time = cometary::timeSincePerihelion(sampleCot[0], sampleCot[1], sampleNu);
		if (!std::isfinite(time) || !std::isfinite(turns))
			return nan;
		if (turns != 0.0) {
			if (!(nominalCot[1] < 1.0) || !(sampleCot[1] < 1.0))
				return nan;
			const double delta = 1.0 - sampleCot[1];
			const double period = Constants::TWOPI * sampleCot[0] * std::sqrt(sampleCot[0]) /
			                      (Constants::GAUSS_K * delta * std::sqrt(delta));
			if (!std::isfinite(period))
				return nan;
			time += turns * period;
		}
		return std::isfinite(time) ? time : nan;
	}

	std::optional<OrbitFitUncertainty::RegularQualityJacobian> OrbitFitUncertainty::regularQualityJacobianCot(
	    const cometary::Vector6& nominalCot)
	{
		const double nominalTime = timeSincePerihelionOnNominalBranch(nominalCot, nominalCot);
		if (!nominalCot.allFinite() || !std::isfinite(nominalTime))
			return std::nullopt;

		const double q = nominalCot[0];
		const double e = nominalCot[1];
		const double nu = nominalCot[5];
		const double cosNu = std::cos(nu);
		const double denominator = 1.0 + e * cosNu;
		const double relativeStep = std::cbrt(std::numeric_limits<double>::epsilon());
		double eccentricityScale = std::max(1.0, e);
		if (cosNu < 0.0)
			eccentricityScale = std::min(eccentricityScale, denominator / -cosNu);
		RegularQualityJacobian out = RegularQualityJacobian::Zero();
		out(0, 0) = 1.0;
		out(1, 0) = -(1.0 - e) / (q * q);
		out(1, 1) = -1.0 / q;
		out(2, 0) = -1.5 * nominalTime / q;
		out(2, 5) =
		    -q * std::sqrt(q) * (1.0 + e) * std::sqrt(1.0 + e) / (Constants::GAUSS_K * denominator * denominator);

		const double h = relativeStep * eccentricityScale;
		auto plus = nominalCot, minus = nominalCot;
		plus[1] += h;
		minus[1] -= h;
		const double width = plus[1] - minus[1];
		if (!(width > 0.0))
			return std::nullopt;
		const double timePlus = timeSincePerihelionOnNominalBranch(plus, nominalCot);
		const double timeMinus = timeSincePerihelionOnNominalBranch(minus, nominalCot);
		if (!std::isfinite(timePlus) || !std::isfinite(timeMinus))
			return std::nullopt;
		out(2, 1) = -(timePlus - timeMinus) / width;
		if (!out.allFinite())
			return std::nullopt;
		return out;
	}

	OrbitFitUncertainty::OutputTransport OrbitFitUncertainty::outputTransport(const std::vector<double>& params,
	                                                                          double fitEpochTdb, double outputEpochTdb,
	                                                                          std::stop_token cancel) const
	{
		checkCancellation(cancel);
		const Eigen::Index nParams = static_cast<Eigen::Index>(params.size());
		if (nParams < 6)
			throw std::invalid_argument("Output covariance requires a six-dimensional state");
		OutputTransport out;
		out.state.assign(params.begin(), params.begin() + 6);
		out.jacobian.setZero(6, nParams);
		out.jacobian.leftCols<6>().setIdentity();
		if (outputEpochTdb == fitEpochTdb)
			return out;

		NonGravModel ng = model_.propagator().nonGrav();
		NonGravProfile::updateFromParams(ng, params, 6);
		const Vec3D r(params[0], params[1], params[2]);
		const Vec3D v = Vec3D(params[3], params[4], params[5]) * Constants::GAUSS_K;
		Vec3D outputR, outputV;
		if (nParams > 6) {
			if (nParams - 6 != ng.numFitParams())
				throw std::invalid_argument("Output covariance non-gravitational parameter mismatch");
			const auto states = model_.propagator().propagateAllWithStmAndSensMjd(r, v, fitEpochTdb - 2400000.5,
			                                                                      {outputEpochTdb - 2400000.5}, ng);
			checkCancellation(cancel);
			if (states.size() != 1)
				throw std::runtime_error("Output covariance propagation returned no state");
			const auto& state = states.front();
			outputR = state.R;
			outputV = state.V;
			out.jacobian.leftCols<6>() = state.Phi;
			for (Eigen::Index col = 6; col < nParams; ++col)
				out.jacobian.col(col) = state.SensA.col(ng.sensitivityColumnForFitParameter(static_cast<int>(col - 6)));
		} else {
			const auto states = model_.propagator().propagateAllWithStmMjd(r, v, fitEpochTdb - 2400000.5,
			                                                               {outputEpochTdb - 2400000.5}, ng);
			checkCancellation(cancel);
			if (states.size() != 1)
				throw std::runtime_error("Output covariance propagation returned no state");
			outputR = states.front().R;
			outputV = states.front().V;
			out.jacobian = states.front().Phi;
		}
		out.jacobian.middleCols<3>(3) *= Constants::GAUSS_K;
		out.jacobian.bottomRows<3>() /= Constants::GAUSS_K;
		outputV /= Constants::GAUSS_K;
		out.state = {outputR.x, outputR.y, outputR.z, outputV.x, outputV.y, outputV.z};
		if (!out.jacobian.allFinite() ||
		    !std::all_of(out.state.begin(), out.state.end(), [](double value) { return std::isfinite(value); }))
			throw std::runtime_error("Output covariance propagation returned a non-finite state");
		return out;
	}

	OrbitFitUncertainty::CometaryDerivatives OrbitFitUncertainty::computeCometaryDerivatives(
	    const std::vector<double>& params, const Orbit& fitOrbit, const Orbit& displayOrbit,
	    std::stop_token cancel) const
	{
		checkCancellation(cancel);
		CometaryDerivatives out;
		const int nParams = static_cast<int>(params.size());
		const double nan = std::numeric_limits<double>::quiet_NaN();
		out.cotJ.setConstant(6, nParams, nan);
		out.regularQualityJ.setConstant(3, nParams, nan);
		out.epochJdTdb = displayOrbit.getEpoch();
		if (nParams < 6)
			return out;

		cometary::Vector6 nominalCot;
		nominalCot << displayOrbit.getPerigee(), displayOrbit.getExcentricity(), displayOrbit.getInclination(),
		    displayOrbit.getLongitudeOfAscendingNode(), displayOrbit.getArgumentOfPericentre(),
		    displayOrbit.getTrueAnomaly();
		const auto localCotJacobian = cometary::cartesianToElementsJacobian(nominalCot);
		const auto qualityJacobian = regularQualityJacobianCot(nominalCot);
		if (!localCotJacobian || !qualityJacobian)
			return out;

		if (out.epochJdTdb == fitOrbit.getEpoch()) {
			out.cotJ.setZero();
			out.cotJ.leftCols<6>() = *localCotJacobian;
			out.cotJ.middleCols<3>(3) *= Constants::GAUSS_K;
			out.regularQualityJ = *qualityJacobian * out.cotJ;
			out.valid = out.cotJ.allFinite() && out.regularQualityJ.allFinite();
			return out;
		}

		const auto transport = outputTransport(params, fitOrbit.getEpoch(), out.epochJdTdb, cancel);
		cometary::Matrix6 localCanonicalJacobian = *localCotJacobian;
		localCanonicalJacobian.rightCols<3>() *= Constants::GAUSS_K;
		out.cotJ = localCanonicalJacobian * transport.jacobian;
		out.regularQualityJ = *qualityJacobian * out.cotJ;
		out.valid = out.cotJ.allFinite() && out.regularQualityJ.allFinite();
		return out;
	}

	void OrbitFitUncertainty::populateCometaryUncertainty(Uncertainty& uncertainty,
	                                                      const CometaryDerivatives& derivatives,
	                                                      const FitCovariance& covariance, const Orbit& displayOrbit)
	{
		const double nan = std::numeric_limits<double>::quiet_NaN();
		uncertainty.elementCovariance.setConstant(nan);
		uncertainty.hasElementCovariance = false;
		uncertainty.cometaryElementCovariance.reset();
		auto& sig = uncertainty.derivedSigmas;
		sig.sigma_a = sig.sigma_e = sig.sigma_i = sig.sigma_Om = sig.sigma_w = sig.sigma_nu = nan;
		sig.sigma_q = sig.sigma_inv_a = sig.sigma_Tp = sig.sigma_M = sig.sigma_Q = nan;
		sig.sigma_n = sig.sigma_P = sig.sigma_Py = nan;

		const Eigen::Index nParams = derivatives.cotJ.cols();
		if (!derivatives.valid || nParams < 6 || derivatives.regularQualityJ.cols() != nParams ||
		    covariance.fullCovariance.rows() != nParams || covariance.fullCovariance.cols() != nParams ||
		    covariance.formalFullCovariance.rows() != nParams || covariance.formalFullCovariance.cols() != nParams ||
		    !std::isfinite(derivatives.epochJdTdb) || derivatives.epochJdTdb != displayOrbit.getEpoch())
			return;

		const Eigen::Matrix<double, 6, 6> cotCovariance =
		    derivatives.cotJ * covariance.fullCovariance * derivatives.cotJ.transpose();
		const Eigen::Matrix<double, 6, 6> formalCotCovariance =
		    derivatives.cotJ * covariance.formalFullCovariance * derivatives.cotJ.transpose();
		const Eigen::Matrix3d qualityCovariance =
		    derivatives.regularQualityJ * covariance.fullCovariance * derivatives.regularQualityJ.transpose();
		if (!cotCovariance.allFinite() || !formalCotCovariance.allFinite() || !qualityCovariance.allFinite())
			return;

		CometaryElementCovariance cot;
		cot.covariance = cotCovariance;
		cot.formalCovariance = formalCotCovariance;
		cot.epochJdTdb = derivatives.epochJdTdb;
		uncertainty.cometaryElementCovariance = cot;

		auto sigma = [nan](double variance) {
			return (variance >= 0.0 && std::isfinite(variance)) ? std::sqrt(variance) : nan;
		};
		sig.sigma_q = sigma(qualityCovariance(0, 0));
		sig.sigma_inv_a = sigma(qualityCovariance(1, 1));
		sig.sigma_Tp = sigma(qualityCovariance(2, 2));
		sig.sigma_e = sigma(cotCovariance(1, 1));
		sig.sigma_i = sigma(cotCovariance(2, 2));
		sig.sigma_Om = sigma(cotCovariance(3, 3));
		sig.sigma_w = sigma(cotCovariance(4, 4));
		sig.sigma_nu = sigma(cotCovariance(5, 5));

		const double q = displayOrbit.getPerigee();
		const double e = displayOrbit.getExcentricity();
		const double invA = (1.0 - e) / q;
		if (!(q > 0.0) || !std::isfinite(invA) || invA == 0.0 || displayOrbit.getType() == OrbitType::Parabolic)
			return;

		Eigen::Matrix<double, 6, Eigen::Dynamic> conventionalJ = derivatives.cotJ;
		conventionalJ.row(0) = -derivatives.regularQualityJ.row(1) / (invA * invA);
		const Eigen::Matrix<double, 6, 6> conventionalCovariance =
		    conventionalJ * covariance.fullCovariance * conventionalJ.transpose();
		if (conventionalCovariance.allFinite()) {
			uncertainty.elementCovariance = conventionalCovariance;
			uncertainty.hasElementCovariance = true;
			sig.sigma_a = sigma(conventionalCovariance(0, 0));
		}

		const double absInvA = std::abs(invA);
		const double n = Constants::GAUSS_K * absInvA * std::sqrt(absInvA);
		const double dn = 1.5 * Constants::GAUSS_K * std::sqrt(absInvA) * std::copysign(1.0, invA);
		const double nu = displayOrbit.getTrueAnomaly();
		const double denominator = 1.0 + e * std::cos(nu);
		const double root = std::sqrt(std::abs(1.0 - e)) * std::sqrt(1.0 + e);
		Eigen::Matrix<double, 1, 6> meanAnomalyGradient = Eigen::Matrix<double, 1, 6>::Zero();
		meanAnomalyGradient[1] =
		    std::copysign(root, e - 1.0) * std::sin(nu) * (2.0 + e * std::cos(nu)) / (denominator * denominator);
		meanAnomalyGradient[5] = std::abs(1.0 - e) * (1.0 + e) * root / (denominator * denominator);
		sig.sigma_M = sigma((meanAnomalyGradient * cotCovariance * meanAnomalyGradient.transpose())(0, 0));
		sig.sigma_n = std::abs(dn) * sig.sigma_inv_a;
		if (e < 1.0) {
			const Eigen::Vector3d aphelionGradient(-1.0, -2.0 / (invA * invA), 0.0);
			sig.sigma_Q = sigma((aphelionGradient.transpose() * qualityCovariance * aphelionGradient)(0, 0));
			sig.sigma_P = Constants::TWOPI * sig.sigma_n / (n * n);
			sig.sigma_Py = sig.sigma_P / 365.25;
		}
	}

	OrbitFitUncertainty::ElementDerivatives OrbitFitUncertainty::computeElementDerivatives(
	    const std::vector<double>& params, const Orbit& orbit, std::stop_token cancel) const
	{
		ElementDerivatives out;
		const int nParams = static_cast<int>(params.size());
		out.elemJ.resize(6, nParams);
		out.qualityElemJ.resize(9, nParams);
		out.elemJ.setZero();
		out.qualityElemJ.setZero();

		const double fitEpochTdb = orbit.getEpoch();
		const double targetEpochTdb =
		    outputEpochTdb_ > 0.0 && std::abs(outputEpochTdb_ - fitEpochTdb) > 1e-6 ? outputEpochTdb_ : fitEpochTdb;
		const auto transport = outputTransport(params, fitEpochTdb, targetEpochTdb, cancel);
		Eigen::Matrix<double, 6, 6> localElementJ;
		Eigen::Matrix<double, 9, 6> localQualityJ;
		auto convert = [&](const std::vector<double>& state) -> std::pair<std::vector<double>, std::vector<double>> {
			checkCancellation(cancel);
			const auto orbAt = targetEpochTdb == fitEpochTdb
			                       ? model_.orbitFromStateVector(state)
			                       : make_orbit_from_states_vectors(StateVectors(Vec3D(state[0], state[1], state[2]),
			                                                                     Vec3D(state[3], state[4], state[5])),
			                                                        targetEpochTdb);
			if (!orbAt)
				throw std::runtime_error("Output covariance element conversion failed");
			return {orbitalElements(*orbAt), qualityCovarianceElements(*orbAt)};
		};

		auto diffElement = [](int row, double plus, double minus) noexcept {
			if (row == 3 || row == 4 || row == 5)
				return Helpers::wrapPi(plus - minus);

			return plus - minus;
		};

		auto diffQualityElement = [](int row, double plus, double minus) noexcept {
			if (row == 6 || row == 7 || row == 8)
				return Helpers::wrapPi(plus - minus);

			return plus - minus;
		};

		for (int j = 0; j < 6; ++j) {
			const double h = std::max(1.0e-7, 1.0e-5 * std::abs(transport.state[j]));

			auto xp = transport.state, xm = transport.state;
			xp[j] += h;
			xm[j] -= h;

			auto [elemP, qualityElemP] = convert(xp);
			auto [elemM, qualityElemM] = convert(xm);

			const double inv2h = 1.0 / (2.0 * h);

			for (int i = 0; i < 6; ++i)
				localElementJ(i, j) = diffElement(i, elemP[i], elemM[i]) * inv2h;

			for (int i = 0; i < 9; ++i)
				localQualityJ(i, j) = diffQualityElement(i, qualityElemP[i], qualityElemM[i]) * inv2h;
		}
		out.elemJ = localElementJ * transport.jacobian;
		out.qualityElemJ = localQualityJ * transport.jacobian;

		return out;
	}

	void OrbitFitUncertainty::populateOrbitElementSigmas(Orbit& orbit,
	                                                     const orbital_core::uncertainty::DerivedSigmas& sig) noexcept
	{
		orbit.setSigmaSemiMajorAxis(sig.sigma_a);
		orbit.setSigmaExcentricity(sig.sigma_e);
		orbit.setSigmaInclination(sig.sigma_i);
		orbit.setSigmaLongitudeOfAscendingNode(sig.sigma_Om);
		orbit.setSigmaArgumentOfPericentre(sig.sigma_w);
		orbit.setSigmaTrueAnomaly(sig.sigma_nu);
		orbit.setSigmaPerigee(sig.sigma_q);
	}

	orbital_core::uncertainty::DerivedSigmas OrbitFitUncertainty::computeDerivedSigmasFromQualityCovariance(
	    const Orbit& orbit, const orbital_core::uncertainty::ElementCov6x6& covOrb,
	    const Eigen::MatrixXd& qualityElementCov)
	{
		auto sig = orbital_core::uncertainty::DerivedSigmasComputer::compute(orbit, covOrb);
		if (qualityElementCov.rows() < 9 || qualityElementCov.cols() < 9)
			return sig;

		auto sigmaAt = [&](int idx) noexcept {
			const double v = qualityElementCov(idx, idx);
			return (v > 0.0 && std::isfinite(v)) ? std::sqrt(v) : 0.0;
		};

		const double a = orbit.getSemiMajorAxis();
		const double absA = std::abs(a);
		const bool elliptical = (orbit.getType() == OrbitType::Elliptical);

		sig.sigma_Tp = sigmaAt(0);
		sig.sigma_e = sigmaAt(1);
		sig.sigma_q = sigmaAt(2);
		sig.sigma_Q = elliptical ? sigmaAt(3) : std::numeric_limits<double>::quiet_NaN();
		sig.sigma_inv_a = sigmaAt(4);
		sig.sigma_i = sigmaAt(5);
		sig.sigma_M = sigmaAt(6);
		sig.sigma_w = sigmaAt(7);
		sig.sigma_Om = sigmaAt(8);

		if (absA > 0.0) {
			sig.sigma_a = sig.sigma_inv_a * a * a;

			const double n_radPerDay = std::sqrt(Constants::GM_SUN / (absA * absA * absA));
			sig.sigma_n = 1.5 * n_radPerDay * sig.sigma_a / absA;

			if (elliptical) {
				const double perDays = 2.0 * Constants::PI * std::sqrt(absA * absA * absA / Constants::GM_SUN);
				sig.sigma_P = 1.5 * perDays * sig.sigma_a / absA;
				sig.sigma_Py = sig.sigma_P / 365.25;
			}

			else {
				sig.sigma_P = std::numeric_limits<double>::quiet_NaN();
				sig.sigma_Py = std::numeric_limits<double>::quiet_NaN();
			}
		}

		return sig;
	}

	OrbitFitUncertainty::OrbitQualityMetrics OrbitFitUncertainty::computeOrbitQualityMetrics(
	    const Orbit& orbit, const orbital_core::uncertainty::DerivedSigmas& sig) noexcept
	{
		OrbitQualityMetrics out;
		if (orbit.getType() != OrbitType::Elliptical)
			return out;

		out.elliptical = true;

		constexpr double kDaysYr = 365.25;
		constexpr double kUconst = 1.49;
		constexpr double kRunoffCoeff = 3.0 * 3600.0 * (180.0 / Constants::PI) * Constants::GAUSS_K;

		const double a = orbit.getSemiMajorAxis();
		const double e = orbit.getExcentricity();
		if (!(a > 0.0) || !std::isfinite(e))
			return out;

		const double periodDays = Constants::TWOPI * std::sqrt(a * a * a / Constants::GM_SUN);
		const double periodYears = periodDays / kDaysYr;
		if (!(periodYears > 0.0) || !std::isfinite(periodYears))
			return out;

		out.termTp = sig.sigma_Tp * e;
		out.termInvA = 15.0 * kDaysYr * sig.sigma_inv_a * a;
		const double termSum = out.termTp + out.termInvA;
		out.runoff = kRunoffCoeff * termSum / periodYears;

		if (!(out.runoff > 0.0) || !std::isfinite(out.runoff))
			return out;

		out.u = std::log(out.runoff) / kUconst + 1.0;
		out.valid = std::isfinite(out.u);
		return out;
	}

	void OrbitFitUncertainty::logOrbitQualityU(const Orbit& orbit, const orbital_core::uncertainty::DerivedSigmas& sig)
	{
		const auto metrics = computeOrbitQualityMetrics(orbit, sig);
		if (!metrics.elliptical) {
			LOG_INFO << "Orbit quality: U not defined for non-elliptical orbits "
			            "(hyperbolic / parabolic — runoff over period is undefined)";
			return;
		}

		if (!metrics.valid) {
			LOG_INFO << "Orbit quality: U not defined (non-positive runoff)";
			return;
		}

		const double termSum = metrics.termTp + metrics.termInvA;
		const double pct_Tp = (termSum > 0.0) ? (100.0 * metrics.termTp / termSum) : 0.0;
		const double pct_invA = (termSum > 0.0) ? (100.0 * metrics.termInvA / termSum) : 0.0;

		LOG_INFO << "Orbit quality runoff decomposition:";
		LOG_INFO << std::scientific << std::setprecision(3) << "  σ_Tp     = " << sig.sigma_Tp << " d";

		LOG_INFO << std::scientific << std::setprecision(3) << "  σ(1/a)   = " << sig.sigma_inv_a << " 1/AU";

		LOG_INFO << std::scientific << std::setprecision(3) << "  term σ_Tp·e          = " << metrics.termTp
		         << std::fixed << std::setprecision(1) << "  (" << pct_Tp << "%)";

		LOG_INFO << std::scientific << std::setprecision(3) << "  term 15·yr·σ(1/a)·a  = " << metrics.termInvA
		         << std::fixed << std::setprecision(1) << "  (" << pct_invA << "%)";

		LOG_INFO << std::scientific << std::setprecision(3) << "Orbit quality: runoff = " << metrics.runoff
		         << " arcsec → U = " << std::fixed << std::setprecision(2) << metrics.u;
	}

	std::vector<double> OrbitFitUncertainty::orbitalElements(const Orbit& orbit)
	{
		return {orbit.getSemiMajorAxis(),        orbit.getExcentricity(),
		        orbit.getInclination(),          orbit.getLongitudeOfAscendingNode(),
		        orbit.getArgumentOfPericentre(), orbit.getTrueAnomaly()};
	}

	std::vector<double> OrbitFitUncertainty::qualityCovarianceElements(const Orbit& orbit)
	{
		return {orbit.getTimeOfPerihelion(),
		        orbit.getExcentricity(),
		        orbit.getPerigee(),
		        orbit.getApogee(),
		        orbit.getReciprocalSemiMajorAxis(),
		        orbit.getInclination(),
		        orbit.getMeanAnomaly(),
		        orbit.getArgumentOfPericentre(),
		        orbit.getLongitudeOfAscendingNode()};
	}

	double OrbitFitUncertainty::orbitQualityU(const Orbit& orbit, const uncertainty::DerivedSigmas& sigmas) noexcept
	{
		return computeOrbitQualityMetrics(orbit, sigmas).u;
	}

}  // namespace orbital_core::fit
