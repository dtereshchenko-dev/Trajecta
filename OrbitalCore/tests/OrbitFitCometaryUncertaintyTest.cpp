#include <gtest/gtest.h>

#include "fit/OrbitFitUncertainty.h"
#include "orbital_core/CometaryCoordinates.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Types.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace {
	using orbital_core::fit::FitCovariance;
	using orbital_core::fit::OrbitFitModel;
	using orbital_core::fit::OrbitFitProblem;
	using orbital_core::fit::OrbitFitUncertainty;
	namespace cometary = orbital_core::cometary;

	constexpr double fitEpoch = 2460000.5;

	Orbit_p cometOrbit(double e, double nu = -0.7, double epoch = fitEpoch, double inclination = 0.41)
	{
		cometary::Vector6 cot;
		cot << 2.3, e, inclination, 0.82, 0.37, nu;
		const auto cartesian = cometary::toCartesian(cot);
		if (!cartesian)
			return nullptr;
		auto orbit = std::make_shared<Orbit>();
		orbit->setEpoch(epoch);
		orbit->setType(e == 1.0 ? OrbitType::Parabolic : (e < 1.0 ? OrbitType::Elliptical : OrbitType::Hyperbolic));
		orbit->setR(Vec3D((*cartesian)[0], (*cartesian)[1], (*cartesian)[2]));
		orbit->setV(Vec3D((*cartesian)[3], (*cartesian)[4], (*cartesian)[5]) / Constants::GAUSS_K);
		orbit->setPerigee(cot[0]);
		orbit->setExcentricity(e);
		orbit->setInclination(cot[2]);
		orbit->setLongitudeOfAscendingNode(cot[3]);
		orbit->setArgumentOfPericentre(cot[4]);
		orbit->setTrueAnomaly(nu);
		orbit->setSemiMajorAxis(e == 1.0 ? std::numeric_limits<double>::infinity() : cot[0] / (1.0 - e));
		orbit->setReciprocalSemiMajorAxis((1.0 - e) / cot[0]);
		orbit->setTimeOfPerihelion(epoch - cometary::timeSincePerihelion(cot[0], e, nu));
		return orbit;
	}

	FitCovariance isotropicCovariance(int size)
	{
		FitCovariance result;
		result.formalFullCovariance = 1e-12 * Eigen::MatrixXd::Identity(size, size);
		result.fullCovariance = 1.04 * result.formalFullCovariance;
		result.varianceFactor = 1.04;
		return result;
	}
}  // namespace

TEST(OrbitFitModelTest, EncounterControlCoversGravityAndNonGravButNotKepler)
{
	struct RestoreEphemeris {
		IEphemerisSource* previous{getDefaultEphemerisSource()};

		~RestoreEphemeris() { setDefaultEphemerisSource(previous); }
	} restore;

	DE440Source ephemeris;
	setDefaultEphemerisSource(&ephemeris);
	for (const bool gravity : {false, true}) {
		for (const bool nonGrav : {false, true}) {
			SCOPED_TRACE(gravity);
			SCOPED_TRACE(nonGrav);
			OrbitOptimizerConfig config;
			config.nbody.enabled = gravity;
			config.nbody.planetMask = 4;
			config.nbody.stepSize = 2.0;
			config.nonGrav.enabled = nonGrav;
			config.nonGrav.model = orbital_core::NonGravModelKind::InverseSquare;
			OrbitFitModel model(cometOrbit(0.999), {}, config);
			EXPECT_EQ(model.propagator().encounterStepPointCount(), gravity || nonGrav ? 100 : 0);
			EXPECT_EQ(model.propagator().planetMask(), gravity ? 4 : 0);
			EXPECT_DOUBLE_EQ(model.propagator().stepSize(), 2.0);
			EXPECT_EQ(model.propagator().nonGrav().enabled(), nonGrav);
		}
	}
}

TEST(OrbitFitCometaryUncertaintyTest, FitEpochJacobianStaysRegularAcrossParabola)
{
	for (const double e : {1.0 - 1e-10, 1.0, 1.0 + 1e-10}) {
		SCOPED_TRACE(e);
		const auto orbit = cometOrbit(e);
		ASSERT_TRUE(orbit);
		OrbitOptimizerConfig config;
		OrbitFitModel model(orbit, {}, config);
		OrbitFitProblem problem(model);
		OrbitFitUncertainty calculator(model, problem, fitEpoch);
		const auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
		const auto derivatives = calculator.computeCometaryDerivatives(params, *orbit, *orbit, {});
		ASSERT_TRUE(derivatives.valid);
		EXPECT_TRUE(derivatives.cotJ.allFinite());
		EXPECT_TRUE(derivatives.regularQualityJ.allFinite());
		EXPECT_DOUBLE_EQ(derivatives.epochJdTdb, fitEpoch);

		cometary::Vector6 cot;
		cot << orbit->getPerigee(), e, orbit->getInclination(), orbit->getLongitudeOfAscendingNode(),
		    orbit->getArgumentOfPericentre(), orbit->getTrueAnomaly();
		const auto analytic = cometary::cartesianToElementsJacobian(cot);
		ASSERT_TRUE(analytic);
		cometary::Matrix6 expected = *analytic;
		expected.rightCols<3>() *= Constants::GAUSS_K;
		for (int row = 0; row < 6; ++row)
			EXPECT_LT((derivatives.cotJ.row(row) - expected.row(row)).norm(), 1e-13 * expected.row(row).norm());

		orbital_core::Uncertainty uncertainty;
		calculator.populateCometaryUncertainty(uncertainty, derivatives, isotropicCovariance(6), *orbit);
		ASSERT_TRUE(uncertainty.cometaryElementCovariance);
		EXPECT_TRUE(uncertainty.cometaryElementCovariance->covariance.allFinite());
		EXPECT_TRUE(uncertainty.cometaryElementCovariance->formalCovariance.allFinite());
		EXPECT_FALSE(uncertainty.empty());
		EXPECT_GT(uncertainty.derivedSigmas.sigma_q, 0.0);
		EXPECT_GT(uncertainty.derivedSigmas.sigma_inv_a, 0.0);
		EXPECT_GT(uncertainty.derivedSigmas.sigma_Tp, 0.0);
		if (e == 1.0) {
			EXPECT_FALSE(uncertainty.hasElementCovariance);
			EXPECT_TRUE(uncertainty.elementCovariance.array().isNaN().all());
			EXPECT_TRUE(std::isnan(uncertainty.derivedSigmas.sigma_a));
			EXPECT_TRUE(std::isnan(uncertainty.derivedSigmas.sigma_M));
			EXPECT_TRUE(std::isnan(uncertainty.derivedSigmas.sigma_P));
			EXPECT_TRUE(std::isnan(uncertainty.derivedSigmas.sigma_Q));
			EXPECT_TRUE(std::isnan(uncertainty.derivedSigmas.sigma_n));
		} else {
			EXPECT_TRUE(uncertainty.hasElementCovariance);
			EXPECT_TRUE(uncertainty.elementCovariance.allFinite());
			EXPECT_TRUE(std::isfinite(uncertainty.derivedSigmas.sigma_a));
		}
	}
}

TEST(OrbitFitCometaryUncertaintyTest, RegularSigmasAreContinuousAcrossParabola)
{
	std::array<Eigen::Vector3d, 3> sigmas;
	for (int index = 0; index < 3; ++index) {
		const auto orbit = cometOrbit(1.0 + (index - 1) * 1e-10);
		ASSERT_TRUE(orbit);
		OrbitFitModel model(orbit, {}, {});
		OrbitFitProblem problem(model);
		OrbitFitUncertainty calculator(model, problem, fitEpoch);
		const auto derivatives =
		    calculator.computeCometaryDerivatives(OrbitFitModel::stateVectorFromOrbit(*orbit), *orbit, *orbit, {});
		orbital_core::Uncertainty uncertainty;
		calculator.populateCometaryUncertainty(uncertainty, derivatives, isotropicCovariance(6), *orbit);
		ASSERT_TRUE(uncertainty.cometaryElementCovariance);
		sigmas[index] << uncertainty.derivedSigmas.sigma_q, uncertainty.derivedSigmas.sigma_inv_a,
		    uncertainty.derivedSigmas.sigma_Tp;
	}
	for (int index : {0, 2})
		for (int row = 0; row < 3; ++row)
			EXPECT_NEAR(sigmas[index][row], sigmas[1][row], 1e-5 * sigmas[1][row]);
}

TEST(OrbitFitCometaryUncertaintyTest, ProjectionRetainsNonGravCrossCovariance)
{
	const auto orbit = cometOrbit(1.0);
	ASSERT_TRUE(orbit);
	OrbitFitUncertainty::CometaryDerivatives derivatives;
	derivatives.cotJ.setZero(6, 7);
	derivatives.cotJ.leftCols<6>().setIdentity();
	derivatives.cotJ(0, 6) = 2.0;
	derivatives.cotJ(1, 6) = 0.4;
	derivatives.regularQualityJ.setZero(3, 7);
	derivatives.regularQualityJ.row(0) = derivatives.cotJ.row(0);
	derivatives.regularQualityJ.row(1) = -derivatives.cotJ.row(1) / orbit->getPerigee();
	derivatives.regularQualityJ.row(2) = 9.0 * derivatives.cotJ.row(5);
	derivatives.epochJdTdb = fitEpoch;
	derivatives.valid = true;
	FitCovariance covariance;
	covariance.formalFullCovariance = Eigen::MatrixXd::Identity(7, 7);
	covariance.formalFullCovariance(0, 6) = covariance.formalFullCovariance(6, 0) = 0.25;
	covariance.fullCovariance = 4.0 * covariance.formalFullCovariance;
	orbital_core::Uncertainty uncertainty;
	uncertainty.cartesianStateCovariance.setConstant(3.0);
	uncertainty.cartesianStateCovarianceFormal.setConstant(5.0);
	uncertainty.fitParameterCovariance = covariance.fullCovariance;
	uncertainty.fitParameterCovarianceFormal = covariance.formalFullCovariance;
	OrbitFitUncertainty::populateCometaryUncertainty(uncertainty, derivatives, covariance, *orbit);
	ASSERT_TRUE(uncertainty.cometaryElementCovariance);
	EXPECT_DOUBLE_EQ(uncertainty.cometaryElementCovariance->formalCovariance(0, 0), 6.0);
	EXPECT_DOUBLE_EQ(uncertainty.cometaryElementCovariance->covariance(0, 0), 24.0);
	EXPECT_DOUBLE_EQ(uncertainty.derivedSigmas.sigma_q, std::sqrt(24.0));
	EXPECT_DOUBLE_EQ(uncertainty.cometaryElementCovariance->epochJdTdb, fitEpoch);
	EXPECT_TRUE((uncertainty.cartesianStateCovariance.array() == 3.0).all());
	EXPECT_TRUE((uncertainty.cartesianStateCovarianceFormal.array() == 5.0).all());
	EXPECT_TRUE(uncertainty.fitParameterCovariance.isApprox(covariance.fullCovariance, 0.0));
	EXPECT_TRUE(uncertainty.fitParameterCovarianceFormal.isApprox(covariance.formalFullCovariance, 0.0));
}

TEST(OrbitFitCometaryUncertaintyTest, PerihelionTimeDerivativeRemainsSignedAcrossAngleWrap)
{
	const auto orbit = cometOrbit(1.0, 0.0);
	ASSERT_TRUE(orbit);
	OrbitFitModel model(orbit, {}, {});
	OrbitFitProblem problem(model);
	OrbitFitUncertainty calculator(model, problem, fitEpoch);
	const auto derivatives =
	    calculator.computeCometaryDerivatives(OrbitFitModel::stateVectorFromOrbit(*orbit), *orbit, *orbit, {});
	ASSERT_TRUE(derivatives.valid);
	const double q = orbit->getPerigee();
	const double timeDerivative = std::sqrt(2.0) * q * std::sqrt(q) / (2.0 * Constants::GAUSS_K);
	const Eigen::RowVectorXd expected = -timeDerivative * derivatives.cotJ.row(5);
	EXPECT_LT((derivatives.regularQualityJ.row(2) - expected).norm(), 2e-6 * expected.norm());
}

TEST(OrbitFitCometaryUncertaintyTest, DisplayTransportRetainsNonGravColumns)
{
	const auto orbit = cometOrbit(1.0);
	ASSERT_TRUE(orbit);
	OrbitOptimizerConfig config;
	config.nonGrav.enabled = true;
	config.nonGrav.model = orbital_core::NonGravModelKind::InverseSquare;
	OrbitFitModel model(orbit, {}, config);
	OrbitFitProblem problem(model);
	constexpr double displayEpoch = fitEpoch + 5.0;
	OrbitFitUncertainty calculator(model, problem, displayEpoch);
	const auto states = model.propagator().propagateAllWithVelocity(orbit->getR(), orbit->getV() * Constants::GAUSS_K,
	                                                                fitEpoch, {displayEpoch});
	ASSERT_EQ(states.size(), 1u);
	auto display =
	    make_orbit_from_states_vectors(StateVectors(states[0].R, states[0].V / Constants::GAUSS_K), displayEpoch);
	ASSERT_TRUE(display);
	auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
	params.resize(6 + model.propagator().nonGrav().numFitParams(), 0.0);
	ASSERT_GT(params.size(), 6u);
	const auto derivatives = calculator.computeCometaryDerivatives(params, *orbit, *display, {});
	ASSERT_TRUE(derivatives.valid);
	EXPECT_EQ(derivatives.cotJ.cols(), static_cast<Eigen::Index>(params.size()));
	EXPECT_GT(derivatives.cotJ.rightCols(params.size() - 6).norm(), 0.0);
	orbital_core::Uncertainty uncertainty;
	calculator.populateCometaryUncertainty(uncertainty, derivatives,
	                                       isotropicCovariance(static_cast<int>(params.size())), *display);
	ASSERT_TRUE(uncertainty.cometaryElementCovariance);
	EXPECT_DOUBLE_EQ(uncertainty.cometaryElementCovariance->epochJdTdb, displayEpoch);
}

TEST(OrbitFitCometaryUncertaintyTest, BorisovOutputTransportIsStableAndResolvesDelayDerivative)
{
	const double epoch = 2400000.5 + 58929.741139240563;
	const double outputEpoch = 2400000.5 + 59122.0;
	const Vec3D r(-1.8708360327800977, -1.1481752691438272, -2.1095460609000192);
	const Vec3D v(-4.6254083702732560e-4, -1.9822243926802881e-2, -1.2197657542355590e-2);
	const auto orbit = make_orbit_from_states_vectors(StateVectors(r, v / Constants::GAUSS_K), epoch);
	ASSERT_TRUE(orbit);
	OrbitOptimizerConfig config;
	config.nonGrav.enabled = true;
	config.nonGrav.model = orbital_core::NonGravModelKind::MarsdenH2O;
	OrbitFitModel model(orbit, {}, config);
	auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
	model.nonGravProfile().extendParams(params, model.differenceSteps());
	params[6] = 7.528364e-8;
	params[7] = 4.183832e-7;
	params[8] = -4.352873e-9;
	params[9] = -92.847801;
	model.nonGravProfile().updateFromParams(params, 6);
	OrbitFitProblem problem(model);
	OrbitFitUncertainty calculator(model, problem, outputEpoch);
	const auto sample = [&](const std::vector<double>& state) {
		NonGravModel ng = model.propagator().nonGrav();
		NonGravProfile::updateFromParams(ng, state, 6);
		const auto propagated = model.propagator().propagateAllWithVelocity(
		    Vec3D(state[0], state[1], state[2]), Vec3D(state[3], state[4], state[5]) * Constants::GAUSS_K, epoch,
		    {outputEpoch}, ng);
		return make_orbit_from_states_vectors(StateVectors(propagated.at(0).R, propagated.at(0).V / Constants::GAUSS_K),
		                                      outputEpoch);
	};
	const auto display = sample(params);
	ASSERT_TRUE(display);
	const auto reference = calculator.computeCometaryDerivatives(params, *orbit, *display, {});
	ASSERT_TRUE(reference.valid);
	Eigen::VectorXd parameterSigmas(10);
	parameterSigmas << 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 3.66e-7, 1.72e-7, 5.64e-8, 8.93;
	const Eigen::VectorXd expectedSigmas = (reference.cotJ * parameterSigmas.asDiagonal()).rowwise().norm();
	for (int repeat = 0; repeat < 4; ++repeat) {
		auto perturbed = params;
		perturbed[repeat] = std::nextafter(params[repeat], std::numeric_limits<double>::infinity());
		const auto derivatives = calculator.computeCometaryDerivatives(perturbed, *orbit, *display, {});
		ASSERT_TRUE(derivatives.valid);
		const Eigen::VectorXd sigmas = (derivatives.cotJ * parameterSigmas.asDiagonal()).rowwise().norm();
		EXPECT_LT((sigmas.array() / expectedSigmas.array() - 1).abs().maxCoeff(), 1e-6);
		for (int col = 6; col < 10; ++col) {
			SCOPED_TRACE(::testing::Message() << "repeat=" << repeat << " col=" << col);
			EXPECT_LT((derivatives.cotJ.col(col) - reference.cotJ.col(col)).norm(),
			          std::max(reference.cotJ.col(col).norm(), col == 9 ? 1e-7 : 1.0) * 2e-5);
		}
	}
	const auto cot = [](const Orbit& value) {
		cometary::Vector6 out;
		out << value.getPerigee(), value.getExcentricity(), value.getInclination(), value.getLongitudeOfAscendingNode(),
		    value.getArgumentOfPericentre(), value.getTrueAnomaly();
		return out;
	};
	for (double step : {0.01, 0.02}) {
		auto plus = params, minus = params;
		plus[9] += step;
		minus[9] -= step;
		cometary::Vector6 difference = cot(*sample(plus)) - cot(*sample(minus));
		for (int row = 3; row < 6; ++row)
			difference[row] = Helpers::wrapPi(difference[row]);
		const cometary::Vector6 derivative = difference / (2 * step);
		EXPECT_LT((reference.cotJ.col(9) - derivative).norm(), std::max(derivative.norm(), 1e-7) * 2e-5);
	}
}

TEST(OrbitFitCometaryUncertaintyTest, ReducedNonGravOutputColumnsMatchFullMaskInBothDirections)
{
	for (const double eccentricity : {0.7, 1.2}) {
		for (const double offset : {-100.0, 100.0}) {
			SCOPED_TRACE(::testing::Message() << "e=" << eccentricity << " offset=" << offset);
			const auto orbit = cometOrbit(eccentricity);
			OrbitOptimizerConfig config;
			config.nonGrav.enabled = true;
			config.nonGrav.model = orbital_core::NonGravModelKind::MarsdenH2O;
			OrbitFitModel model(orbit, {}, config);
			auto& ng = model.propagator().nonGrav();
			ng.setParams(7.5e-8, 0.0, -4.4e-9);
			ng.setDT(-92.85);
			auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
			model.nonGravProfile().extendParams(params, model.differenceSteps());
			OrbitFitProblem problem(model);
			OrbitFitUncertainty calculator(model, problem, fitEpoch + offset);
			const auto states = model.propagator().propagateAllWithVelocity(
			    orbit->getR(), orbit->getV() * Constants::GAUSS_K, fitEpoch, {fitEpoch + offset});
			ASSERT_EQ(states.size(), 1u);
			const auto display = make_orbit_from_states_vectors(
			    StateVectors(states[0].R, states[0].V / Constants::GAUSS_K), fitEpoch + offset);
			ASSERT_TRUE(display);
			const auto cotFull = calculator.computeCometaryDerivatives(params, *orbit, *display, {});
			const auto elemFull = calculator.computeElementDerivatives(params, *orbit, {});
			for (int column : {7, 9}) {
				const double h = column == 9 ? 0.01 : 1e-9;
				const auto sample = [&](double delta) {
					auto state = params;
					state[column] += delta;
					auto force = ng;
					NonGravProfile::updateFromParams(force, state, 6);
					const auto propagated = model.propagator().propagateAllWithVelocity(
					    orbit->getR(), orbit->getV() * Constants::GAUSS_K, fitEpoch, {fitEpoch + offset}, force);
					cometary::Vector6 car;
					const auto& value = propagated.at(0);
					car << value.R.x, value.R.y, value.R.z, value.V.x, value.V.y, value.V.z;
					return cometary::fromCartesian(car).value();
				};
				cometary::Vector6 difference = sample(h) - sample(-h);
				for (int row = 3; row < 6; ++row)
					difference[row] = Helpers::wrapPi(difference[row]);
				const cometary::Vector6 derivative = difference / (2 * h);
				EXPECT_LT((cotFull.cotJ.col(column) - derivative).norm(),
				          std::max(derivative.norm(), column == 9 ? 1e-7 : 1.0) * 2e-5);
			}
			ng.setFitParameterMask(NgParameterMask::A1 | NgParameterMask::A2);
			params.resize(6);
			model.differenceSteps().resize(6);
			model.nonGravProfile().extendParams(params, model.differenceSteps());
			ASSERT_EQ(params.size(), 8u);
			const auto cotSparse = calculator.computeCometaryDerivatives(params, *orbit, *display, {});
			const auto elemSparse = calculator.computeElementDerivatives(params, *orbit, {});
			ASSERT_TRUE(cotFull.valid && cotSparse.valid);
			EXPECT_TRUE(elemSparse.elemJ.allFinite());
			if (eccentricity < 1.0)
				EXPECT_TRUE(elemSparse.qualityElemJ.allFinite());
			for (int k = 0; k < 2; ++k) {
				const int fullColumn = 6 + k;
				const double floor = 1.0;
				EXPECT_LT((cotSparse.cotJ.col(6 + k) - cotFull.cotJ.col(fullColumn)).norm(),
				          std::max(cotFull.cotJ.col(fullColumn).norm(), floor) * 1e-7);
				EXPECT_LT((elemSparse.elemJ.col(6 + k) - elemFull.elemJ.col(fullColumn)).norm(),
				          std::max(elemFull.elemJ.col(fullColumn).norm(), floor) * 1e-7);
				if (eccentricity < 1.0)
					EXPECT_LT((elemSparse.qualityElemJ.col(6 + k) - elemFull.qualityElemJ.col(fullColumn)).norm(),
					          std::max(elemFull.qualityElemJ.col(fullColumn).norm(), floor) * 1e-7);
			}
		}
	}
}

TEST(OrbitFitCometaryUncertaintyTest, LongArcYarkovskyTransportMatchesSmallPerturbations)
{
	const auto orbit = cometOrbit(0.2);
	OrbitOptimizerConfig config;
	config.nonGrav.enabled = true;
	config.nonGrav.model = orbital_core::NonGravModelKind::Yarkovsky;
	OrbitFitModel model(orbit, {}, config);
	auto& ng = model.propagator().nonGrav();
	ng.setFitParameterMask(NgParameterMask::A2);
	ng.setParams(0.0, -1.5e-14, 0.0);
	auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
	model.nonGravProfile().extendParams(params, model.differenceSteps());
	ASSERT_EQ(params.size(), 7u);
	OrbitFitProblem problem(model);
	for (double duration : {-3652.5, 3652.5}) {
		SCOPED_TRACE(duration);
		const double epoch = fitEpoch + duration;
		OrbitFitUncertainty calculator(model, problem, epoch);
		const auto sample = [&](double delta) {
			auto perturbed = params;
			perturbed[6] += delta;
			auto force = ng;
			NonGravProfile::updateFromParams(force, perturbed, 6);
			const auto states = model.propagator().propagateAllWithVelocity(
			    orbit->getR(), orbit->getV() * Constants::GAUSS_K, fitEpoch, {epoch}, force);
			return make_orbit_from_states_vectors(StateVectors(states.at(0).R, states.at(0).V / Constants::GAUSS_K),
			                                      epoch);
		};
		const auto display = sample(0.0);
		ASSERT_TRUE(display);
		const auto cot = calculator.computeCometaryDerivatives(params, *orbit, *display, {});
		const auto elem = calculator.computeElementDerivatives(params, *orbit, {});
		ASSERT_TRUE(cot.valid);
		ASSERT_TRUE(elem.elemJ.allFinite());
		for (double h : {1e-11, 3e-12}) {
			const auto plus = sample(h), minus = sample(-h);
			cometary::Vector6 difference;
			difference << plus->getPerigee() - minus->getPerigee(), plus->getExcentricity() - minus->getExcentricity(),
			    plus->getInclination() - minus->getInclination(),
			    Helpers::wrapPi(plus->getLongitudeOfAscendingNode() - minus->getLongitudeOfAscendingNode()),
			    Helpers::wrapPi(plus->getArgumentOfPericentre() - minus->getArgumentOfPericentre()),
			    Helpers::wrapPi(plus->getTrueAnomaly() - minus->getTrueAnomaly());
			const cometary::Vector6 derivative = difference / (2 * h);
			EXPECT_LT((cot.cotJ.col(6) - derivative).norm(), derivative.norm() * 2e-5);
			const double da = (plus->getSemiMajorAxis() - minus->getSemiMajorAxis()) / (2 * h);
			EXPECT_NEAR(elem.elemJ(0, 6), da, std::abs(da) * 2e-5);
		}
		auto covariance = isotropicCovariance(7);
		covariance.fullCovariance.setZero();
		covariance.formalFullCovariance.setZero();
		covariance.fullCovariance(6, 6) = covariance.formalFullCovariance(6, 6) = 1e-28;
		orbital_core::Uncertainty uncertainty;
		calculator.populateCometaryUncertainty(uncertainty, cot, covariance, *display);
		ASSERT_TRUE(uncertainty.cometaryElementCovariance);
		EXPECT_NEAR(uncertainty.derivedSigmas.sigma_q, std::abs(cot.cotJ(0, 6)) * 1e-14,
		            std::abs(cot.cotJ(0, 6)) * 1e-25);
	}
}

TEST(OrbitFitCometaryUncertaintyTest, YarkovskyTransportSelectsA2SensitivityColumn)
{
	const auto orbit = cometOrbit(0.7);
	OrbitOptimizerConfig config;
	config.nonGrav.enabled = true;
	config.nonGrav.model = orbital_core::NonGravModelKind::InverseSquare;
	OrbitFitModel model(orbit, {}, config);
	auto& ng = model.propagator().nonGrav();
	ng.setParams(0.0, -1.5e-14, 0.0);
	auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
	model.nonGravProfile().extendParams(params, model.differenceSteps());
	OrbitFitProblem problem(model);
	const double epoch = fitEpoch + 200.0;
	const auto states = model.propagator().propagateAllWithVelocity(orbit->getR(), orbit->getV() * Constants::GAUSS_K,
	                                                                fitEpoch, {epoch});
	const auto display =
	    make_orbit_from_states_vectors(StateVectors(states.at(0).R, states.at(0).V / Constants::GAUSS_K), epoch);
	OrbitFitUncertainty calculator(model, problem, epoch);
	const auto full = calculator.computeCometaryDerivatives(params, *orbit, *display, {});
	ASSERT_TRUE(full.valid);
	ng.setModel(NgModel::Yarkovsky);
	params.resize(6);
	model.differenceSteps().resize(6);
	model.nonGravProfile().extendParams(params, model.differenceSteps());
	ASSERT_EQ(params.size(), 7u);
	const auto sparse = calculator.computeCometaryDerivatives(params, *orbit, *display, {});
	ASSERT_TRUE(sparse.valid);
	EXPECT_TRUE(sparse.cotJ.leftCols<6>().isApprox(full.cotJ.leftCols<6>(), 1e-12));
	EXPECT_TRUE(sparse.cotJ.col(6).isApprox(full.cotJ.col(7), 1e-12));
}

TEST(OrbitFitCometaryUncertaintyTest, GravityTransportMatchesIndependentPropagation)
{
	const auto orbit = cometOrbit(0.2);
	OrbitFitModel model(orbit, {}, {});
	OrbitFitProblem problem(model);
	const auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
	for (double duration : {-3652.5, 3652.5}) {
		const double epoch = fitEpoch + duration;
		const auto transported =
		    OrbitFitUncertainty(model, problem, epoch).computeElementDerivatives(params, *orbit, {});
		Eigen::Matrix<double, 6, 6> numerical;
		for (int col = 0; col < 6; ++col) {
			const double h = 1e-4 * std::max(1.0, std::abs(params[col]));
			const auto sample = [&](double delta) {
				auto state = params;
				state[col] += delta;
				const auto states = model.propagator().propagateAllWithVelocity(
				    Vec3D(state[0], state[1], state[2]), Vec3D(state[3], state[4], state[5]) * Constants::GAUSS_K,
				    fitEpoch, {epoch});
				const auto propagated = make_orbit_from_states_vectors(
				    StateVectors(states.at(0).R, states.at(0).V / Constants::GAUSS_K), epoch);
				Eigen::Matrix<double, 6, 1> elements;
				elements << propagated->getSemiMajorAxis(), propagated->getExcentricity(), propagated->getInclination(),
				    propagated->getLongitudeOfAscendingNode(), propagated->getArgumentOfPericentre(),
				    propagated->getTrueAnomaly();
				return elements;
			};
			auto near = (sample(h) - sample(-h)).eval();
			auto far = (sample(2 * h) - sample(-2 * h)).eval();
			for (int row = 3; row < 6; ++row) {
				near[row] = Helpers::wrapPi(near[row]);
				far[row] = Helpers::wrapPi(far[row]);
			}
			numerical.col(col) = (8 * near - far) / (12 * h);
		}
		for (int row = 0; row < 6; ++row) {
			SCOPED_TRACE(::testing::Message() << duration << " row=" << row);
			EXPECT_LT((transported.elemJ.row(row) - numerical.row(row)).norm(), numerical.row(row).norm() * 2e-7);
		}
	}
}

TEST(OrbitFitCometaryUncertaintyTest, InvalidOrMismatchedEpochDoesNotPublishCovariance)
{
	const auto orbit = cometOrbit(1.0);
	ASSERT_TRUE(orbit);
	OrbitFitUncertainty::CometaryDerivatives derivatives;
	derivatives.cotJ.setIdentity(6, 6);
	derivatives.regularQualityJ.setZero(3, 6);
	derivatives.epochJdTdb = fitEpoch + 1.0;
	derivatives.valid = true;
	orbital_core::Uncertainty uncertainty;
	OrbitFitUncertainty::populateCometaryUncertainty(uncertainty, derivatives, isotropicCovariance(6), *orbit);
	EXPECT_FALSE(uncertainty.cometaryElementCovariance);
	EXPECT_FALSE(uncertainty.hasElementCovariance);
	EXPECT_TRUE(uncertainty.empty());
	EXPECT_TRUE(std::isnan(uncertainty.derivedSigmas.sigma_q));
	derivatives.epochJdTdb = fitEpoch;
	derivatives.valid = false;
	OrbitFitUncertainty::populateCometaryUncertainty(uncertainty, derivatives, isotropicCovariance(6), *orbit);
	EXPECT_FALSE(uncertainty.cometaryElementCovariance);
	EXPECT_TRUE(uncertainty.elementCovariance.array().isNaN().all());
}

TEST(OrbitFitCometaryUncertaintyTest, SingularAngularChartIsUnavailable)
{
	const auto orbit = cometOrbit(1.0, -0.7, fitEpoch, 0.0);
	ASSERT_TRUE(orbit);
	OrbitFitModel model(orbit, {}, {});
	OrbitFitProblem problem(model);
	OrbitFitUncertainty calculator(model, problem, fitEpoch);
	const auto derivatives =
	    calculator.computeCometaryDerivatives(OrbitFitModel::stateVectorFromOrbit(*orbit), *orbit, *orbit, {});
	EXPECT_FALSE(derivatives.valid);
}

TEST(OrbitFitCometaryUncertaintyTest, FitEpochAphelionJacobianHasNoDirectNonGravColumns)
{
	const auto orbit = cometOrbit(1.0 - 1e-7, Constants::PI);
	ASSERT_TRUE(orbit);
	EXPECT_LT(orbit->getV().magnitude(), 1e-7);
	OrbitOptimizerConfig config;
	config.nonGrav.enabled = true;
	config.nonGrav.model = orbital_core::NonGravModelKind::InverseSquare;
	OrbitFitModel model(orbit, {}, config);
	OrbitFitProblem problem(model);
	OrbitFitUncertainty calculator(model, problem, fitEpoch);
	auto params = OrbitFitModel::stateVectorFromOrbit(*orbit);
	params.resize(6 + model.propagator().nonGrav().numFitParams(), 0.0);
	ASSERT_GT(params.size(), 6u);
	const auto derivatives = calculator.computeCometaryDerivatives(params, *orbit, *orbit, {});
	ASSERT_TRUE(derivatives.valid);
	cometary::Vector6 cot;
	cot << orbit->getPerigee(), orbit->getExcentricity(), orbit->getInclination(), orbit->getLongitudeOfAscendingNode(),
	    orbit->getArgumentOfPericentre(), orbit->getTrueAnomaly();
	const auto local = cometary::cartesianToElementsJacobian(cot);
	ASSERT_TRUE(local);
	cometary::Matrix6 expected = *local;
	expected.rightCols<3>() *= Constants::GAUSS_K;
	EXPECT_TRUE(derivatives.cotJ.leftCols<6>().isApprox(expected, 1e-13));
	EXPECT_TRUE(derivatives.cotJ.rightCols(params.size() - 6).isZero(0.0));
	EXPECT_TRUE(derivatives.regularQualityJ.rightCols(params.size() - 6).isZero(0.0));
	const auto quality = OrbitFitUncertainty::regularQualityJacobianCot(cot);
	ASSERT_TRUE(quality);
	const Eigen::Matrix<double, 3, Eigen::Dynamic> expectedQuality = *quality * derivatives.cotJ;
	EXPECT_TRUE(derivatives.regularQualityJ.isApprox(expectedQuality, 1e-13))
	    << "actual:\n"
	    << derivatives.regularQualityJ << "\nexpected:\n"
	    << expectedQuality;
}

TEST(OrbitFitCometaryUncertaintyTest, PerihelionTimeJacobianPinsBothAphelionBranches)
{
	for (const double e : {0.999, 1.0 - 1e-7})
		for (const double nu : {Constants::PI - 1e-10, Constants::PI, Constants::PI + 1e-10, -Constants::PI}) {
			SCOPED_TRACE(e);
			SCOPED_TRACE(nu);
			constexpr double q = 2.3;
			cometary::Vector6 cot;
			cot << q, e, 0.41, 0.82, 0.37, nu;
			const auto quality = OrbitFitUncertainty::regularQualityJacobianCot(cot);
			ASSERT_TRUE(quality);
			const double time = cometary::timeSincePerihelion(q, e, nu);
			const double denominator = 1.0 + e * std::cos(nu);
			const double p = q * (1.0 + e);
			const double expectedNu = -p * std::sqrt(p) / (Constants::GAUSS_K * denominator * denominator);
			const double expectedQ = -1.5 * time / q;
			const double angularTerm = q * std::sqrt(q) * std::sqrt(1.0 + e) * std::sin(nu) * (2.0 + e * std::cos(nu)) /
			                           (Constants::GAUSS_K * denominator * denominator);
			const double expectedE = -(1.5 * time - angularTerm) / (1.0 - e);
			EXPECT_NEAR((*quality)(2, 0), expectedQ, 2e-7 * std::abs(expectedQ));
			EXPECT_NEAR((*quality)(2, 1), expectedE, 2e-7 * std::abs(expectedE));
			EXPECT_NEAR((*quality)(2, 5), expectedNu, 2e-7 * std::abs(expectedNu));
			EXPECT_LT((*quality)(2, 5), 0.0);
			EXPECT_DOUBLE_EQ((*quality)(0, 0), 1.0);
			EXPECT_DOUBLE_EQ((*quality)(1, 0), -(1.0 - e) / (q * q));
			EXPECT_DOUBLE_EQ((*quality)(1, 1), -1.0 / q);
		}
}

TEST(OrbitFitCometaryUncertaintyTest, AphelionSigmasDoNotContainAPeriodJump)
{
	for (const double e : {0.999, 1.0 - 1e-7})
		for (const double nu : {Constants::PI, -Constants::PI}) {
			SCOPED_TRACE(e);
			SCOPED_TRACE(nu);
			const auto orbit = cometOrbit(e, nu);
			ASSERT_TRUE(orbit);
			cometary::Vector6 cot;
			cot << orbit->getPerigee(), e, 0.41, 0.82, 0.37, nu;
			const auto quality = OrbitFitUncertainty::regularQualityJacobianCot(cot);
			ASSERT_TRUE(quality);
			OrbitFitUncertainty::CometaryDerivatives derivatives;
			derivatives.cotJ.setIdentity(6, 6);
			derivatives.regularQualityJ = *quality;
			derivatives.epochJdTdb = fitEpoch;
			derivatives.valid = true;
			FitCovariance covariance;
			covariance.formalFullCovariance = Eigen::MatrixXd::Zero(6, 6);
			covariance.formalFullCovariance.diagonal() << 1e-18, 1e-24, 1e-18, 1e-18, 1e-18, 1e-18;
			covariance.fullCovariance = covariance.formalFullCovariance;
			orbital_core::Uncertainty uncertainty;
			OrbitFitUncertainty::populateCometaryUncertainty(uncertainty, derivatives, covariance, *orbit);
			ASSERT_TRUE(uncertainty.cometaryElementCovariance);
			const auto& sig = uncertainty.derivedSigmas;
			const double expectedM = (1.0 + e) * std::sqrt(1.0 + e) / std::sqrt(1.0 - e) * 1e-9;
			EXPECT_NEAR(sig.sigma_M, expectedM, 1e-12 * expectedM);
			EXPECT_LT(sig.sigma_M, 1e-3);
			const double invA = (1.0 - e) / orbit->getPerigee();
			const double n = Constants::GAUSS_K * invA * std::sqrt(invA);
			const double period = Constants::TWOPI / n;
			EXPECT_TRUE(std::isfinite(sig.sigma_Tp));
			EXPECT_GT(sig.sigma_Tp, 0.0);
			EXPECT_LT(sig.sigma_Tp, 1e-3 * period);
			const double expectedN = 1.5 * Constants::GAUSS_K * std::sqrt(invA) * sig.sigma_inv_a;
			EXPECT_NEAR(sig.sigma_n, expectedN, 1e-14 * expectedN);
			const double expectedPeriod = Constants::TWOPI * expectedN / (n * n);
			EXPECT_NEAR(sig.sigma_P, expectedPeriod, 1e-14 * expectedPeriod);
		}
}

TEST(OrbitFitCometaryUncertaintyTest, UnboundAphelionGeometryIsUnavailable)
{
	for (const double e : {1.0, 1.0001}) {
		cometary::Vector6 cot;
		cot << 2.3, e, 0.41, 0.82, 0.37, Constants::PI;
		EXPECT_FALSE(OrbitFitUncertainty::regularQualityJacobianCot(cot));
	}
}
