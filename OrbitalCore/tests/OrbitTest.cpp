#include <gtest/gtest.h>
#include "orbital_core/Orbit.h"
#include "orbital_core/FitPipeline.h"
#include "orbital_core/Types.h"
#include "orbital_core/CometaryCoordinates.h"
#include "TwoBodyPropagation.h"
#include <cmath>

TEST(Orbit, OrbitalGroupsAndObjectDesignation)
{
	using Class = orbital_core::OrbitalClass;
	Orbit orbit;
	orbit.setType(OrbitType::Elliptical);
	const auto classify = [&](double a, double e, std::string_view name = {}) {
		orbit.setSemiMajorAxis(a);
		orbit.setExcentricity(e);
		return orbital_core::orbitalClassification(orbit, name);
	};
	EXPECT_EQ(classify(0.7, 0.1), Class::Atira);
	EXPECT_EQ(classify(0.9, 0.2), Class::Aten);
	EXPECT_EQ(classify(1.2, 0.3, "K20V00V"), Class::Apollo);
	EXPECT_EQ(classify(1.5, 0.25), Class::Amor);
	EXPECT_EQ(classify(2.0, 0.3), Class::MarsCrossing);
	EXPECT_EQ(classify(1.9, 0.05), Class::InnerMainBelt);
	EXPECT_EQ(classify(2.5, 0.1), Class::MainBelt);
	EXPECT_EQ(classify(4.0, 0.1), Class::OuterMainBelt);
	EXPECT_EQ(classify(5.28, 0.0243), Class::JupiterTrojanCandidate);
	EXPECT_EQ(classify(5.21, 0.139), Class::JupiterTrojanCandidate);
	EXPECT_EQ(classify(6.0, 0.1), Class::Centaur);
	EXPECT_EQ(classify(39.6, 0.252), Class::Tno);
	EXPECT_EQ(classify(544.0, 0.86), Class::Tno);
	EXPECT_EQ(classify(1.0, 0.3), Class::Apollo);
	EXPECT_EQ(classify(0.983, 0.0), Class::Aten);
	EXPECT_EQ(classify(4.6, 0.0), Class::JupiterTrojanCandidate);
	EXPECT_EQ(classify(4.6, 0.3), Class::Unknown);
	EXPECT_EQ(classify(5.5, 0.0), Class::Centaur);
	EXPECT_EQ(classify(30.1, 0.0), Class::Tno);
	EXPECT_EQ(classify(8.0, 0.9), Class::Apollo);
	EXPECT_EQ(classify(1.2, 0.3, "0001PI35P010"), Class::Comet);
	EXPECT_EQ(classify(1.2, 0.3, "    CK20F030"), Class::Comet);
	EXPECT_EQ(classify(1.2, 0.3, "C/2020 F3"), Class::Comet);
	EXPECT_EQ(classify(1.2, 0.3, "2I"), Class::Interstellar);
	EXPECT_EQ(classify(1.2, 0.3, "A/2017 U1"), Class::Apollo);
	orbit.setType(OrbitType::Hyperbolic);
	EXPECT_EQ(orbital_core::orbitalClassification(orbit, "tracklet"), Class::Unknown);
	orbit.setExcentricity(1.2);
	EXPECT_EQ(orbital_core::orbitalClassification(orbit, "tracklet"), Class::Hyperbolic);
	EXPECT_EQ(orbital_core::orbitalClassification(orbit, "C/2019 Q4"), Class::HyperbolicComet);
	EXPECT_EQ(orbital_core::orbitalClassification(orbit, "2I"), Class::Interstellar);
	orbit.setType(OrbitType::Parabolic);
	orbit.setExcentricity(1.0);
	EXPECT_EQ(orbital_core::orbitalClassification(orbit, "tracklet"), Class::Parabolic);
	EXPECT_EQ(orbital_core::orbitalClassification(orbit, "C/2020 F3"), Class::ParabolicComet);
}

TEST(Orbit, CometPeriodAndTisserandClassificationsAreIndependent)
{
	using Period = orbital_core::CometPeriodClass;
	using Tisserand = orbital_core::CometTisserandClass;
	Orbit orbit;
	orbit.setType(OrbitType::Elliptical);
	const auto classify = [&](double a, double e, double iDegrees, double periodDays) {
		orbit.setSemiMajorAxis(a);
		orbit.setExcentricity(e);
		orbit.setInclination(Convertors::degreesToRadians(iDegrees));
		orbit.setOrbitalPeriod(periodDays);
		return orbital_core::cometOrbitClasses(orbit, "1P");
	};
	const auto halley = classify(17.9, 0.968, 162.0, 27700.0);
	EXPECT_EQ(halley.period, Period::HalleyType);
	EXPECT_EQ(halley.tisserand, Tisserand::Unknown);
	const auto encke = classify(2.22, 0.847, 11.3, 1210.0);
	EXPECT_EQ(encke.period, Period::JupiterFamily);
	EXPECT_EQ(encke.tisserand, Tisserand::EnckeType);
	const auto churyumov = classify(3.46, 0.641, 7.04, 2350.0);
	EXPECT_EQ(churyumov.period, Period::JupiterFamily);
	EXPECT_EQ(churyumov.tisserand, Tisserand::JupiterFamily);
	const auto chironType = classify(10.0, 0.0, 0.0, 100.0 * 365.25);
	EXPECT_EQ(chironType.period, Period::HalleyType);
	EXPECT_EQ(chironType.tisserand, Tisserand::ChironType);
	EXPECT_EQ(classify(10.0, 0.0, 0.0, 20.0 * 365.25).period, Period::HalleyType);
	EXPECT_EQ(classify(10.0, 0.0, 0.0, 200.0 * 365.25).period, Period::HalleyType);
	EXPECT_EQ(classify(10.0, 0.0, 0.0, 201.0 * 365.25).period, Period::Unknown);
	const auto notComet = orbital_core::cometOrbitClasses(orbit, "K20V00V");
	EXPECT_EQ(notComet.period, Period::Unknown);
	EXPECT_EQ(notComet.tisserand, Tisserand::Unknown);
}

TEST(Orbit, FitResultClassRequiresConvergenceForOrbitalGroups)
{
	using Class = orbital_core::OrbitalClass;
	orbital_core::FitPipelineResult result;
	result.fit.orbit = std::make_shared<Orbit>();
	result.fit.orbit->setType(OrbitType::Elliptical);
	result.fit.orbit->setSemiMajorAxis(1.2);
	result.fit.orbit->setExcentricity(0.3);
	EXPECT_EQ(result.orbitalClass(), Class::Unknown);
	result.fit.converged = true;
	EXPECT_EQ(result.orbitalClass(), Class::Apollo);
	result.observations.classificationDesignation = "1P";
	result.fit.converged = false;
	EXPECT_EQ(result.orbitalClass(), Class::Comet);
	result.observations.classificationDesignation.clear();
	result.observations.format = "ADES";
	result.observations.objectName = "2I";
	EXPECT_EQ(result.orbitalClass(), Class::Unknown);
}

TEST(Orbit, DefaultConstructor)
{
	Orbit orb;
	EXPECT_EQ(orb.getType(), OrbitType::Undef);
	EXPECT_EQ(orb.getEpoch(), 0.0);
	EXPECT_EQ(orb.getExcentricity(), 0.0);
	EXPECT_EQ(orb.getSemiMajorAxis(), 0.0);
}

TEST(Orbit, SettersAndGetters)
{
	Orbit orb;
	orb.setType(OrbitType::Elliptical);
	orb.setEpoch(Constants::J2000_JD);
	orb.setExcentricity(0.5);
	orb.setSemiMajorAxis(2.0);
	orb.setInclination(0.3);
	orb.setLongitudeOfAscendingNode(1.2);
	orb.setArgumentOfPericentre(0.8);
	orb.setTrueAnomaly(1.5);
	orb.setMeanAnomaly(1.0);

	EXPECT_EQ(orb.getType(), OrbitType::Elliptical);
	EXPECT_EQ(orb.getEpoch(), Constants::J2000_JD);
	EXPECT_EQ(orb.getExcentricity(), 0.5);
	EXPECT_EQ(orb.getSemiMajorAxis(), 2.0);
	EXPECT_EQ(orb.getInclination(), 0.3);
	EXPECT_EQ(orb.getLongitudeOfAscendingNode(), 1.2);
	EXPECT_EQ(orb.getArgumentOfPericentre(), 0.8);
	EXPECT_EQ(orb.getTrueAnomaly(), 1.5);
	EXPECT_EQ(orb.getMeanAnomaly(), 1.0);
}

TEST(Orbit, StateVectors)
{
	Orbit orb;
	Vec3D R(1.0, 2.0, 3.0);
	Vec3D V(0.1, 0.2, 0.3);
	orb.setR(R);
	orb.setV(V);

	EXPECT_EQ(orb.getR().x, 1.0);
	EXPECT_EQ(orb.getR().y, 2.0);
	EXPECT_EQ(orb.getR().z, 3.0);
	EXPECT_EQ(orb.getV().x, 0.1);
	EXPECT_EQ(orb.getV().y, 0.2);
	EXPECT_EQ(orb.getV().z, 0.3);
}

TEST(Orbit, CopyConstructor)
{
	Orbit orb;
	orb.setType(OrbitType::Elliptical);
	orb.setEpoch(Constants::J2000_JD);
	orb.setExcentricity(0.5);
	orb.setSemiMajorAxis(2.0);
	orb.setInclination(0.3);
	orb.setR(Vec3D(1.0, 2.0, 3.0));
	orb.setV(Vec3D(0.1, 0.2, 0.3));

	Orbit copy(orb);
	EXPECT_EQ(copy.getType(), OrbitType::Elliptical);
	EXPECT_EQ(copy.getEpoch(), Constants::J2000_JD);
	EXPECT_EQ(copy.getExcentricity(), 0.5);
	EXPECT_EQ(copy.getSemiMajorAxis(), 2.0);
	EXPECT_EQ(copy.getR().x, 1.0);
	EXPECT_EQ(copy.getV().z, 0.3);
}

TEST(Orbit, SigmaFields)
{
	Orbit orb;
	orb.setSigmaExcentricity(0.01);
	orb.setSigmaSemiMajorAxis(0.1);
	orb.setSigmaInclination(0.001);
	orb.setSigmaLongitudeOfAscendingNode(0.002);
	orb.setSigmaArgumentOfPericentre(0.003);
	orb.setSigmaTrueAnomaly(0.004);

	EXPECT_EQ(orb.getSigmaExcentricity(), 0.01);
	EXPECT_EQ(orb.getSigmaSemiMajorAxis(), 0.1);
	EXPECT_EQ(orb.getSigmaInclination(), 0.001);
	EXPECT_EQ(orb.getSigmaLongitudeOfAscendingNode(), 0.002);
	EXPECT_EQ(orb.getSigmaArgumentOfPericentre(), 0.003);
	EXPECT_EQ(orb.getSigmaTrueAnomaly(), 0.004);
}

TEST(Orbit, TypeString)
{
	Orbit orb;
	EXPECT_EQ(orb.getTypeString(), "Undefined");
	orb.setType(OrbitType::Elliptical);
	EXPECT_EQ(orb.getTypeString(), "Elliptical");
	orb.setType(OrbitType::Hyperbolic);
	EXPECT_EQ(orb.getTypeString(), "Hyperbolic");
}

TEST(MakeOrbit, PerfectlyCircularReturnsNull)
{
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, 1.0, 0.0);
	auto orb = make_orbit_from_states_vectors(StateVectors(R, V), Constants::J2000_JD);
	EXPECT_EQ(orb, nullptr);
}

TEST(MakeOrbit, NearCircularOrbit)
{
	double v = 1.001;
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, v, 0.0);
	auto orb = make_orbit_from_states_vectors(StateVectors(R, V), Constants::J2000_JD);

	ASSERT_NE(orb, nullptr);
	EXPECT_EQ(orb->getType(), OrbitType::Elliptical);
	EXPECT_NEAR(orb->getSemiMajorAxis(), 1.0, 0.1);
	EXPECT_LT(orb->getExcentricity(), 0.01);
}

TEST(MakeOrbit, EllipticalOrbit)
{
	double a = 2.0, e = 0.5;
	double r = a * (1.0 - e);
	double v = std::sqrt(2.0 / r - 1.0 / a);

	auto orb =
	    make_orbit_from_states_vectors(StateVectors(Vec3D(r, 0.0, 0.0), Vec3D(0.0, v, 0.0)), Constants::J2000_JD);

	ASSERT_NE(orb, nullptr);
	EXPECT_EQ(orb->getType(), OrbitType::Elliptical);
	EXPECT_NEAR(orb->getSemiMajorAxis(), a, 1e-10);
	EXPECT_NEAR(orb->getExcentricity(), e, 1e-10);
	EXPECT_NEAR(orb->getPerigee(), r, 1e-10);
}

TEST(MakeOrbit, HyperbolicOrbit)
{
	double v_hyp = 2.0;
	auto orb =
	    make_orbit_from_states_vectors(StateVectors(Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, v_hyp, 0.0)), Constants::J2000_JD);

	ASSERT_NE(orb, nullptr);
	EXPECT_EQ(orb->getType(), OrbitType::Hyperbolic);
	EXPECT_GT(orb->getExcentricity(), 1.0);
	EXPECT_LT(orb->getSemiMajorAxis(), 0.0);
}

TEST(MakeOrbit, PreservesStateVectors)
{
	Vec3D R(0.5, -1.2, 0.8);
	Vec3D V(-0.3, 0.6, 0.1);
	auto orb = make_orbit_from_states_vectors(StateVectors(R, V), 2458000.0);

	ASSERT_NE(orb, nullptr);
	EXPECT_NEAR(orb->getR().x, R.x, 1e-15);
	EXPECT_NEAR(orb->getR().y, R.y, 1e-15);
	EXPECT_NEAR(orb->getR().z, R.z, 1e-15);
	EXPECT_NEAR(orb->getV().x, V.x, 1e-15);
	EXPECT_NEAR(orb->getV().y, V.y, 1e-15);
	EXPECT_NEAR(orb->getV().z, V.z, 1e-15);
}

TEST(MakeOrbit, InclinedOrbit)
{
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.01, 0.0, 1.0);
	auto orb = make_orbit_from_states_vectors(StateVectors(R, V), Constants::J2000_JD);

	ASSERT_NE(orb, nullptr);
	EXPECT_NEAR(orb->getInclination(), Constants::PI / 2.0, 0.02);
}

TEST(MakeOrbit, NearParabolicElementsRemainRegular)
{
	using namespace orbital_core;
	for (double e : {1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10}) {
		cometary::Vector6 cot;
		cot << 10.9, e, 1.6, 3.2, 5.7, -1.2;
		const auto state = cometary::toCartesian(cot);
		ASSERT_TRUE(state);
		const Vec3D r((*state)(0), (*state)(1), (*state)(2));
		const Vec3D v = Vec3D((*state)(3), (*state)(4), (*state)(5)) / Constants::GAUSS_K;
		const auto orbit = make_orbit_from_states_vectors(StateVectors(r, v), Constants::J2000_JD);
		ASSERT_TRUE(orbit);
		EXPECT_NEAR(orbit->getExcentricity(), e, 2.0e-15);
		EXPECT_NEAR(orbit->getPerigee(), cot(0), 5.0e-14);
		EXPECT_NEAR(orbit->getReciprocalSemiMajorAxis(), (1.0 - e) / cot(0), 2.0e-16);
		EXPECT_TRUE(std::isfinite(orbit->getTimeOfPerihelion()));
		EXPECT_DOUBLE_EQ(orbit->getR().x, r.x);
		EXPECT_DOUBLE_EQ(orbit->getV().z, v.z);
		const Orbit copy(*orbit);
		EXPECT_DOUBLE_EQ(copy.getTimeOfPerihelion(), orbit->getTimeOfPerihelion());
		EXPECT_DOUBLE_EQ(copy.getReciprocalSemiMajorAxis(), orbit->getReciprocalSemiMajorAxis());
		const auto advanced = Calculations::propagateKeplerState(StateVectors(r, v), 100.0);
		ASSERT_TRUE(advanced.isValid());
		const auto recovered = Calculations::propagateKeplerState(advanced, -100.0);
		ASSERT_TRUE(recovered.isValid());
		EXPECT_LT((recovered.first - r).magnitude(), 1.0e-11);
		EXPECT_LT((recovered.second - v).magnitude(), 1.0e-11);
	}
}

TEST(MakeOrbit, ExactParabolaHasNoFiniteSemiMajorAxis)
{
	const auto orbit =
	    make_orbit_from_states_vectors(StateVectors(Vec3D(2.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0)), Constants::J2000_JD);
	ASSERT_TRUE(orbit);
	EXPECT_EQ(orbit->getType(), OrbitType::Parabolic);
	EXPECT_EQ(orbit->getTypeString(), "Parabolic");
	EXPECT_DOUBLE_EQ(orbit->getExcentricity(), 1.0);
	EXPECT_DOUBLE_EQ(orbit->getPerigee(), 2.0);
	EXPECT_DOUBLE_EQ(orbit->getReciprocalSemiMajorAxis(), 0.0);
	EXPECT_DOUBLE_EQ(orbit->getTimeOfPerihelion(), Constants::J2000_JD);
	EXPECT_TRUE(std::isinf(orbit->getSemiMajorAxis()));
	EXPECT_TRUE(std::isnan(orbit->getOrbitalPeriod()));
	EXPECT_TRUE(std::isnan(orbit->getMeanAnomaly()));
}

TEST(StateVectors, DefaultIsInvalid)
{
	StateVectors sv;
	EXPECT_FALSE(sv.isValid());
}

TEST(StateVectors, ValidVectors)
{
	StateVectors sv(Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0));
	EXPECT_TRUE(sv.isValid());
	EXPECT_EQ(sv.first.x, 1.0);
	EXPECT_EQ(sv.second.y, 1.0);
}

TEST(StateVectors, CopyAndMove)
{
	StateVectors orig(Vec3D(1.0, 2.0, 3.0), Vec3D(4.0, 5.0, 6.0));

	StateVectors copy(orig);
	EXPECT_EQ(copy.first.x, 1.0);
	EXPECT_EQ(copy.second.z, 6.0);

	StateVectors moved(std::move(copy));
	EXPECT_EQ(moved.first.y, 2.0);
}
