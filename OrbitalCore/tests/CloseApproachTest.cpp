#include <gtest/gtest.h>

#include "orbital_core/Constants.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/Orbit.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/CloseApproach.h"
#include "orbital_core/hazard/CloseApproachSearch.h"
#include "NBodyPropagator.h"
#include "hazard/CloseApproachTargetPlaneDifferential.h"
#include "hazard/TargetPlaneDistanceBounds.h"
#include "orbital_core/hazard/PlanetData.h"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <stdexcept>

TEST(TargetPlaneDistanceBoundsTest, PointCircleAndFilledInterior)
{
	using orbital_core::hazard::detail::ellipseDistanceBounds;
	const auto point = ellipseDistanceBounds({3, 4}, {0, 0});
	ASSERT_TRUE(point);
	EXPECT_DOUBLE_EQ(point->minimum, 5);
	EXPECT_DOUBLE_EQ(point->maximum, 5);
	const auto circle = ellipseDistanceBounds({3, 4}, {2, 2});
	ASSERT_TRUE(circle);
	EXPECT_DOUBLE_EQ(circle->minimum, 3);
	EXPECT_DOUBLE_EQ(circle->maximum, 7);
	const auto interior = ellipseDistanceBounds({1, 1}, {5, 2});
	ASSERT_TRUE(interior);
	EXPECT_DOUBLE_EQ(interior->minimum, 0);
	const auto centered = ellipseDistanceBounds({0, 0}, {5, 2});
	ASSERT_TRUE(centered);
	EXPECT_DOUBLE_EQ(centered->minimum, 0);
	EXPECT_DOUBLE_EQ(centered->maximum, 5);
}

TEST(TargetPlaneDistanceBoundsTest, DegenerateSegmentDoesNotClipOffsetToZero)
{
	using orbital_core::hazard::detail::ellipseDistanceBounds;
	for (double minor : {0.0, 1e-12}) {
		const auto bounds = ellipseDistanceBounds({5, 2}, {10, minor});
		ASSERT_TRUE(bounds);
		EXPECT_NEAR(bounds->minimum, 2, 1e-11);
		EXPECT_NEAR(bounds->maximum, std::hypot(15.0, 2.0), 1e-11);
	}
	const auto external = ellipseDistanceBounds({12, 2}, {10, 0});
	ASSERT_TRUE(external);
	EXPECT_DOUBLE_EQ(external->minimum, std::hypot(2.0, 2.0));
	EXPECT_DOUBLE_EQ(external->maximum, std::hypot(22.0, 2.0));
}

TEST(TargetPlaneDistanceBoundsTest, MaximumHardCaseOnMinorAxis)
{
	const auto bounds = orbital_core::hazard::detail::ellipseDistanceBounds({0, 1}, {5, 2});
	ASSERT_TRUE(bounds);
	const double u = 2.0 / 21.0;
	EXPECT_NEAR(bounds->maximum, std::hypot(5 * std::sqrt(1 - u * u), 1 + 2 * u), 1e-14);
	const auto outside = orbital_core::hazard::detail::ellipseDistanceBounds({0, 20}, {5, 2});
	ASSERT_TRUE(outside);
	EXPECT_NEAR(outside->minimum, 18, 1e-13);
	EXPECT_NEAR(outside->maximum, 22, 1e-13);
}

TEST(TargetPlaneDistanceBoundsTest, MatchesIndependentDenseBoundarySearch)
{
	using orbital_core::hazard::detail::ellipseDistanceBounds;
	for (const Eigen::Vector2d center : {Eigen::Vector2d(7, 3), Eigen::Vector2d(-2, 8), Eigen::Vector2d(0, .9)}) {
		for (const Eigen::Vector2d axes : {Eigen::Vector2d(5, 2), Eigen::Vector2d(10, .01), Eigen::Vector2d(1, .99)}) {
			const auto bounds = ellipseDistanceBounds(center, axes);
			ASSERT_TRUE(bounds);
			double low = std::numeric_limits<double>::infinity(), high = 0;
			for (int i = 0; i < 200000; ++i) {
				const double theta = 2 * Constants::PI * i / 200000;
				const double distance =
				    (center + Eigen::Vector2d(axes.x() * std::cos(theta), axes.y() * std::sin(theta))).norm();
				low = std::min(low, distance);
				high = std::max(high, distance);
			}
			if ((center.array() / axes.array()).matrix().norm() <= 1)
				low = 0;
			EXPECT_NEAR(bounds->minimum, low, 2e-7);
			EXPECT_NEAR(bounds->maximum, high, 2e-7);
		}
	}
}

TEST(TargetPlaneDistanceBoundsTest, RescalingAndThinEllipse)
{
	using orbital_core::hazard::detail::ellipseDistanceBounds;
	const Eigen::Vector2d center(5, 1e-9), axes(10, 1e-12);
	const auto original = ellipseDistanceBounds(center, axes);
	ASSERT_TRUE(original);
	EXPECT_GT(original->minimum, 9.98e-10);
	EXPECT_LT(original->minimum, 1e-9);
	for (double scale : {1e-100, 1e-6, 1e6, 1e100}) {
		const auto scaled = ellipseDistanceBounds(center * scale, axes * scale);
		ASSERT_TRUE(scaled);
		EXPECT_NEAR(scaled->minimum / scale, original->minimum, 1e-23);
		EXPECT_NEAR(scaled->maximum / scale, original->maximum, 1e-13);
	}
}

TEST(TargetPlaneDistanceBoundsTest, RejectsInvalidInput)
{
	using orbital_core::hazard::detail::ellipseDistanceBounds;
	EXPECT_FALSE(ellipseDistanceBounds({0, 0}, {1, 2}));
	EXPECT_FALSE(ellipseDistanceBounds({0, 0}, {1, -1}));
	EXPECT_FALSE(ellipseDistanceBounds({NAN, 0}, {1, 1}));
	EXPECT_FALSE(ellipseDistanceBounds({0, 0}, {INFINITY, 1}));
}

TEST(TargetPlaneDistanceBoundsTest, InverseFocusingRecoversPericentre)
{
	using orbital_core::hazard::detail::pericentreFromImpactParameter;
	const double mu = orbital_core::hazard::planetMuKm3PerS2(orbital_core::hazard::Planet::Earth);
	for (double speed : {1.0, 5.0, 40.0}) {
		const double alpha = mu / (speed * speed);
		EXPECT_DOUBLE_EQ(pericentreFromImpactParameter(0, alpha), 0);
		for (double q : {1e-6, 1000.0, 6378.137, 1e9}) {
			const double b = std::sqrt(q * (q + 2 * alpha));
			EXPECT_NEAR(pericentreFromImpactParameter(b, alpha), q, q * 1e-14);
		}
	}
}

TEST(CloseApproachPassageTest, ReportsUnfinishedPassageAtRequestedEndpoint)
{
	using namespace orbital_core::hazard;

	class FixedEarth final : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, Vec3D(1, 0, 0), Vec3D()}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "fixed-earth-passage-test"; }
	} earth;

	constexpr double epoch = 2460000.5;
	for (double direction : {-1.0, 1.0}) {
		SCOPED_TRACE(direction);
		NBodyPropagator propagator;
		propagator.setAllPlanetsEnabled(false);
		propagator.setUseSb441(false);
		CloseApproachConfig cfg;
		cfg.horizonYears = direction > 0 ? 2.75 / 365.25 : 1e-9;
		cfg.lookbackYears = direction < 0 ? 2.75 / 365.25 : 0;
		cfg.thresholdAu = cfg.encounterRadiusAu = .02;
		cfg.sampleStepDay = 1;
		std::vector<CloseApproachPassage> passages{{99, 0, 0}};
		const auto events = findCloseApproachesFromState(Vec3D(1, -.03, 0), Vec3D(0, direction * .004, 0), epoch,
		                                                 propagator, propagator.nonGrav(), cfg, &earth, passages);
		EXPECT_TRUE(events.empty());
		ASSERT_EQ(passages.size(), 1U);
		EXPECT_EQ(passages.front().encounterIndex, 0);
		EXPECT_NEAR(passages.front().firstInsideJd, epoch + direction * 2.75, 1e-8);
		EXPECT_NEAR(passages.front().lastInsideJd, epoch + direction * 2.75, 1e-8);

		cfg.horizonYears = 0;
		EXPECT_TRUE(findCloseApproachesFromState(Vec3D(1, -.03, 0), Vec3D(0, direction * .004, 0), epoch, propagator,
		                                         propagator.nonGrav(), cfg, &earth, passages)
		                .empty());
		EXPECT_TRUE(passages.empty());
	}
}

using namespace orbital_core::hazard;

TEST(CloseApproachDifferentialTest, FindsMinimumOfVariationalTrajectory)
{
	class FixedEarth final : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, Vec3D(1, 0, 0), Vec3D()}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "fixed-earth-minimum-test"; }
	} earth;

	NBodyPropagator propagator;
	propagator.setAllPlanetsEnabled(false);
	propagator.setUseSb441(false);
	constexpr double epoch = 2460000.5;
	CloseApproachEvent event;
	event.planet = Planet::Earth;
	event.jd = epoch;
	event.rPlanet = Vec3D(1, 0, 0);
	event.rTarget = Vec3D(1.0001, 0, 0);
	event.vTarget = Vec3D(0.000046, 0.01, 0.001);
	event.distAu = 0.0001;
	CloseApproachTargetPlaneDifferentialInput input;
	input.event = &event;
	input.parameters.resize(6);
	input.parameters << event.rTarget.x, event.rTarget.y, event.rTarget.z, event.vTarget.x, event.vTarget.y,
	    event.vTarget.z;
	input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6) * 1e-16;
	input.parameterStepPerSigma = Eigen::VectorXd::Zero(6);
	input.parameterStepPerSigma(0) = 1e-8;
	input.nominalNonGrav = &propagator.nonGrav();
	input.propagator = &propagator;
	input.eph = &earth;
	input.planeKind = EncounterPlaneKind::OpikTargetPlane;
	input.epochTdb = epoch;
	const auto trace = makeCloseApproachTargetPlaneTracePoint(input, {});
	ASSERT_TRUE(trace);
	const Vec3D r = trace->event.rTarget - trace->event.rPlanet;
	const Vec3D v = trace->event.vTarget - trace->event.vPlanet;
	EXPECT_NEAR(r.dot(v) / r.magnitude(), 0.0, 5e-10);
	EXPECT_LT(trace->event.jd, epoch);
	EXPECT_GT(trace->event.jd, epoch - 0.0001);
	EXPECT_DOUBLE_EQ(trace->event.jd, trace->trace.jd);
	EXPECT_TRUE(trace->dCoordsDParams.array().isFinite().all());
}

TEST(CloseApproachDifferentialTest, ExternalEventDerivativeIncludesMovingMinimum)
{
	class MovingEarth final : public IEphemerisSource {
	   public:
		EphemerisRecord queryJD(EphemerisPlanet, double jd) override
		{
			const double t = jd - 2460000.5;
			const Vec3D v(0.001, 0.002, 0.0003), a(0.0001, -0.0002, 0.00003);
			return {jd, Vec3D(1, 0, 0) + v * t + a * (0.5 * t * t), v + a * t};
		}

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "moving-earth-event-test"; }
	} earth;

	NBodyPropagator propagator(IntegratorKind::RA15);
	propagator.setAllPlanetsEnabled(false);
	propagator.setUseSb441(false);
	constexpr double epoch = 2460000.5;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	CloseApproachEvent event;
	event.planet = Planet::Earth;
	event.jd = epoch;
	const auto planet = earth.queryJD(EphemerisPlanet::Earth, epoch);
	event.rPlanet = planet.R;
	event.vPlanet = planet.V;
	event.rTarget = planet.R + Vec3D(0.0001, 0, 0);
	event.vTarget = planet.V + Vec3D(0, 0.01, 0.001);
	event.distAu = 0.0001;
	CloseApproachTargetPlaneDifferentialInput input;
	input.event = &event;
	input.parameters.resize(6);
	input.parameters << event.rTarget.x, event.rTarget.y, event.rTarget.z, event.vTarget.x, event.vTarget.y,
	    event.vTarget.z;
	input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6) * 1e-16;
	input.parameterStepPerSigma = Eigen::VectorXd::Zero(6);
	input.parameterStepPerSigma(0) = 1e-8;
	input.nominalNonGrav = &propagator.nonGrav();
	input.propagator = &propagator;
	input.eph = &earth;
	input.planeKind = EncounterPlaneKind::OpikTargetPlane;
	input.epochTdb = epoch;
	const auto center = makeCloseApproachTargetPlaneTracePoint(input, {});
	ASSERT_TRUE(center);
	EXPECT_GT(center->event.distAu * auKm, planetRadiusKm(Planet::Earth));
	for (int column = 0; column < 6; ++column) {
		SCOPED_TRACE(column);
		const double step = column < 3 ? 10.0 / auKm : 0.01 * 86400.0 / auKm;
		auto plus = input, minus = input;
		plus.parameters(column) += step;
		minus.parameters(column) -= step;
		const auto p = makeCloseApproachTargetPlaneTracePoint(plus, {});
		const auto m = makeCloseApproachTargetPlaneTracePoint(minus, {});
		ASSERT_TRUE(p);
		ASSERT_TRUE(m);
		const Eigen::Vector2d finiteDifference = (p->trace.coordsKm - m->trace.coordsKm) / (2.0 * step);
		const Eigen::Vector2d analytic = center->dCoordsDParams.col(column);
		EXPECT_LT((analytic - finiteDifference).norm(), 5e-4 * std::max(1.0, finiteDifference.norm()));
	}
}

TEST(CloseApproachConfigTest, EncounterRadiiAreIndependent)
{
	CloseApproachConfig cfg;
	EXPECT_DOUBLE_EQ(cfg.encounterRadiusAu, Constants::EARTH_TARGET_PLANE_AU);
	EXPECT_DOUBLE_EQ(cfg.thresholdAu, Constants::EARTH_CLOSE_APPROACH_REPORT_AU);
	EXPECT_DOUBLE_EQ(Constants::EARTH_LOCAL_DYNAMICS_AU, 0.1);

	cfg.thresholdAu = 0.05;
	EXPECT_DOUBLE_EQ(cfg.encounterRadiusAu, 0.1);
}

TEST(CloseApproachConfigTest, NominalCatalogueUsesPerBodyReportingLimits)
{
	EXPECT_DOUBLE_EQ(defaultCloseApproachLimitAu(Planet::Mercury), 0.10);
	EXPECT_DOUBLE_EQ(defaultCloseApproachLimitAu(Planet::Venus), 0.15);
	EXPECT_DOUBLE_EQ(defaultCloseApproachLimitAu(Planet::Earth), 0.50);
	EXPECT_DOUBLE_EQ(defaultCloseApproachLimitAu(Planet::Mars), 0.10);
	EXPECT_DOUBLE_EQ(defaultCloseApproachLimitAu(Planet::Moon), 0.01);
	for (Planet body : {Planet::Jupiter, Planet::Saturn, Planet::Uranus, Planet::Neptune})
		EXPECT_DOUBLE_EQ(defaultCloseApproachLimitAu(body), 2.00);
}

class CloseApproachTest : public ::testing::Test {
   protected:
	void SetUp() override
	{
		if (getExecutableDir().empty()) {
			const auto srcDir = std::filesystem::path(__FILE__).parent_path().parent_path();
			const std::filesystem::path candidates[] = {
			    srcDir.parent_path() / "data",
			    srcDir / "data",
			    srcDir / "release",
			    srcDir / "debug",
			};
			for (const auto& d : candidates) {
				if (std::filesystem::exists(d / "linux_p1550p2650.440")) {
					getExecutableDir() = d;
					break;
				}
			}
		}

		m_de440 = std::make_unique<DE440Source>();
		if (m_de440->isAvailable())
			setDefaultEphemerisSource(m_de440.get());

		hasDE440 = m_de440->isAvailable();
	}

	void TearDown() override { setDefaultEphemerisSource(nullptr); }

	Orbit_p makeApophis() const
	{
		const Vec3D R(-9.817715004807913e-01, +4.911298271092577e-01, -4.945859325816258e-02);
		const Vec3D V(-3.762046598819735e-01, -7.703087434762487e-01, +3.211201769300336e-02);
		return make_orbit_from_states_vectors({R, V}, 2459266.348);
	}

	Orbit_p makeCeres() const
	{
		const Vec3D R(-2.609107, -0.5926373, +0.4751728);
		const Vec3D V(+1.660797e-3, -1.071961e-2, +5.094957e-5);
		return make_orbit_from_states_vectors({R, V * (1.0 / Constants::GAUSS_K)}, 2459266.348);
	}

	std::unique_ptr<DE440Source> m_de440;
	NBodyPropagator m_propagator;
	bool hasDE440{false};
};

TEST_F(CloseApproachTest, Apophis2029_FlybyDetected)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	auto apophis = makeApophis();
	ASSERT_TRUE(apophis);

	CloseApproachConfig cfg;
	cfg.horizonYears = 20.0;
	cfg.thresholdAu = 0.01;
	cfg.sampleStepDay = 1.0;
	cfg.planet = Planet::Earth;

	const auto events = findCloseApproaches(*apophis, m_propagator, cfg);
	ASSERT_FALSE(events.empty()) << "no Earth close approaches found in 20-yr window";

	auto closest =
	    std::min_element(events.begin(), events.end(),
	                     [](const CloseApproachEvent& a, const CloseApproachEvent& b) { return a.distAu < b.distAu; });

	EXPECT_NEAR(closest->jd, 2462240.4, 30.0) << "TCA off by more than 30 days from 2029-04-13";

	constexpr double kKmToAu = 1.0 / 149597870.7;
	EXPECT_GT(closest->distAu, 10000.0 * kKmToAu);
	EXPECT_LT(closest->distAu, 100000.0 * kKmToAu);

	EXPECT_GT(closest->vInfKmS, 4.5);
	EXPECT_LT(closest->vInfKmS, 7.5);

	EXPECT_GT(closest->rTarget.magnitude(), 0.5);
	EXPECT_GT(closest->vTarget.magnitude(), 0.0);
	EXPECT_GT(closest->rPlanet.magnitude(), 0.5);
	EXPECT_GT(closest->vPlanet.magnitude(), 0.0);
	EXPECT_EQ(closest->planet, Planet::Earth);
}

TEST_F(CloseApproachTest, Ceres_NoEarthEncounter)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	auto ceres = makeCeres();
	ASSERT_TRUE(ceres);

	CloseApproachConfig cfg;
	cfg.horizonYears = 50.0;
	cfg.thresholdAu = 0.05;
	cfg.sampleStepDay = 1.0;
	cfg.planet = Planet::Earth;

	const auto events = findCloseApproaches(*ceres, m_propagator, cfg);
	EXPECT_TRUE(events.empty()) << "Ceres should have no Earth approaches < 0.05 AU; got " << events.size()
	                            << " events";
}

TEST_F(CloseApproachTest, OutboundStartInsideThresholdIsNotAnEvent)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	const Vec3D r(1.80270860448509540e-1, 1.08308230387968929e+0, 5.79029086249183231e-3);
	const Vec3D vPhysical(-1.59033209065004774e-2, 4.61757822571171473e-3, 5.85213757466040150e-5);
	auto orbit = make_orbit_from_states_vectors({r, vPhysical * (1.0 / Constants::GAUSS_K)}, 2459200.5);
	ASSERT_TRUE(orbit);

	CloseApproachConfig cfg;
	cfg.horizonYears = 0.1;
	cfg.thresholdAu = 0.2;
	cfg.sampleStepDay = 1.0;
	cfg.planet = Planet::Earth;

	const auto events = findCloseApproaches(*orbit, m_propagator, cfg);
	EXPECT_TRUE(events.empty());
}

TEST_F(CloseApproachTest, Defensive_BadConfigReturnsEmpty)
{
	auto apophis = makeApophis();
	ASSERT_TRUE(apophis);

	{
		CloseApproachConfig cfg;
		cfg.horizonYears = 0.0;
		const auto events = findCloseApproaches(*apophis, m_propagator, cfg);
		EXPECT_TRUE(events.empty());
	}

	{
		CloseApproachConfig cfg;
		cfg.sampleStepDay = 0.0;
		const auto events = findCloseApproaches(*apophis, m_propagator, cfg);
		EXPECT_TRUE(events.empty());
	}

	{
		auto* saved = getDefaultEphemerisSource();
		setDefaultEphemerisSource(nullptr);
		CloseApproachConfig cfg;
		const auto events = findCloseApproaches(*apophis, m_propagator, cfg, nullptr);
		EXPECT_TRUE(events.empty());
		setDefaultEphemerisSource(saved);
	}
}

TEST_F(CloseApproachTest, SurfaceMonitoredMissPreservesGrid)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";
	const auto orbit = makeApophis();
	ASSERT_TRUE(orbit);
	const double epoch = orbit->getEpoch();
	const std::vector<double> times{epoch, epoch + 1.0, epoch + 2.0};
	const Vec3D velocity = orbit->getV() * Constants::GAUSS_K;
	const auto reference = m_propagator.propagateAllWithVelocity(orbit->getR(), velocity, epoch, times);
	std::vector<NBodyPropagator::StateVel> samples;
	const auto terminal = m_propagator.propagateUntilSurface(
	    orbit->getR(), velocity, epoch, times, samples, EphemerisPlanet::Earth, planetRadiusKm(Planet::Earth), *m_de440,
	    m_propagator.nonGrav(), NBodyPropagator::SurfaceShape::Sphere);
	EXPECT_FALSE(terminal.hitSurface);
	EXPECT_DOUBLE_EQ(terminal.jd, times.back());
	ASSERT_EQ(samples.size(), reference.size());
	for (std::size_t i = 0; i < samples.size(); ++i) {
		EXPECT_DOUBLE_EQ((samples[i].R - reference[i].R).magnitude(), 0.0);
		EXPECT_DOUBLE_EQ((samples[i].V - reference[i].V).magnitude(), 0.0);
	}
}

TEST_F(CloseApproachTest, AA2014SameVaImpactPericentersMatchReference)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	struct Case {
		int index;
		Vec3D r;
		Vec3D v;
		double pericenterKm;
	};

	const Case cases[] = {{2358,
	                       {-0.18023364170396783, 0.96923362054015716, -0.00044597763318709587},
	                       {-0.017572464202686478, -0.0061273995052632665, 0.00046249666221673453},
	                       621.744},
	                      {3504,
	                       {-0.18021463971838217, 0.96939218859832244, -0.00047210301668478607},
	                       {-0.017595728050108658, -0.0064245762240321838, 0.0005102535943700547},
	                       285.612},
	                      {3505,
	                       {-0.18021461679735187, 0.96939237986297888, -0.00047213452917804137},
	                       {-0.017595752903954421, -0.0064249079796618591, 0.00051030679897850215},
	                       285.253},
	                      {3506,
	                       {-0.18021459385513716, 0.96939257130439727, -0.00047216607079450843},
	                       {-0.017595777772334304, -0.0064252399716607303, 0.00051036004118444871},
	                       284.894},
	                      {3507,
	                       {-0.1802145708919691, 0.96939276292063126, -0.00047219764121354184},
	                       {-0.017595802655433274, -0.0064255722003075978, 0.00051041332104764391},
	                       284.550},
	                      {4000,
	                       {-0.18019921713867235, 0.96952087895479244, -0.0004933058955252981},
	                       {-0.017610885896776601, -0.0066347675099103984, 0.00054390500637628125},
	                       90.148},
	                      {4422,
	                       {-0.18015622264398037, 0.9698796140469379, -0.0005524109437415654},
	                       {-0.017642740930680401, -0.0071341246773716609, 0.00062344537137972225},
	                       261.647},
	                      {4423,
	                       {-0.18015589353299519, 0.96988235997782668, -0.00055286336316625868},
	                       {-0.017642946091268787, -0.0071376250493113882, 0.00062400115777709586},
	                       268.513},
	                      {4424,
	                       {-0.18015555933216701, 0.96988514837468276, -0.00055332277929312501},
	                       {-0.01764315396527626, -0.0071411757302944729, 0.00062456490922240795},
	                       275.589}};
	constexpr double epoch = 2456658.7874675924;
	constexpr double auKm = Constants::AU_IN_METERS / 1000.0;
	m_propagator.setUseSb441(true);
	CloseApproachConfig cfg;
	cfg.horizonYears = 0.008;
	cfg.thresholdAu = 0.2;
	cfg.sampleStepDay = 1.0;
	for (const auto& c : cases) {
		SCOPED_TRACE(c.index);
		const auto events = findCloseApproachesFromState(c.r, c.v, epoch, m_propagator, cfg, m_de440.get());
		ASSERT_EQ(events.size(), 1u);
		const auto& event = events.front();
		EXPECT_NEAR(event.distAu * auKm, c.pericenterKm, 0.2);
		const Vec3D dr = event.rTarget - event.rPlanet;
		const Vec3D dv = event.vTarget - event.vPlanet;
		EXPECT_NEAR(dr.dot(dv) / dr.magnitude(), 0.0, 1e-10);
		EXPECT_NEAR(dr.magnitude(), event.distAu, 1e-16);

		CloseApproachTargetPlaneDifferentialInput input;
		input.event = &event;
		input.parameters.resize(6);
		input.parameters << c.r.x, c.r.y, c.r.z, c.v.x, c.v.y, c.v.z;
		input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6) * 1e-14;
		input.parameterStepPerSigma = Eigen::VectorXd::Zero(6);
		input.parameterStepPerSigma(0) = 1e-7;
		input.nominalNonGrav = &m_propagator.nonGrav();
		input.propagator = &m_propagator;
		input.eph = m_de440.get();
		input.closeApproach = cfg;
		input.planeKind = EncounterPlaneKind::OpikTargetPlane;
		input.epochTdb = epoch;
		const auto trace = makeCloseApproachTargetPlaneTracePoint(input, {});
		ASSERT_TRUE(trace);
		EXPECT_DOUBLE_EQ(trace->event.jd, trace->trace.jd);
		EXPECT_NEAR(trace->event.distAu * auKm, event.distAu * auKm, 0.01);
		EXPECT_NEAR(trace->event.jd, event.jd, 1e-7);
		EXPECT_DOUBLE_EQ(trace->event.vInfKmS, trace->vInfKmS);
		EXPECT_NEAR(event.vInfKmS, trace->vInfKmS, 0.001);
		EXPECT_GT(event.vInfKmS, 4.0);
		EXPECT_LT(event.vInfKmS, 9.0);
	}
}

TEST_F(CloseApproachTest, MonitoringKeepsEntryStepCapAcrossMinimumOutputs)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";

	class BoundedEphemeris final : public IEphemerisSource {
	   public:
		explicit BoundedEphemeris(IEphemerisSource& source) : source(source) {}

		EphemerisRecord queryJD(EphemerisPlanet planet, double jd) override
		{
			if (planet == EphemerisPlanet::Earth && ++earthQueries > 30000)
				throw std::runtime_error("Encounter step cap was reset at a minimum");
			return source.queryJD(planet, jd);
		}

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet planet, double jd, int days) override
		{
			return source.fetchWindow(planet, jd, days);
		}

		bool isAvailable() const override { return source.isAvailable(); }

		std::string name() const override { return "bounded-encounter-regression"; }

		IEphemerisSource& source;
		std::size_t earthQueries = 0;
	} bounded(*m_de440);

	auto propagator = m_propagator.forEncounterMonitoring(bounded, 600, .2);
	const std::vector<double> times{2463886.8451113841, 2467899.1643204414, 2468151.9443960623};
	const double distances[]{.006316220789579, .000048045660579, .114899190270227};
	std::vector<CloseApproachEvent> events;
	for (std::size_t i = 0; i < times.size(); ++i) {
		CloseApproachEvent event;
		event.planet = Planet::Earth;
		event.jd = times[i];
		event.distAu = distances[i];
		const auto earth = m_de440->queryJD(EphemerisPlanet::Earth, event.jd);
		event.rPlanet = earth.R;
		event.vPlanet = earth.V;
		events.push_back(event);
	}
	CloseApproachTargetPlaneDifferentialInput input;
	input.parameters.resize(6);
	input.parameters << .62433518638649632, .83391558445757552, .0034703537785752041, -.01285383170372115,
	    .011750243759764011, .000093034442710722442;
	input.parameterCovariance = Eigen::MatrixXd::Identity(6, 6) * 1e-16;
	input.parameterStepPerSigma = Eigen::VectorXd::Constant(6, 1e-8);
	input.nominalNonGrav = &propagator.nonGrav();
	input.propagator = &propagator;
	input.eph = &bounded;
	input.event = &events.front();
	input.epochTdb = 2459170.2715327665;
	input.planeKind = EncounterPlaneKind::OpikTargetPlane;
	for (const bool fitA2 : {false, true}) {
		SCOPED_TRACE(fitA2);
		bounded.earthQueries = 0;
		if (fitA2) {
			propagator.nonGrav().setModel(NgModel::Yarkovsky);
			input.parameters.conservativeResize(7);
			input.parameters(6) = 0.0;
			input.parameterCovariance = Eigen::MatrixXd::Identity(7, 7) * 1e-16;
			input.parameterStepPerSigma = Eigen::VectorXd::Constant(7, 1e-8);
		}
		const auto traces = makeCloseApproachTargetPlaneTracePoints(input, events);
		ASSERT_EQ(traces.size(), times.size());
		for (const auto& trace : traces) {
			ASSERT_TRUE(trace);
			EXPECT_EQ(trace->dCoordsDParams.cols(), fitA2 ? 7 : 6);
			EXPECT_TRUE(trace->dCoordsDParams.array().isFinite().all());
		}
		EXPECT_LT(bounded.earthQueries, 30000U);
	}
}

TEST_F(CloseApproachTest, AA2014SurfaceStopKeepsOnlyCompletedPrefix)
{
	if (!hasDE440)
		GTEST_SKIP() << "DE440 kernel not present";
	constexpr double epoch = 2456658.7874675924;
	const Vec3D r(-0.18023364170396783, 0.96923362054015716, -0.00044597763318709587);
	const Vec3D v(-0.017572464202686478, -0.0061273995052632665, 0.00046249666221673453);
	const std::vector<double> times{epoch, epoch + 0.5, epoch + 1.0, epoch + 2.0};
	std::vector<NBodyPropagator::StateVel> samples;
	const auto terminal = m_propagator.propagateUntilSurface(
	    r, v, epoch, times, samples, EphemerisPlanet::Earth, planetRadiusKm(Planet::Earth), *m_de440,
	    m_propagator.nonGrav(), NBodyPropagator::SurfaceShape::Sphere);
	ASSERT_TRUE(terminal.hitSurface);
	ASSERT_EQ(samples.size(), 2u);
	EXPECT_GT(terminal.jd, times[1]);
	EXPECT_LT(terminal.jd, times[2]);
	const double distanceKm = (terminal.R - m_de440->queryJD(EphemerisPlanet::Earth, terminal.jd).R).magnitude() *
	                          Constants::AU_IN_METERS / 1000.0;
	EXPECT_NEAR(distanceKm, planetRadiusKm(Planet::Earth), 0.02);

	CloseApproachConfig cfg;
	cfg.horizonYears = 0.008;
	cfg.lookbackYears = 0.001;
	const auto events = findCloseApproachesFromState(r, v, epoch, m_propagator, cfg);
	ASSERT_EQ(events.size(), 1u);
	EXPECT_GT(events.front().jd, epoch);
	EXPECT_NEAR(events.front().distAu * Constants::AU_IN_METERS / 1000.0, 621.744, 0.2);
}
