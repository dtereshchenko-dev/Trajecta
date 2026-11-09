#include <gtest/gtest.h>
#include "obs/LightTime.h"
#include "obs/LightTimeSeries.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Types.h"

#include <optional>
#include <string>
#include <vector>

namespace {
	class UniformlyMovingSun final : public IEphemerisSource {
	   public:
		UniformlyMovingSun(const Vec3D& origin, const Vec3D& velocity) : origin(origin), velocity(velocity) {}

		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, {}, {}}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "UniformlyMovingSun"; }

		std::optional<EphemerisRecord> barycentricSunMjd(double mjd) override
		{
			return EphemerisRecord{2400000.5 + mjd, origin + velocity * mjd, velocity};
		}

		Vec3D origin;
		Vec3D velocity;
	};

	class ScopedEphemeris {
	   public:
		explicit ScopedEphemeris(IEphemerisSource* source) : previous(getDefaultEphemerisSource())
		{
			setDefaultEphemerisSource(source);
		}

		~ScopedEphemeris() { setDefaultEphemerisSource(previous); }

		IEphemerisSource* previous;
	};

	const Vec3D kObserver(1.0, 0.0, 0.0);

	void expectNear(const Vec3D& actual, const Vec3D& expected, double tolerance)
	{
		EXPECT_NEAR(actual.x, expected.x, tolerance);
		EXPECT_NEAR(actual.y, expected.y, tolerance);
		EXPECT_NEAR(actual.z, expected.z, tolerance);
	}
}  // namespace

TEST(LightTime, BatchedEmissionSatisfiesLightTimeAndLeavesRadarAtReception)
{
	using namespace orbital_core;
	UniformlyMovingSun sun({0.003, -0.004, 0.0001}, {4e-6, -7e-6, 1e-7});
	ScopedEphemeris scoped(&sun);
	constexpr double epoch = 60000.0;
	const Vec3D position(60.340206, -943.345623, -517.355543);
	const Vec3D velocity(0.001080, -0.016335, -0.008933);
	const std::vector<double> times{epoch, epoch};
	const std::vector<NBodyPropagator::StateVel> reception{{position, velocity}, {position, velocity}};
	const std::vector<Vec3D> observers{kObserver, kObserver};
	const NBodyPropagator propagator;
	const NonGravModel nonGrav;
	const auto result = obs::solveLightTimeSeries(position, velocity, epoch, times, reception, observers,
	                                              {true, false}, false, propagator, nonGrav);
	ASSERT_EQ(result.states.size(), 2u);
	const Vec3D apparent = obs::apparentEmissionPosition(result.states[0], result.emissionMjd[0], epoch,
	                                                    result.emissionOffsetDays[0]);
	const double delay = (epoch - result.emissionMjd[0]) - result.emissionOffsetDays[0];
	EXPECT_NEAR(delay, (apparent - kObserver).magnitude() / Constants::AU_PER_DAY, 1e-11);
	EXPECT_DOUBLE_EQ(result.emissionMjd[1], epoch);
	EXPECT_DOUBLE_EQ(result.emissionOffsetDays[1], 0.0);
	expectNear(result.states[1].R, position, 1e-14);
	const auto dense = propagator.propagateDenseMjd(position, velocity, epoch, epoch - 10.0, epoch + 1.0,
	                                                 nonGrav, true, true);
	const auto interpolated = obs::solveDenseLightTimeSeries([&](double t) { return dense.at(t); }, times, observers,
	                                                         {true, false});
	expectNear(interpolated.light.states[0].R, result.states[0].R, 1e-10);
	expectNear(interpolated.light.states[1].R, position, 1e-14);
	EXPECT_DOUBLE_EQ(interpolated.light.emissionMjd[1], epoch);
}

TEST(LightTime, DenseStateAndVariationsMatchIndependentPropagations)
{
	ScopedEphemeris none(nullptr);
	constexpr double epoch = 60000.0;
	const Vec3D position(1.2, 0.4, 0.1);
	const Vec3D velocity(-0.005, 0.013, 0.001);
	const std::vector<double> times{epoch - 60.0, epoch - 17.3, epoch, epoch + 3.5, epoch + 55.0};
	NonGravModel ng;
	ng.setModel(NgModel::InverseSquare);
	ng.setEnabled(true);
	ng.setParams(1e-9, -2e-10, 3e-11);
	for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
		NBodyPropagator propagator(kind);
		propagator.setPlanetMask(0);
		auto dense = propagator.propagateDenseMjd(position, velocity, epoch, epoch - 10.0, times.back(), ng,
		                                          true, true);
		EXPECT_FALSE(dense.at(times.front()));
		propagator.extendDenseBackwardMjd(dense, times.front(), ng, true, true);
		const auto reference = propagator.propagateAllWithStmAndSensMjd(position, velocity, epoch, times, ng);
		ASSERT_EQ(reference.size(), times.size());
		for (std::size_t i = 0; i < times.size(); ++i) {
			const auto state = dense.at(times[i]);
			ASSERT_TRUE(state) << "integrator=" << static_cast<int>(kind) << " epoch=" << times[i];
			expectNear(state->R, reference[i].R, 1e-10);
			expectNear(state->V, reference[i].V, 1e-11);
			EXPECT_TRUE(state->Phi.isApprox(reference[i].Phi, 1e-9));
			EXPECT_TRUE(state->Sens.isApprox(reference[i].SensA, 1e-8));
		}
		EXPECT_FALSE(dense.at(epoch - 1000.0));
	}
}

TEST(LightTime, DenseSeriesMatchesConvergedFullPropagation)
{
	using namespace orbital_core;
	ScopedEphemeris none(nullptr);
	constexpr double epoch = 60000.0;
	const Vec3D position(1.2, 0.4, 0.1);
	const Vec3D velocity(-0.005, 0.013, 0.001);
	const std::vector<double> times{epoch - 60.0, epoch - 17.3, epoch, epoch + 3.5, epoch + 55.0};
	const std::vector<Vec3D> observers(times.size(), kObserver);
	const std::vector<bool> optical{true, true, false, true, true};
	NonGravModel ng;
	ng.setModel(NgModel::InverseSquare);
	ng.setEnabled(true);
	ng.setParams(1e-9, -2e-10, 3e-11);
	for (const auto kind : {IntegratorKind::RA15, IntegratorKind::IAS15}) {
		NBodyPropagator p(kind);
		p.setPlanetMask(0);
		const auto reception = p.propagateAllWithVelocityMjd(position, velocity, epoch, times, ng);
		const auto reference = obs::solveLightTimeSeries(position, velocity, epoch, times, reception, observers,
		                                                  optical, true, p, ng);
		const auto dense = obs::propagateDenseLightTimeSeries(position, velocity, epoch, times, observers, optical,
		                                                      p, ng, false, false);
		ASSERT_EQ(dense.light.states.size(), reference.states.size());
		for (std::size_t i = 0; i < times.size(); ++i) {
			expectNear(dense.light.states[i].R, reference.states[i].R, 1e-10);
			expectNear(dense.light.states[i].V, reference.states[i].V, 1e-11);
			EXPECT_NEAR(dense.light.emissionMjd[i], reference.emissionMjd[i], 1e-10);
		}
	}
}

TEST(LightTime, EmissionPointIncludesTheSunsBarycentricMotion)
{
	const Vec3D sunVelocity(4.0e-6, -7.0e-6, 1.0e-7);
	UniformlyMovingSun sun(Vec3D(0.003, -0.004, 0.0001), sunVelocity);
	ScopedEphemeris scoped(&sun);
	const double delay = 0.125;
	const auto shift = orbital_core::obs::solarDisplacementDuringDelay(2400000.5, delay);
	expectNear(shift, sunVelocity * -delay, 1e-15);
}

TEST(LightTime, DistantHyperbolicStateAgreesWithNumericalIntegration)
{
	ScopedEphemeris none(nullptr);
	const Vec3D position(60.340206, -943.345623, -517.355543);
	const Vec3D velocity(0.001080, -0.016335, -0.008933);
	const auto propagate = [&](double dt) {
		Vec3D r = position, v = velocity;
		const double h = dt / 128.0;
		const auto acceleration = [](const Vec3D& p) { return p * (-Constants::GM_SUN / std::pow(p.magnitude(), 3)); };
		for (int i = 0; i < 128; ++i) {
			const Vec3D k1r = v, k1v = acceleration(r);
			const Vec3D k2r = v + k1v * (h / 2), k2v = acceleration(r + k1r * (h / 2));
			const Vec3D k3r = v + k2v * (h / 2), k3v = acceleration(r + k2r * (h / 2));
			const Vec3D k4r = v + k3v * h, k4v = acceleration(r + k3r * h);
			r += (k1r + k2r * 2 + k3r * 2 + k4r) * (h / 6);
			v += (k1v + k2v * 2 + k3v * 2 + k4v) * (h / 6);
		}
		return std::pair{r, v};
	};
	const double delay = (position - kObserver).magnitude() / Constants::AU_PER_DAY;
	const auto reference = propagate(-delay);
	const auto state = Calculations::propagateKeplerState({position, velocity / Constants::GAUSS_K}, -delay);
	ASSERT_TRUE(state.isValid());
	expectNear(state.first, reference.first, 2e-10);
	expectNear(state.second * Constants::GAUSS_K, reference.second, 2e-13);
}
