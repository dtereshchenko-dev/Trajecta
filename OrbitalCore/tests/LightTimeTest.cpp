#include <gtest/gtest.h>
#include "obs/LightTime.h"
#include "orbital_core/IEphemerisSource.h"
#include "orbital_core/Types.h"

#include <optional>
#include <string>
#include <vector>

namespace {
	using orbital_core::obs::solveLightTimeTwoPass;

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

	constexpr double kReceptionJd = 2400000.5;
	const Vec3D kObserver(1.0, 0.0, 0.0);
	const Vec3D kTarget(1.2, 1.5, 0.1);
	const Vec3D kTargetVelocity(-0.004, 0.009, 0.0005);

	struct ReferenceSolution {
		double delayDays;
		Vec3D position;
	};

	ReferenceSolution barycentricAber2(const Vec3D& sunVelocity)
	{
		const StateVectors reception(kTarget, kTargetVelocity / Constants::GAUSS_K);
		const double first = (kTarget - kObserver).magnitude() / Constants::AU_PER_DAY;
		const Vec3D firstEmission = Calculations::propagateKeplerState(reception, -first).first + sunVelocity * -first;
		const double refined = (firstEmission - kObserver).magnitude() / Constants::AU_PER_DAY;
		return {refined, Calculations::propagateKeplerState(reception, -refined).first + sunVelocity * -refined};
	}

	void expectNear(const Vec3D& actual, const Vec3D& expected, double tolerance)
	{
		EXPECT_NEAR(actual.x, expected.x, tolerance);
		EXPECT_NEAR(actual.y, expected.y, tolerance);
		EXPECT_NEAR(actual.z, expected.z, tolerance);
	}
}  // namespace

TEST(LightTime, WithoutBarycentricSunKeepsHeliocentricAber2)
{
	ScopedEphemeris none(nullptr);
	const auto solution = solveLightTimeTwoPass(kReceptionJd, kTarget, kTargetVelocity, kObserver);
	ASSERT_TRUE(solution.valid);
	const auto reference = barycentricAber2(Vec3D());
	expectNear(solution.position, reference.position, 1e-15);
	EXPECT_NEAR(kReceptionJd - solution.emissionJd, reference.delayDays, 1e-9);
}

TEST(LightTime, EmissionPointIncludesTheSunsBarycentricMotion)
{
	const Vec3D sunVelocity(4.0e-6, -7.0e-6, 1.0e-7);
	UniformlyMovingSun sun(Vec3D(0.003, -0.004, 0.0001), sunVelocity);
	ScopedEphemeris scoped(&sun);
	const auto solution = solveLightTimeTwoPass(kReceptionJd, kTarget, kTargetVelocity, kObserver);
	ASSERT_TRUE(solution.valid);
	const auto reference = barycentricAber2(sunVelocity);
	expectNear(solution.position, reference.position, 1e-15);

	const Vec3D heliocentric = barycentricAber2(Vec3D()).position;
	const Vec3D shift = solution.position - heliocentric;
	const Vec3D expectedShift = sunVelocity * -reference.delayDays;
	EXPECT_LT((shift - expectedShift).magnitude(), 1e-3 * expectedShift.magnitude());
}

TEST(LightTime, ReceptionStateJacobianMatchesIndependentPerturbations)
{
	UniformlyMovingSun sun({0.003, -0.004, 0.0001}, {4e-6, -7e-6, 1e-7});
	for (auto* source : {static_cast<IEphemerisSource*>(nullptr), static_cast<IEphemerisSource*>(&sun)}) {
		ScopedEphemeris scoped(source);
		for (const Vec3D r : {Vec3D(1.001, .003, .001), kTarget, Vec3D(4, -2, .1)}) {
			for (const Vec3D v : {kTargetVelocity, Vec3D(-.06, .03, .01)}) {
				Eigen::Matrix<double, 3, 6> jacobian;
				const auto solution = solveLightTimeTwoPass(kReceptionJd, r, v, kObserver, &jacobian);
				ASSERT_TRUE(solution.valid);
				const auto nominal = solveLightTimeTwoPass(kReceptionJd, r, v, kObserver);
				expectNear(solution.position, nominal.position, 0);
				for (int column = 0; column < 6; ++column) {
					SCOPED_TRACE(column);
					const double h = column < 3 ? 1e-7 : 1e-5;
					const Vec3D delta(column % 3 == 0 ? h : 0, column % 3 == 1 ? h : 0, column % 3 == 2 ? h : 0);
					const Vec3D dr = column < 3 ? delta : Vec3D{}, dv = column >= 3 ? delta : Vec3D{};
					const auto plus = solveLightTimeTwoPass(kReceptionJd, r + dr, v + dv, kObserver);
					const auto minus = solveLightTimeTwoPass(kReceptionJd, r - dr, v - dv, kObserver);
					ASSERT_TRUE(plus.valid && minus.valid);
					const Vec3D derivative = (plus.position - minus.position) / (2 * h);
					expectNear({jacobian(0, column), jacobian(1, column), jacobian(2, column)}, derivative, 1e-8);
				}
			}
		}
	}
}
