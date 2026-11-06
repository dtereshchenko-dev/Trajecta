#include <gtest/gtest.h>
#include "earth_orientation/EarthRotation.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/Constants.h"

#include <cmath>

namespace {
	constexpr double ASEC_TO_RAD = 4.84813681109535994e-6;
	constexpr double J2000 = Constants::J2000_JD;

	inline double dot(const Vec3D& a, const Vec3D& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

	inline double norm(const Vec3D& v) { return std::sqrt(dot(v, v)); }
}  // namespace

TEST(Nutation, MeanObliquityAtJ2000)
{
	const double eps = EarthRotation::Iau2000B::meanObliquity(J2000);
	const double expected = 84381.406 * ASEC_TO_RAD;
	EXPECT_NEAR(eps, expected, 1e-12);
}

TEST(Nutation, MeanObliquityDecreasesWithTime)
{
	const double eps0 = EarthRotation::Iau2000B::meanObliquity(J2000);
	const double eps100 = EarthRotation::Iau2000B::meanObliquity(J2000 + 36525.0);
	const double delta_asec = (eps0 - eps100) / ASEC_TO_RAD;
	EXPECT_NEAR(delta_asec, 46.836769, 0.01);
}

TEST(Nutation, AnglesMagnitudeBoundsAtJ2000)
{
	double dpsi = 0.0, deps = 0.0;
	EarthRotation::Iau2000B::nutationAngles(J2000, dpsi, deps);

	EXPECT_LT(std::fabs(dpsi), 20.0 * ASEC_TO_RAD);
	EXPECT_LT(std::fabs(deps), 10.0 * ASEC_TO_RAD);

	EXPECT_GT(std::fabs(dpsi), 1.0 * ASEC_TO_RAD);
	EXPECT_GT(std::fabs(deps), 1.0 * ASEC_TO_RAD);
}

TEST(Nutation, AnglesAcrossACentury)
{
	for (double years : {-100.0, -50.0, -10.0, 0.0, 10.0, 50.0, 100.0}) {
		const double jd = J2000 + years * 365.25;
		double dpsi = 0.0, deps = 0.0;
		EarthRotation::Iau2000B::nutationAngles(jd, dpsi, deps);
		EXPECT_LT(std::fabs(dpsi), 25.0 * ASEC_TO_RAD) << "year offset " << years;
		EXPECT_LT(std::fabs(deps), 15.0 * ASEC_TO_RAD) << "year offset " << years;
	}
}

TEST(Nutation, PrincipalPeriodDominatesAt18_6Year)
{
	double dpsi0 = 0.0, deps0 = 0.0;
	EarthRotation::Iau2000B::nutationAngles(J2000, dpsi0, deps0);

	double dpsi1 = 0.0, deps1 = 0.0;
	const double jd_half = J2000 + 9.3 * 365.25;
	EarthRotation::Iau2000B::nutationAngles(jd_half, dpsi1, deps1);

	EXPECT_GT(std::fabs(dpsi0 - dpsi1), 10.0 * ASEC_TO_RAD);
}

TEST(Nutation, EquationOfEquinoxesIsDpsiCosEps)
{
	const double jd = J2000 + 5.0 * 365.25;
	double dpsi = 0.0, deps = 0.0;
	EarthRotation::Iau2000B::nutationAngles(jd, dpsi, deps);
	const double eps = EarthRotation::Iau2000B::meanObliquity(jd);
	const double ee = EarthRotation::Iau2000B::equationOfEquinoxes(jd);
	EXPECT_NEAR(ee, dpsi * std::cos(eps), 1e-15);
}

TEST(Nutation, ApplyNutationPreservesMagnitude)
{
	Vec3D v(0.7, -0.3, 0.5);
	const double mag0 = norm(v);
	EarthRotation::Iau2000B::applyNutation(J2000 + 1000.0, v, +1);
	EXPECT_NEAR(norm(v), mag0, 1e-14);
}

TEST(Nutation, ApplyNutationRoundtrip)
{
	for (double years : {-50.0, 0.0, 5.0, 50.0}) {
		const double jd = J2000 + years * 365.25;
		Vec3D v(1.0, 2.0, 3.0);
		Vec3D v0 = v;
		EarthRotation::Iau2000B::applyNutation(jd, v, +1);
		EarthRotation::Iau2000B::applyNutation(jd, v, -1);
		EXPECT_NEAR(v.x, v0.x, 1e-13) << "year " << years;
		EXPECT_NEAR(v.y, v0.y, 1e-13) << "year " << years;
		EXPECT_NEAR(v.z, v0.z, 1e-13) << "year " << years;
	}
}

TEST(Nutation, ApplyPrecessionNutationRoundtrip)
{
	for (double years : {-50.0, 0.0, 5.0, 50.0}) {
		const double jd = J2000 + years * 365.25;
		Vec3D v(0.5, -0.2, 0.8);
		Vec3D v0 = v;
		EarthRotation::Iau2000B::applyPrecessionNutation(jd, v, +1);
		EarthRotation::Iau2000B::applyPrecessionNutation(jd, v, -1);
		EXPECT_NEAR(v.x, v0.x, 1e-13) << "year " << years;
		EXPECT_NEAR(v.y, v0.y, 1e-13) << "year " << years;
		EXPECT_NEAR(v.z, v0.z, 1e-13) << "year " << years;
	}
}

TEST(Nutation, ApplyPrecessionNutationAtJ2000IsSmallRotation)
{
	Vec3D v(1.0, 0.0, 0.0);
	Vec3D v0 = v;
	EarthRotation::Iau2000B::applyPrecessionNutation(J2000, v, +1);

	const double dx = v.x - v0.x;
	const double dy = v.y - v0.y;
	const double dz = v.z - v0.z;
	const double disp = std::sqrt(dx * dx + dy * dy + dz * dz);

	EXPECT_LT(disp, 1.2e-4);
	EXPECT_GT(disp, 1.0e-6);
	EXPECT_NEAR(norm(v), 1.0, 1e-14);
}

TEST(Nutation, ForwardAndReverseAreInverses)
{
	const double jd = J2000 + 1234.5;
	Vec3D a(1.7, -0.9, 0.4);
	Vec3D b = a;
	EarthRotation::Iau2000B::applyNutation(jd, b, -1);
	EarthRotation::Iau2000B::applyNutation(jd, b, +1);
	EXPECT_NEAR(b.x, a.x, 1e-13);
	EXPECT_NEAR(b.y, a.y, 1e-13);
	EXPECT_NEAR(b.z, a.z, 1e-13);
}

TEST(NutationReferenceModel, MatchesObserverTraceAtFirst2020VVObservation)
{
	constexpr double jdTt = 2459139.71531574056;
	constexpr double jdUt1 = 2459139.71451301407;

	EXPECT_NEAR(EarthRotation::Iau1976_1980::gmstRadians(jdUt1), 1.80459500639612003, 2e-9);
	EXPECT_NEAR(EarthRotation::Iau1976_1980::equationOfEquinoxes(jdTt), -8.17278227348924873e-5, 1e-11);

	Vec3D v(2.40975897709992899e-5, -3.18222036023614867e-5, 1.49881941172809766e-5);
	EarthRotation::Iau1976_1980::applyPrecessionNutation(jdTt, v, -1);

	constexpr double eclX = 2.39816710640224586e-5;
	constexpr double eclY = -2.33538647314036919e-5;
	constexpr double eclZ = 2.64096363637052111e-5;
	const double cosE = std::cos(Constants::EARTH_INCL);
	const double sinE = std::sin(Constants::EARTH_INCL);
	EXPECT_NEAR(v.x, eclX, 2e-14);
	EXPECT_NEAR(v.y, eclY * cosE - eclZ * sinE, 2e-14);
	EXPECT_NEAR(v.z, eclY * sinE + eclZ * cosE, 2e-14);
}

TEST(NutationReferenceModel, PrecessionNutationRoundtrip)
{
	for (double years : {-50.0, 0.0, 20.0, 50.0}) {
		const double jd = J2000 + years * 365.25;
		Vec3D v(0.5, -0.2, 0.8);
		const Vec3D original = v;
		EarthRotation::Iau1976_1980::applyPrecessionNutation(jd, v, +1);
		EarthRotation::Iau1976_1980::applyPrecessionNutation(jd, v, -1);
		EXPECT_NEAR(v.x, original.x, 1e-13) << "year " << years;
		EXPECT_NEAR(v.y, original.y, 1e-13) << "year " << years;
		EXPECT_NEAR(v.z, original.z, 1e-13) << "year " << years;
	}
}
