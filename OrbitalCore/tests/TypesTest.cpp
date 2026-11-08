#include <gtest/gtest.h>
#include "orbital_core/Types.h"

TEST(Convertors, DegreesToRadians)
{
	EXPECT_NEAR(Convertors::degreesToRadians(0.0), 0.0, 1e-15);
	EXPECT_NEAR(Convertors::degreesToRadians(90.0), Constants::PI / 2.0, 1e-15);
	EXPECT_NEAR(Convertors::degreesToRadians(180.0), Constants::PI, 1e-15);
	EXPECT_NEAR(Convertors::degreesToRadians(360.0), 2.0 * Constants::PI, 1e-15);
	EXPECT_NEAR(Convertors::degreesToRadians(-45.0), -Constants::PI / 4.0, 1e-15);
}

TEST(Convertors, RadiansToDegrees)
{
	EXPECT_NEAR(Convertors::radiansToDegrees(0.0), 0.0, 1e-12);
	EXPECT_NEAR(Convertors::radiansToDegrees(Constants::PI / 2.0), 90.0, 1e-12);
	EXPECT_NEAR(Convertors::radiansToDegrees(Constants::PI), 180.0, 1e-12);
	EXPECT_NEAR(Convertors::radiansToDegrees(2.0 * Constants::PI), 360.0, 1e-12);
}

TEST(Convertors, DegreesRadiansRoundtrip)
{
	for (double deg : {0.0, 30.0, 45.0, 90.0, 123.456, -270.0, 359.99}) {
		EXPECT_NEAR(Convertors::radiansToDegrees(Convertors::degreesToRadians(deg)), deg, 1e-10);
	}
}

TEST(Convertors, ConvertKmAURoundtrip)
{
	EXPECT_NEAR(Convertors::convertAUToKm(1.0), 149597870.7, 1e-1);
	EXPECT_NEAR(Convertors::convertKmToAU(149597870.7), 1.0, 1e-12);

	for (double au : {0.5, 1.0, 5.2, 30.0}) {
		EXPECT_NEAR(Convertors::convertKmToAU(Convertors::convertAUToKm(au)), au, 1e-12);
	}
}

TEST(Convertors, ConvertGregorianToJD)
{
	EXPECT_NEAR(Convertors::convertGregorianToJD(2000, 1, 1.5), Constants::J2000_JD, 1e-6);

	EXPECT_NEAR(Convertors::convertGregorianToJD(1970, 1, 1.0), 2440587.5, 1e-6);

	EXPECT_NEAR(Convertors::convertGregorianToJD(2019, 5, 31.0), 2458634.5, 1e-6);
}

TEST(Convertors, ConvertJDToGregorian)
{
	const auto start = Convertors::convertGregorianToJD(2024, 1, 1.0);
	const auto atSeconds = [start](double seconds) { return Convertors::convertJDToGregorian(start + seconds / 86400.0); };
	EXPECT_EQ(atSeconds(21 * 3600 + 50 * 60 + 1.54), "2024-01-01 21:50:01.54");
	EXPECT_EQ(atSeconds(21 * 3600 + 50 * 60 + 0.01), "2024-01-01 21:50:00.01");
	EXPECT_EQ(atSeconds(59.999), "2024-01-01 00:01:00.00");
	EXPECT_EQ(atSeconds(86399.999), "2024-01-02 00:00:00.00");
	EXPECT_EQ(Convertors::convertJDToGregorian(Constants::J2000_JD), "2000-01-01 12:00:00.00");
}

TEST(Convertors, StringToInt)
{
	EXPECT_EQ(Convertors::stringToInt("123"), 123);
	EXPECT_EQ(Convertors::stringToInt("-5"), -5);
	EXPECT_EQ(Convertors::stringToInt("abc"), 0);
	EXPECT_EQ(Convertors::stringToInt(""), 0);
}

TEST(Convertors, StringToDouble)
{
	EXPECT_NEAR(Convertors::stringToDouble("3.14"), 3.14, 1e-12);
	EXPECT_NEAR(Convertors::stringToDouble("-1.5"), -1.5, 1e-12);
	EXPECT_NEAR(Convertors::stringToDouble("abc"), 0.0, 1e-12);
	EXPECT_NEAR(Convertors::stringToDouble(""), 0.0, 1e-12);
}

TEST(Convertors, RaToDegrees)
{
	EXPECT_NEAR(Convertors::raToDegrees(0, 0, 0.0), 0.0, 1e-10);
	EXPECT_NEAR(Convertors::raToDegrees(6, 0, 0.0), 90.0, 1e-10);
	EXPECT_NEAR(Convertors::raToDegrees(12, 0, 0.0), 180.0, 1e-10);
	EXPECT_NEAR(Convertors::raToDegrees(24, 0, 0.0), 360.0, 1e-10);
	EXPECT_NEAR(Convertors::raToDegrees(1, 30, 0.0), 22.5, 1e-10);
}

TEST(Convertors, RaToRadians)
{
	EXPECT_NEAR(Convertors::raToRadians(6, 0, 0.0), Constants::PI / 2.0, 1e-10);
	EXPECT_NEAR(Convertors::raToRadians(12, 0, 0.0), Constants::PI, 1e-10);
}

TEST(Convertors, DecToDegrees)
{
	EXPECT_NEAR(Convertors::decToDegrees(45, 30, 0.0), 45.5, 1e-10);
	EXPECT_NEAR(Convertors::decToDegrees(-45, 30, 0.0), -45.5, 1e-10);
	EXPECT_NEAR(Convertors::decToDegrees(0, 0, 0.0), 0.0, 1e-10);
	EXPECT_NEAR(Convertors::decToDegrees(23, 26, 12.0), 23.43666667, 1e-5);
}

TEST(Convertors, RadiansToHMS)
{
	int h, m;
	double s;
	Convertors::radiansToHMS(Constants::PI / 2.0, h, m, s);
	EXPECT_EQ(h, 6);
	EXPECT_EQ(m, 0);
	EXPECT_NEAR(s, 0.0, 1e-6);

	Convertors::radiansToHMS(Constants::PI, h, m, s);
	EXPECT_EQ(h, 12);
	EXPECT_EQ(m, 0);
	EXPECT_NEAR(s, 0.0, 1e-6);
}

TEST(Convertors, RadiansToDMS)
{
	int d, m;
	double s;
	Convertors::radiansToDMS(Convertors::degreesToRadians(45.5), d, m, s);
	EXPECT_EQ(d, 45);
	EXPECT_EQ(m, 30);
	EXPECT_NEAR(s, 0.0, 1e-4);
}

TEST(Convertors, EcvToEcl)
{
	Vec3D xAxis(1.0, 0.0, 0.0);
	Vec3D result = Convertors::ecvToEcl(xAxis, Constants::EARTH_INCL);
	EXPECT_NEAR(result.x, 1.0, 1e-15);
	EXPECT_NEAR(result.y, 0.0, 1e-15);
	EXPECT_NEAR(result.z, 0.0, 1e-15);

	Vec3D v(0.3, 0.7, -0.5);
	Vec3D ecl = Convertors::ecvToEcl(v, Constants::EARTH_INCL);
	double ce = std::cos(Constants::EARTH_INCL), se = std::sin(Constants::EARTH_INCL);
	Vec3D back(ecl.x, ecl.y * ce - ecl.z * se, ecl.y * se + ecl.z * ce);
	EXPECT_NEAR(back.x, v.x, 1e-14);
	EXPECT_NEAR(back.y, v.y, 1e-14);
	EXPECT_NEAR(back.z, v.z, 1e-14);
}

TEST(Helpers, WrapPi)
{
	EXPECT_NEAR(Helpers::wrapPi(0.0), 0.0, 1e-15);
	EXPECT_NEAR(Helpers::wrapPi(Constants::PI), Constants::PI, 1e-15);
	EXPECT_NEAR(Helpers::wrapPi(-Constants::PI), -Constants::PI, 1e-15);
	EXPECT_NEAR(Helpers::wrapPi(2.0 * Constants::PI), 0.0, 1e-15);
	EXPECT_NEAR(std::abs(Helpers::wrapPi(3.0 * Constants::PI)), Constants::PI, 1e-15);
	EXPECT_NEAR(std::abs(Helpers::wrapPi(-3.0 * Constants::PI)), Constants::PI, 1e-15);
	EXPECT_NEAR(Helpers::wrapPi(100.0 * Constants::PI), 0.0, 1e-12);
}

TEST(Helpers, Matrix3x1ToVec3DRoundtrip)
{
	Vec3D v(1.5, -2.3, 4.7);
	auto mat = Helpers::vec3DToMatrix3x1(v);
	Vec3D back = Helpers::matrix3x1ToVec3D(mat);
	EXPECT_NEAR(back.x, v.x, 1e-15);
	EXPECT_NEAR(back.y, v.y, 1e-15);
	EXPECT_NEAR(back.z, v.z, 1e-15);
}

TEST(Helpers, EclEcvRoundtrip)
{
	Vec3D v(0.5, -0.3, 0.8);
	Vec3D ecl = Helpers::ecvToEcl(v);
	Vec3D back = Helpers::eclToEcv(ecl);
	EXPECT_NEAR(back.x, v.x, 1e-12);
	EXPECT_NEAR(back.y, v.y, 1e-12);
	EXPECT_NEAR(back.z, v.z, 1e-12);
}

TEST(Helpers, EcvToEclPreservesXAxis)
{
	Vec3D xAxis(1.0, 0.0, 0.0);
	Vec3D result = Helpers::ecvToEcl(xAxis);
	EXPECT_NEAR(result.x, 1.0, 1e-15);
	EXPECT_NEAR(result.y, 0.0, 1e-15);
	EXPECT_NEAR(result.z, 0.0, 1e-15);
}

TEST(Helpers, EcvToEclPreservesMagnitude)
{
	Vec3D v(3.0, 4.0, 5.0);
	Vec3D ecl = Helpers::ecvToEcl(v);
	EXPECT_NEAR(ecl.magnitude(), v.magnitude(), 1e-12);
}

TEST(Helpers, ApplyPrecessionMatrixAtJ2000)
{
	Vec3D v(1.0, 2.0, 3.0);
	Vec3D copy = v;
	Helpers::applyPrecessionMatrix(Constants::J2000_JD, copy);
	EXPECT_NEAR(copy.x, v.x, 1e-10);
	EXPECT_NEAR(copy.y, v.y, 1e-10);
	EXPECT_NEAR(copy.z, v.z, 1e-10);
}

TEST(Helpers, ApplyPrecessionMatrixPreservesMagnitude)
{
	Vec3D v(1.0, 0.0, 0.0);
	Helpers::applyPrecessionMatrix(2458634.5, v);
	EXPECT_NEAR(v.magnitude(), 1.0, 1e-4);
}

TEST(Calculations, KeplerSolveECircular)
{
	for (double M : {0.0, 0.5, 1.0, Constants::PI / 2.0, Constants::PI, 3.0, 2.0 * Constants::PI - 0.1}) {
		EXPECT_NEAR(Calculations::kepler_solve_e(M, 0.0), M, 1e-10);
	}
}

TEST(Calculations, KeplerSolveETrivial)
{
	EXPECT_NEAR(Calculations::kepler_solve_e(0.0, 0.0), 0.0, 1e-10);
	EXPECT_NEAR(Calculations::kepler_solve_e(0.0, 0.5), 0.0, 1e-10);
	EXPECT_NEAR(Calculations::kepler_solve_e(0.0, 0.9), 0.0, 1e-10);
}

TEST(Calculations, KeplerSolveEConsistency)
{
	double e = 0.5;
	for (double M : {0.1, 0.5, 1.0, 2.0, Constants::PI, 5.0}) {
		double E = Calculations::kepler_solve_e(M, e);
		double M_check = E - e * std::sin(E);
		EXPECT_NEAR(M_check, M, 1e-10) << "Failed for M=" << M << " e=" << e;
	}
}

TEST(Calculations, KeplerSolveEHighEccentricity)
{
	double e = 0.99;
	for (double M : {0.01, 0.1, 1.0, Constants::PI}) {
		double E = Calculations::kepler_solve_e(M, e);
		double M_check = E - e * std::sin(E);
		EXPECT_NEAR(M_check, M, 1e-9) << "Failed for M=" << M << " e=" << e;
	}
}

TEST(Calculations, KeplerSolveHConsistency)
{
	double e = 2.0;
	for (double M : {0.1, 0.5, 1.0, 3.0, 10.0}) {
		double F = Calculations::kepler_solve_h(M, e);
		double M_check = e * std::sinh(F) - F;
		EXPECT_NEAR(M_check, M, 1e-8) << "Failed for M=" << M << " e=" << e;
	}
}

TEST(Calculations, AngMomentum)
{
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, 1.0, 0.0);
	Vec3D H = Calculations::angMomentum(R, V);
	EXPECT_NEAR(H.x, 0.0, 1e-15);
	EXPECT_NEAR(H.y, 0.0, 1e-15);
	EXPECT_NEAR(H.z, 1.0, 1e-15);
}

TEST(Calculations, SpecificEnergy)
{
	double mu = 1.0;
	double r = 1.0;
	double v = std::sqrt(mu / r);
	Vec3D R(r, 0.0, 0.0);
	Vec3D V(0.0, v, 0.0);
	double E = Calculations::specificEnergy(R, V, mu);
	EXPECT_NEAR(E, -mu / (2.0 * r), 1e-12);
}

TEST(Calculations, SemiMajorAxis)
{
	double mu = 1.0;
	double a_expected = 2.0;
	double E = -mu / (2.0 * a_expected);
	EXPECT_NEAR(Calculations::semiMajorAxis(mu, E), a_expected, 1e-12);
}

TEST(Calculations, InclinationPrograde)
{
	Vec3D H(0.0, 0.0, 1.0);
	EXPECT_NEAR(Calculations::inclination(H), 0.0, 1e-12);
}

TEST(Calculations, InclinationPolar)
{
	Vec3D H(1.0, 0.0, 0.0);
	EXPECT_NEAR(Calculations::inclination(H), Constants::PI / 2.0, 1e-12);
}

TEST(Calculations, InclinationRetrograde)
{
	Vec3D H(0.0, 0.0, -1.0);
	EXPECT_NEAR(Calculations::inclination(H), Constants::PI, 1e-12);
}

TEST(Calculations, CircularOrbitElements)
{
	Vec3D R(1.0, 0.0, 0.0);
	Vec3D V(0.0, 1.0, 0.0);

	double E = Calculations::specificEnergy(R, V, Constants::MU_IN_UNIT);
	double a = Calculations::semiMajorAxis(Constants::MU_IN_UNIT, E);
	EXPECT_NEAR(a, 1.0, 1e-12);

	Vec3D eVec = Calculations::eccentricityVec(R, V);
	EXPECT_NEAR(eVec.magnitude(), 0.0, 1e-12);

	Vec3D H = Calculations::angMomentum(R, V);
	EXPECT_NEAR(Calculations::inclination(H), 0.0, 1e-12);
}

TEST(Calculations, EllipticalOrbitConsistency)
{
	double a = 2.0, e_val = 0.5;
	double r_peri = a * (1.0 - e_val);
	double v_peri = std::sqrt(Constants::MU_IN_UNIT * (2.0 / r_peri - 1.0 / a));

	Vec3D R(r_peri, 0.0, 0.0);
	Vec3D V(0.0, v_peri, 0.0);

	Vec3D eVec = Calculations::eccentricityVec(R, V);
	EXPECT_NEAR(eVec.magnitude(), e_val, 1e-10);

	double E = Calculations::specificEnergy(R, V, Constants::MU_IN_UNIT);
	double a_calc = Calculations::semiMajorAxis(Constants::MU_IN_UNIT, E);
	EXPECT_NEAR(a_calc, a, 1e-10);

	EXPECT_NEAR(Calculations::perigee(eVec, a_calc), r_peri, 1e-10);
}

TEST(Calculations, MeanAnomalyAtPerihelion)
{
	Vec3D eVec(0.5, 0.0, 0.0);
	double E_anom = Calculations::excentricAnomaly(eVec, 0.0, OrbitType::Elliptical);
	EXPECT_NEAR(E_anom, 0.0, 1e-10);

	double M = Calculations::meanAnomaly(eVec, E_anom, OrbitType::Elliptical);
	EXPECT_NEAR(M, 0.0, 1e-10);
}

TEST(Calculations, KeplerMeanAnomalyRoundtrip)
{
	double e_val = 0.3;
	double M_orig = 1.5;
	Vec3D eVec(e_val, 0.0, 0.0);

	double E = Calculations::kepler_solve_e(M_orig, e_val);
	double M_back = E - e_val * std::sin(E);
	EXPECT_NEAR(M_back, M_orig, 1e-10);
}

TEST(Calculations, Perigee)
{
	Vec3D eVec(0.5, 0.0, 0.0);
	double a = 2.0;
	EXPECT_NEAR(Calculations::perigee(eVec, a), 1.0, 1e-12);
}

TEST(Calculations, LongitudeOfAscendingNode)
{
	Vec3D H(0.0, -1.0, 1.0);
	double Omega = Calculations::longitudeOfAscendingNode(H);
	EXPECT_NEAR(std::fmod(Omega, 2.0 * Constants::PI), 0.0, 1e-10);
}
