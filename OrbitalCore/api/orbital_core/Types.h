#pragma once

#include <algorithm>
#include <cmath>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include "orbital_core/Vec3D.h"
#include "orbital_core/Matrix2D.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/Constants.h"

ORBITAL_CORE_API std::filesystem::path& getExecutableDir();

class IEphemerisSource;
ORBITAL_CORE_API IEphemerisSource* getDefaultEphemerisSource();
ORBITAL_CORE_API void setDefaultEphemerisSource(IEphemerisSource* source);

enum class OrbitType : std::uint8_t { Undef, Elliptical, Hyperbolic, Parabolic };

namespace Helpers {
	inline int mpcBase62Value(char ch) noexcept
	{
		if (ch >= '0' && ch <= '9')
			return ch - '0';
		if (ch >= 'A' && ch <= 'Z')
			return 10 + ch - 'A';
		if (ch >= 'a' && ch <= 'z')
			return 36 + ch - 'a';
		return -1;
	}

	inline std::string unpackMpcPackedDesignation(std::string_view name)
	{
		if (name.size() == 7 && (name[0] == 'I' || name[0] == 'J' || name[0] == 'K') &&
		    std::isdigit(static_cast<unsigned char>(name[1])) && std::isdigit(static_cast<unsigned char>(name[2])) &&
		    std::isalpha(static_cast<unsigned char>(name[3])) && std::isdigit(static_cast<unsigned char>(name[5])) &&
		    std::isalpha(static_cast<unsigned char>(name[6]))) {
			const int century = (name[0] == 'I') ? 1800 : (name[0] == 'J') ? 1900 : 2000;
			const int year = century + (name[1] - '0') * 10 + (name[2] - '0');
			const int cycle = mpcBase62Value(name[4]) * 10 + (name[5] - '0');
			if (cycle >= 0) {
				std::ostringstream os;
				os << year << name[3] << name[6];
				if (cycle > 0)
					os << cycle;
				return os.str();
			}
		}

		return std::string(name);
	}

	inline Vec3D matrix3x1ToVec3D(const Matrix2D<double, 3, 1>& mat) noexcept
	{
		return Vec3D(mat(0, 0), mat(1, 0), mat(2, 0));
	}

	inline Matrix2D<double, 3, 1> vec3DToMatrix3x1(const Vec3D& vec) noexcept
	{
		return Matrix2D<double, 3, 1>({{vec.x}, {vec.y}, {vec.z}});
	}

	inline Vec3D eclToEcv(const Vec3D& ecl) noexcept
	{
		const Matrix2D<double, 3, 3> rotateMatrix(
		    {{1, 0, 0},
		     {0, std::cos(Constants::EARTH_INCL), std::sin(Constants::EARTH_INCL)},
		     {0, -std::sin(Constants::EARTH_INCL), std::cos(Constants::EARTH_INCL)}});

		return matrix3x1ToVec3D(rotateMatrix.inverse() * vec3DToMatrix3x1(ecl));
	}

	inline Vec3D ecvToEcl(const Vec3D& ecv) noexcept
	{
		const Matrix2D<double, 3, 3> rotateMatrix(
		    {{1, 0, 0},
		     {0, std::cos(Constants::EARTH_INCL), std::sin(Constants::EARTH_INCL)},
		     {0, -std::sin(Constants::EARTH_INCL), std::cos(Constants::EARTH_INCL)}});

		return matrix3x1ToVec3D(rotateMatrix * vec3DToMatrix3x1(ecv));
	}

	inline void applyPrecessionMatrix(const double& JD, Vec3D& obsVector)
	{
		const double T = (JD - Constants::J2000_JD) / Constants::JULIAN_CENTURY_DAYS;

		// Precession angles (arcseconds → radians), Horner-form polynomials.
		constexpr double asec2rad = Constants::ARCSEC_RAD;
		const double zeta =
		    asec2rad *
		    (2.650545 + T * (2306.083227 + T * (0.2988499 + T * (0.01801828 + T * (-5.971e-6 + T * (-3.173e-7))))));
		const double z =
		    asec2rad *
		    (-2.650545 + T * (2306.077181 + T * (1.0927348 + T * (0.01826837 + T * (-28.596e-6 + T * (-2.904e-7))))));
		const double th =
		    asec2rad * T * (2004.191903 + T * (-0.4294934 + T * (-0.04182264 + T * (-7.089e-6 + T * (-1.274e-7)))));

		// Compose P = R3(−z_A) · R2(θ_A) · R3(−ζ_A).
		const double cZ = std::cos(zeta), sZ = std::sin(zeta);
		const double cz = std::cos(z), sz = std::sin(z);
		const double ct = std::cos(th), st = std::sin(th);

		// Matrix elements of P (J2000 → date).
		const double P11 = cz * ct * cZ - sz * sZ;
		const double P12 = -cz * ct * sZ - sz * cZ;
		const double P13 = -cz * st;
		const double P21 = sz * ct * cZ + cz * sZ;
		const double P22 = -sz * ct * sZ + cz * cZ;
		const double P23 = -sz * st;
		const double P31 = st * cZ;
		const double P32 = -st * sZ;
		const double P33 = ct;

		// Apply P^T (date → J2000): v_new[i] = Σ_j P_ji · v_old[j].
		const double nx = P11 * obsVector.x + P21 * obsVector.y + P31 * obsVector.z;
		const double ny = P12 * obsVector.x + P22 * obsVector.y + P32 * obsVector.z;
		const double nz = P13 * obsVector.x + P23 * obsVector.y + P33 * obsVector.z;
		obsVector.x = nx;
		obsVector.y = ny;
		obsVector.z = nz;
	}

	inline double wrapPi(double x) noexcept
	{
		x = std::remainder(x, Constants::TWOPI);
		return x;
	}

	inline double wrapTwoPi(double x) noexcept
	{
		x = std::fmod(x, Constants::TWOPI);
		return x < 0.0 ? x + Constants::TWOPI : x;
	}
};  // namespace Helpers

namespace Convertors {
	inline double degreesToRadians(const double& degrees) { return degrees * Constants::DEG_TO_RAD; }

	// Seconds to add to UTC to obtain TT/TDT.  Reference treats TDT, TT and ET
	// as equivalent in this conversion path; no periodic TDB-TT term is added.
	// Before leap seconds, UTC has a frequency drift. Pre-1960 UTC is undefined.
	// Historical observation UT1 is handled separately through TimeScale.
	inline double utcToTdtSeconds(double jd_utc) noexcept
	{
		if (!std::isfinite(jd_utc) || jd_utc < 2436934.5)
			return std::numeric_limits<double>::quiet_NaN();
		if (jd_utc < 2441317.5) {
			static constexpr struct {
				double mjd0, offset, reference, drift;
			} utcDrift[] = {
			    {39887, 4.2131700, 39126, 0.0025920}, {39126, 4.3131700, 39126, 0.0025920},
			    {39004, 3.8401300, 38761, 0.0012960}, {38942, 3.7401300, 38761, 0.0012960},
			    {38820, 3.6401300, 38761, 0.0012960}, {38761, 3.5401300, 38761, 0.0012960},
			    {38639, 3.4401300, 38761, 0.0012960}, {38486, 3.3401300, 38761, 0.0012960},
			    {38395, 3.2401300, 38761, 0.0012960}, {38334, 1.9458580, 37665, 0.0011232},
			    {37665, 1.8458580, 37665, 0.0011232}, {37512, 1.3728180, 37300, 0.0012960},
			    {37300, 1.4228180, 37300, 0.0012960}, {36934, 1.4178180, 37300, 0.0012960},
			};

			const double mjd = jd_utc - 2400000.5;
			for (const auto& entry : utcDrift)
				if (mjd >= entry.mjd0)
					return 32.184 + entry.offset + (mjd - entry.reference) * entry.drift;
		}

		static constexpr struct {
			double jd0;
			int dat;
		} leapSec[] = {
		    {2457754.5, 37}, {2457204.5, 36}, {2456109.5, 35}, {2454832.5, 34}, {2453736.5, 33}, {2451179.5, 32},
		    {2450630.5, 31}, {2450083.5, 30}, {2449534.5, 29}, {2449169.5, 28}, {2448804.5, 27}, {2448257.5, 26},
		    {2447892.5, 25}, {2447161.5, 24}, {2446247.5, 23}, {2445516.5, 22}, {2445151.5, 21}, {2444786.5, 20},
		    {2444239.5, 19}, {2443874.5, 18}, {2443509.5, 17}, {2443144.5, 16}, {2442778.5, 15}, {2442413.5, 14},
		    {2442048.5, 13}, {2441683.5, 12}, {2441499.5, 11}, {2441317.5, 10},
		};

		int dat = 10;
		for (const auto& ls : leapSec) {
			if (jd_utc >= ls.jd0) {
				dat = ls.dat;
				break;
			}
		}
		return static_cast<double>(dat) + 32.184;
	}

	inline double radiansToDegrees(const double& radians) { return radians * 180.0 / Constants::PI; }

	inline double convertKmToAU(const double& km) { return km / 149597870.7; }

	inline double convertAUToKm(const double& au) { return au * 149597870.7; }

	inline double convertGregorianToJD(int year, int month, double day) noexcept
	{
		if (month <= 2) {
			year -= 1;
			month += 12;
		}

		const int A = static_cast<int>(year / 100);
		const int B = 2 - A + static_cast<int>(A / 4);

		const double JD = std::floor(365.25 * (year + 4716)) + std::floor(30.6001 * (month + 1)) + day + B - 1524.5;

		return JD;
	}

	inline void radiansToHMS(const double& rad, int& hours, int& minutes, double& seconds)
	{
		double totalHours = rad * (12.0 / Constants::PI);

		int sign = (totalHours < 0) ? -1 : 1;

		totalHours = std::fabs(totalHours);
		hours = static_cast<int>(totalHours);

		double fractionalHours = totalHours - hours;
		minutes = static_cast<int>(fractionalHours * 60.0);

		double fractionalMinutes = fractionalHours * 60.0 - minutes;
		seconds = fractionalMinutes * 60.0;

		hours *= sign;
	}

	inline void radiansToDMS(const double& rad, int& degrees, int& minutes, double& seconds)
	{
		double totalDegrees = radiansToDegrees(rad);
		int sign = (totalDegrees < 0) ? -1 : 1;

		double absDegrees = std::fabs(totalDegrees);
		degrees = static_cast<int>(absDegrees);

		double fractionalDegrees = absDegrees - degrees;

		double totalMinutes = fractionalDegrees * 60.0;
		minutes = static_cast<int>(totalMinutes);
		seconds = (totalMinutes - minutes) * 60.0;

		degrees *= sign;
	}

	inline void jdToYMD(double jd, int& year, int& month, double& dayFrac)
	{
		const double jdh = jd + 0.5;
		const int Z = static_cast<int>(std::floor(jdh));
		const double F = jdh - Z;

		int A;
		if (Z >= 2299161) {
			const int alpha = static_cast<int>((Z - 1867216.25) / 36524.25);
			A = Z + 1 + alpha - alpha / 4;
		} else {
			A = Z;
		}

		const int B = A + 1524;
		const int C = static_cast<int>((B - 122.1) / 365.25);
		const int D = static_cast<int>(365.25 * C);
		const int E = static_cast<int>((B - D) / 30.6001);

		dayFrac = (B - D - static_cast<int>(30.6001 * E)) + F;
		month = (E < 14) ? E - 1 : E - 13;
		year = (month > 2) ? C - 4716 : C - 4715;
	}

	inline std::string convertJDToGregorian(const double& jd)
	{
		int year, month;
		double dayFrac;
		jdToYMD(jd, year, month, dayFrac);

		constexpr int centisecondsPerDay = 8640000;
		int day = static_cast<int>(std::floor(dayFrac));
		int centiseconds = static_cast<int>(std::llround((dayFrac - day) * centisecondsPerDay));
		if (centiseconds == centisecondsPerDay) {
			jdToYMD(std::floor(jd + 0.5) + 0.5, year, month, dayFrac);
			day = static_cast<int>(std::floor(dayFrac));
			centiseconds = 0;
		}
		const int hour = centiseconds / 360000;
		const int minute = centiseconds / 6000 % 60;
		const int second = centiseconds / 100 % 60;
		const int fraction = centiseconds % 100;

		std::ostringstream oss;
		oss << year << "-" << std::setfill('0') << std::setw(2) << month << "-" << std::setfill('0') << std::setw(2)
		    << day << " " << std::setfill('0') << std::setw(2) << hour << ":" << std::setfill('0') << std::setw(2)
		    << minute << ":" << std::setw(2) << second << "." << std::setw(2) << fraction;

		return oss.str();
	}

	inline std::string jdToPackedDate(double jd)
	{
		int year, month;
		double dayFrac;
		jdToYMD(jd, year, month, dayFrac);

		std::ostringstream oss;
		oss << std::setfill('0') << std::setw(2) << (year % 100) << std::setw(2) << month << " " << std::fixed
		    << std::setprecision(5) << std::setw(8) << dayFrac;
		return oss.str();
	}

	// The one string-to-number conversion in the project.
	//
	// std::stod / std::stoi answer "is this a number?" by throwing, which in a
	// parser that probes fields is control flow: every non-numeric token costs
	// a throw and a catch, and the debugger reports each one.  That noise
	// trains you to scroll past "Exception thrown", which is how a real
	// exception gets missed.  std::from_chars reports failure in its return
	// value instead.
	//
	// Accepts exactly what std::stod and std::stoi accepted, because the MPC
	// readers depend on it: leading whitespace and a leading '+' are skipped,
	// parsing stops at the first character that cannot continue the number,
	// and success means at least one character was consumed.  `consumed`
	// reports how far it got, counted from the start of `text`, so a caller
	// that requires the whole token can compare it against text.size().
	template <class T>
	[[nodiscard]] inline bool tryParseNumber(std::string_view text, T& out, std::size_t* consumed = nullptr) noexcept
	{
		static_assert(std::is_arithmetic_v<T>, "tryParseNumber handles numbers");

		std::size_t begin = 0;
		while (begin < text.size() && (text[begin] == ' ' || text[begin] == '	')) {
			++begin;
		}

		// from_chars rejects a leading '+'; stod and stoi accepted one.
		if (begin < text.size() && text[begin] == '+')
			++begin;

		if (begin >= text.size())
			return false;

		const char* const first = text.data() + begin;
		const char* const last = text.data() + text.size();

		T value{};
		const auto result = std::from_chars(first, last, value);
		if (result.ec != std::errc{})
			return false;

		out = value;
		if (consumed != nullptr)
			*consumed = static_cast<std::size_t>(result.ptr - text.data());
		return true;
	}

	// Value or zero, for the fixed-width readers where a blank field legitimately
	// means zero.  Callers that must tell "absent" from "zero" use
	// tryParseNumber directly.
	inline int stringToInt(const std::string& s) noexcept
	{
		int value = 0;
		return tryParseNumber(s, value) ? value : 0;
	}

	inline double stringToDouble(const std::string& s) noexcept
	{
		double value = 0.0;
		return tryParseNumber(s, value) ? value : 0.0;
	}

	inline double raToDegrees(const int& hours, const int& minutes, const double& seconds) noexcept
	{
		return hours * 15.0 + ((minutes * 15.0) / 60.0) + ((seconds * 15.0) / 3600.0);
	}

	inline double raToRadians(const int& hours, const int& minutes, const double& seconds) noexcept
	{
		return degreesToRadians(raToDegrees(hours, minutes, seconds));
	}

	inline double decToDegrees(int degrees, const int& minutes, const double& seconds)
	{
		// Проверяем знак Dec (+ или -)
		const double sign = (degrees < 0) ? -1.0 : 1.0;
		degrees = std::abs(degrees);

		double totalDegrees = degrees + minutes / 60.0 + seconds / 3600.0;
		totalDegrees *= sign;

		return totalDegrees;
	}

	inline double decToRadians(const int& degrees, const int& minutes, const double& seconds)
	{
		return degreesToRadians(decToDegrees(degrees, minutes, seconds));
	}

	inline Vec3D ecvToEcl(const Vec3D& vector, const double& incl) noexcept
	{
		return Vec3D(vector.x, vector.y * std::cos(incl) + vector.z * std::sin(incl),
		             vector.z * std::cos(incl) - vector.y * std::sin(incl));
	}
};  // namespace Convertors

namespace Calculations {
	inline double kepler_solve_e(const double& M, const double& e, double tolerance = 1e-10)
	{
		double E = M;
		double E_prev;

		auto kepler_eq = [&](const double& E, const double& M, const double& e) -> double {
			return E - e * std::sin(E) - M;
		};

		do {
			E_prev = E;
			E = E - kepler_eq(E, M, e) / (1 - e * std::cos(E));
		} while (std::fabs(E - E_prev) > tolerance);

		return E;
	}

	inline double kepler_solve_h(const double& M, const double& e, double tolerance = 1e-10)
	{
		if (!std::isfinite(M) || !std::isfinite(e) || e < 1.0 || !std::isfinite(tolerance) || tolerance <= 0.0)
			return std::numeric_limits<double>::quiet_NaN();
		if (M == 0.0)
			return 0.0;

		const double mean = std::abs(M);
		double lower = 0.0;
		double upper = std::max(1.0, std::asinh(mean / e) + 1.0);
		while (e * std::sinh(upper) - upper < mean)
			upper *= 2.0;
		double anomaly = std::asinh(mean / e);
		for (int i = 0; i < 100; ++i) {
			const double residual = e * std::sinh(anomaly) - anomaly - mean;
			if (residual == 0.0)
				break;
			if (residual > 0.0)
				upper = anomaly;
			else
				lower = anomaly;
			const double next = anomaly - residual / (e * std::cosh(anomaly) - 1.0);
			const double candidate = std::isfinite(next) && next > lower && next < upper ? next : (lower + upper) / 2.0;
			if (std::abs(candidate - anomaly) <= tolerance * std::max(1.0, std::abs(candidate))) {
				anomaly = candidate;
				break;
			}
			anomaly = candidate;
		}
		return std::copysign(anomaly, M);
	}

	inline Vec3D angMomentum(const Vec3D& R, const Vec3D& V) noexcept { return R.cross(V); }

	inline Vec3D eccentricityVec(const Vec3D& R, const Vec3D& V) noexcept
	{
		const Vec3D H = angMomentum(R, V);
		return V.cross(H) - (R / R.magnitude());
	}

	inline double specificEnergy(const Vec3D& R, const Vec3D& V, const double& mu) noexcept
	{
		return (V.magnitude() * V.magnitude()) / 2. - mu / R.magnitude();
	}

	inline double semiMajorAxis(const double& mu, const double& E) noexcept { return -mu / (2 * E); }

	inline double semiMinorAxis(const double& a, const double& e) noexcept
	{
		return std::abs(a) * std::sqrt(std::abs(1.0 - e * e));
	}

	inline double inclination(const Vec3D& H) noexcept
	{
		const double hm = H.magnitude();
		if (hm <= 0.0)
			return 0.0;

		const double c = std::clamp(H.z / hm, -1.0, 1.0);
		return std::acos(c);
	}

	inline double longitudeOfAscendingNode(const Vec3D& H) noexcept
	{
		const Vec3D k = Vec3D(0., 0., 1.);
		const Vec3D node = k.cross(H);
		const double& l = std::acos(node.x / node.magnitude());

		if (node.y > 0)
			return l;

		return Convertors::degreesToRadians(360.) - l;
	}

	inline double argumentOfPericentre(const Vec3D& E, const Vec3D& H) noexcept
	{
		const Vec3D k = Vec3D(0., 0., 1.);
		const Vec3D node = k.cross(H);

		const double& w = std::acos(node.dot(E) / (node.magnitude() * E.magnitude()));

		if (E.z > 0)
			return w;

		return Convertors::degreesToRadians(360.) - w;
	}

	inline double trueAnomaly(const Vec3D& E, const Vec3D& R, const Vec3D& V, const OrbitType orbitType) noexcept
	{
		const double& scalar = R.dot(V);

		if (orbitType == OrbitType::Elliptical) {

			if (scalar > 0)
				return std::acos((E.dot(R)) / (E.magnitude() * R.magnitude()));

			return Convertors::degreesToRadians(360.) - std::acos((E.dot(R)) / (E.magnitude() * R.magnitude()));
		}

		if (scalar > 0)
			return std::acos((E.dot(R)) / (E.magnitude() * R.magnitude()));

		return -std::acos((E.dot(R)) / (E.magnitude() * R.magnitude()));
	}

	inline double excentricAnomaly(const Vec3D& E, const double& realAnomaly, const OrbitType orbitType) noexcept
	{
		if (orbitType == OrbitType::Elliptical) {
			double anom =
			    2 * std::atan(std::sqrt((1 - E.magnitude()) / (1 + E.magnitude())) * std::tan(realAnomaly / 2));

			if (anom < 0)
				return anom + 2 * Constants::PI;

			return anom;
		}

		return 2 * std::atanh(std::sqrt((E.magnitude() - 1) / (E.magnitude() + 1)) * std::tan(realAnomaly / 2));
	}

	inline double meanAnomaly(const Vec3D& E, const double& e, const OrbitType orbitType) noexcept
	{
		if (orbitType == OrbitType::Elliptical)
			return e - E.magnitude() * std::sin(e);

		return E.magnitude() * std::sinh(e) - e;
	}

	// Time of perihelion passage Tp.
	//
	// Units: pass `mu` in AU³/day² (e.g. `Constants::GM_SUN`); `t0` is JD,
	// `meanAnom` is rad, `a` is AU.  Returns Tp as JD in the same time
	// scale as `t0`.
	//
	// Formula:  factor = |a|^(3/2) / sqrt(µ)   ⇒ days/rad (mean motion⁻¹).
	//   Elliptical, r·v>0:  body receding ⇒ past perihelion  ⇒  Tp = t0 − factor·M
	//   Elliptical, r·v<0:  body approaching ⇒ next perihelion ⇒ Tp = t0 + factor·(2π − M)
	//   Hyperbolic:         Tp = t0 − factor·M_h  (sign of M_h tells past/future)
	//
	// NB: previous version divided `factor` by `Constants::SECONDS_PER_DAY` —
	// a unit-conversion artefact that left over from a draft assuming
	// µ in CGS (cm³/s²).  With `µ = GM_SUN` (AU³/day²) no conversion is
	// needed; the `1/SECONDS_PER_DAY` divisor under-shrunk Tp by 5e6× and
	// pinned the displayed Tp to ≈t0 regardless of the orbit's geometry.
	inline double timeOfPerihelion(const double& t0, const Vec3D& R, const Vec3D& V, const double& meanAnom,
	                               const double& a, const double& mu, const OrbitType orbitType) noexcept
	{
		const double factor = (std::fabs(a) * std::sqrt(std::fabs(a))) / std::sqrt(mu);
		const double scalar = R.dot(V);

		if (orbitType == OrbitType::Elliptical) {
			if (scalar > 0)
				return t0 - factor * meanAnom;  // past perihelion
			else
				return t0 + factor * (Constants::TWOPI - meanAnom);  // next perihelion
		}

		// Hyperbolic — single passage, sign of M_h chooses past/future.
		return t0 - factor * meanAnom;
	}

	inline double perigee(const Vec3D& E, const double& a) noexcept { return a * (1 - E.magnitude()); }

	// Apogee Q = a(1+e).  Defined only for elliptical orbits — hyperbolic
	// bodies asymptotically escape and have no apogee.  Returns NaN for
	// a ≤ 0 so UI / consumers can detect "not applicable" without sign
	// gymnastics.
	inline double apogee(const Vec3D& E, const double& a) noexcept
	{
		if (a <= 0.0)
			return std::numeric_limits<double>::quiet_NaN();
		return a * (1.0 + E.magnitude());
	}

	// Mean motion n = sqrt(µ / |a|³).  Physically meaningful for both
	// regimes (rate of mean-anomaly progression — for hyperbolic that's
	// the hyperbolic mean anomaly M_h).  Uses |a| so result is always
	// positive.
	inline double meanMotion(const double& a, const double& mu) noexcept
	{
		const double absA = std::fabs(a);
		return std::sqrt(mu / (absA * absA * absA));
	}

	// Orbital period P = 2π/n.  Defined only for elliptical orbits —
	// hyperbolic bodies never return.  Returns NaN for a ≤ 0.
	inline double orbitalPeriod(const double& a, const double& mu) noexcept
	{
		if (a <= 0.0)
			return std::numeric_limits<double>::quiet_NaN();
		return Constants::TWOPI / meanMotion(a, mu);
	}

	// 1/a — heliocentric "energy" parameter.  Sign distinguishes bound
	// (>0, elliptical) from unbound (<0, hyperbolic); both are physically
	// meaningful.
	inline double reciprocalSemiMajorAxis(const double& a) noexcept { return 1.0 / a; }

	inline double focalParam(const Vec3D& R, const Vec3D& V, const double& mu) noexcept
	{
		return (std::pow(R.dot(V), 2)) / mu;
	}

	/**
	 * Heliocentric distance from geocentric distance and elongation.
	 * Law of cosines on Sun-Earth-Object triangle:
	 *   r² = d_E² + Δ² - 2·d_E·Δ·cos(elongation)
	 *
	 * @param delta          Geocentric distance to object (AU)
	 * @param elongRad       Elongation angle (radians, at Earth between object and Sun)
	 * @param earthSunDist   Earth-Sun distance (AU, default 1.0 — good to ~2%)
	 * @return Heliocentric distance of object (AU)
	 */
	inline double heliocentricDistance(double delta, double elongRad, double earthSunDist = 1.0) noexcept
	{
		const double r2 = earthSunDist * earthSunDist + delta * delta - 2.0 * earthSunDist * delta * std::cos(elongRad);
		return r2 > 0.0 ? std::sqrt(r2) : 0.0;
	}

	/**
	 * Heliocentric (or body-centric) position from Keplerian elements at a
	 * given mean anomaly.  Handles both elliptical (a > 0, e < 1) and
	 * hyperbolic (a < 0, e > 1) orbits.  Replaces the former copies of the
	 * same math in OrbitOptimizer and NBodyPropagator.
	 *
	 * Returns position in the same length units as `a`.  Orientation is the
	 * ecliptic if Ω, ω, i are the standard heliocentric-ecliptic Keplerian
	 * angles.
	 */
	inline Vec3D positionFromElements(double a, double e, double i, double Omega, double omega, double meanAnomaly,
	                                  OrbitType type) noexcept
	{
		double r;
		double trueAnomaly;

		if (type == OrbitType::Elliptical) {
			const double modM =
			    std::fmod(std::fmod(meanAnomaly, Constants::TWOPI) + Constants::TWOPI, Constants::TWOPI);
			const double E = kepler_solve_e(modM, e);

			r = a * (1.0 - e * std::cos(E));
			trueAnomaly =
			    2.0 * std::atan2(std::sqrt(1.0 + e) * std::sin(E * 0.5), std::sqrt(1.0 - e) * std::cos(E * 0.5));
		}

		else {
			const double F = kepler_solve_h(meanAnomaly, e);
			r = a * (1.0 - e * std::cosh(F));  // a<0, cosh>=1 → r > 0
			trueAnomaly =
			    2.0 * std::atan2(std::sqrt(e + 1.0) * std::sinh(F * 0.5), std::sqrt(e - 1.0) * std::cosh(F * 0.5));
		}

		const double u = trueAnomaly + omega;
		const double cu = std::cos(u), su = std::sin(u);
		const double cO = std::cos(Omega), sO = std::sin(Omega);
		const double ci = std::cos(i), si = std::sin(i);
		return Vec3D((cu * cO - su * sO * ci) * r, (cu * sO + su * cO * ci) * r, (su * si) * r);
	}

	/**
	 * Propagate Keplerian mean anomaly from epoch to jd and return the
	 * corresponding heliocentric (ecliptic) position.  Time units: (jd, epoch)
	 * and `mu` must be consistent — e.g. jd/epoch in days + mu in AU³/day²
	 * (Constants::GAUSS_K2) gives position in AU.
	 */
	inline Vec3D positionAt(double a, double e, double i, double Omega, double omega, double M0, double epoch,
	                        double jd, double mu, OrbitType type) noexcept
	{
		const double absA = std::abs(a);
		const double n = std::sqrt(mu / (absA * absA * absA));
		const double M = M0 + n * (jd - epoch);
		return positionFromElements(a, e, i, Omega, omega, M, type);
	}

	/**
	 * TDB - TT in seconds at the given MJD (TT): Earth's orbital term through
	 * its eccentric anomaly plus the Jupiter term, ~1.7 ms amplitude. Used by
	 * TimeScale conversions and the radar clock model. Dynamics runs on TDB.
	 */
	inline double tdbMinusTtSeconds(double mjd) noexcept
	{
		const double days = mjd - 51544.5;
		const auto angle = [days](double phase, double rate) {
			return std::remainder(phase + rate * days, 360.0) * Constants::DEG_TO_RAD;
		};
		const double mean = angle(357.51716, 0.985599987);
		constexpr double e = 0.01671022;
		const double eccentric = mean + e * std::sin(mean) + 0.5 * e * e * std::sin(2 * mean);
		return 0.001658 * std::sin(eccentric) +
		       0.00002073 * std::sin(angle(280.46435, 0.985609100) - angle(34.40438, 0.083086762));
	}

	/**
	 * H-G phase function (1 - G) Phi1 + G Phi2 (Bowell et al. 1989) for the
	 * phase angle in radians.  Magnitudes use -2.5 log10 of this value.
	 */
	inline double hgPhaseFunction(double phaseRad, double slopeG) noexcept
	{
		const double tanHalf = std::tan(0.5 * phaseRad);
		const double phi1 = std::exp(-3.33 * std::pow(tanHalf, 0.63));
		const double phi2 = std::exp(-1.87 * std::pow(tanHalf, 1.22));
		return (1.0 - slopeG) * phi1 + slopeG * phi2;
	}
}  // namespace Calculations
