#include "EarthRotation.h"
#include "earth_orientation/EarthOrientation.h"
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "erfa.h"

#include <cmath>
#include <utility>

namespace {

	constexpr double U2R = 4.84813681109535994e-6 * 1.0e-7;
	constexpr double ASEC_TO_RAD = 4.84813681109535994e-6;

	inline double juliancentury(double jd_tt) noexcept
	{
		return (jd_tt - Constants::J2000_JD) / Constants::JULIAN_CENTURY_DAYS;
	}

	inline double wrap(double x) noexcept { return std::remainder(x, Constants::TWOPI); }

	inline double meanAnomalyMoon(double t) noexcept
	{
		return wrap(ASEC_TO_RAD * std::fmod(485868.249036 + t * (1717915923.2178 +
		                                                         t * (31.8792 + t * (0.051635 + t * (-0.00024470)))),
		                                    1296000.0));
	}

	inline double meanAnomalySun(double t) noexcept
	{
		return wrap(ASEC_TO_RAD *
		            std::fmod(1287104.79305 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))),
		                      1296000.0));
	}

	inline double meanArgLatMoon(double t) noexcept
	{
		return wrap(ASEC_TO_RAD * std::fmod(335779.526232 + t * (1739527262.8478 +
		                                                         t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))),
		                                    1296000.0));
	}

	inline double meanElongMoon(double t) noexcept
	{
		return wrap(ASEC_TO_RAD * std::fmod(1072260.70369 + t * (1602961601.2090 +
		                                                         t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))),
		                                    1296000.0));
	}

	inline double meanLongAscNode(double t) noexcept
	{
		return wrap(ASEC_TO_RAD *
		            std::fmod(450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939)))),
		                      1296000.0));
	}

	struct LuniSolarTerm {
		int nl, nlp, nF, nD, nOm;
		long Ap, Apt, Aps;
		long Bp, Bpt, Bps;
	};

	constexpr LuniSolarTerm kLunisolar[] = {
	    {0, 0, 0, 0, 1, -172064161L, -174666L, 33386L, 92052331L, 9086L, 15377L},
	    {0, 0, 2, -2, 2, -13170906L, -1675L, -13696L, 5730336L, -3015L, -4587L},
	    {0, 0, 2, 0, 2, -2276413L, -234L, 2796L, 978459L, -485L, 1374L},
	    {0, 0, 0, 0, 2, 2074554L, 207L, -698L, -897492L, 470L, -291L},
	    {0, 1, 0, 0, 0, 1475877L, -3633L, 11817L, 73871L, -184L, -1924L},
	    {0, 1, 2, -2, 2, -516821L, 1226L, -524L, 224386L, -677L, -174L},
	    {1, 0, 0, 0, 0, 711159L, 73L, -872L, -6750L, 0L, 358L},
	    {0, 0, 2, 0, 1, -387298L, -367L, 380L, 200728L, 18L, 318L},
	    {1, 0, 2, 0, 2, -301461L, -36L, 816L, 129025L, -63L, 367L},
	    {0, -1, 2, -2, 2, 215829L, -494L, 111L, -95929L, 299L, 132L},
	    {0, 0, 2, -2, 1, 128227L, 137L, 181L, -68982L, -9L, 39L},
	    {-1, 0, 2, 0, 2, 123457L, 11L, 19L, -53311L, 32L, -4L},
	    {-1, 0, 0, 2, 0, 156994L, 10L, -168L, -1235L, 0L, 82L},
	    {1, 0, 0, 0, 1, 63110L, 63L, 27L, -33228L, 0L, -9L},
	    {-1, 0, 0, 0, 1, -57976L, -63L, -189L, 31429L, 0L, -75L},
	    {-1, 0, 2, 2, 2, -59641L, -11L, 149L, 25543L, -11L, 66L},
	    {1, 0, 2, 0, 1, -51613L, -42L, 129L, 26366L, 0L, 78L},
	    {-2, 0, 2, 0, 1, 45893L, 50L, 31L, -24236L, -10L, 20L},
	    {0, 0, 0, 2, 0, 63384L, 11L, -150L, -1220L, 0L, 29L},
	    {0, 0, 2, 2, 2, -38571L, -1L, 158L, 16452L, -11L, 68L},
	    {0, -2, 2, -2, 2, 32481L, 0L, 0L, -13870L, 0L, 0L},
	    {-2, 0, 0, 2, 0, -47722L, 0L, -18L, 477L, 0L, -25L},
	    {2, 0, 2, 0, 2, -31046L, -1L, 131L, 13238L, -11L, 59L},
	    {1, 0, 2, -2, 2, 28593L, 0L, -1L, -12338L, 10L, -3L},
	    {-1, 0, 2, 0, 1, 20441L, 21L, 10L, -10758L, 0L, -3L},
	    {2, 0, 0, 0, 0, 29243L, 0L, -74L, -609L, 0L, 13L},
	    {0, 0, 2, 0, 0, 25887L, 0L, -66L, -550L, 0L, 11L},
	    {0, 1, 0, 0, 1, -14053L, -25L, 79L, 8551L, -2L, -45L},
	    {-1, 0, 0, 2, 1, 15164L, 10L, 11L, -8001L, 0L, -1L},
	    {0, 2, 2, -2, 2, -15794L, 72L, -16L, 6850L, -42L, -5L},
	    {0, 0, -2, 2, 0, 21783L, 0L, 13L, -167L, 0L, 13L},
	    {1, 0, 0, -2, 1, -12873L, -10L, -37L, 6953L, 0L, -14L},
	    {0, -1, 0, 0, 1, -12654L, 11L, 63L, 6415L, 0L, 26L},
	    {-1, 0, 2, 2, 1, -10204L, 0L, 25L, 5222L, 0L, 15L},
	    {0, 2, 0, 0, 0, 16707L, -85L, -10L, 168L, -1L, 10L},
	    {1, 0, 2, 2, 2, -7691L, 0L, 44L, 3268L, 0L, 19L},
	    {-2, 0, 2, 0, 0, -11024L, 0L, -14L, 104L, 0L, 2L},
	    {0, 1, 2, 0, 2, 7566L, -21L, -11L, -3250L, 0L, -5L},
	    {0, 0, 2, 2, 1, -6637L, -11L, 25L, 3353L, 0L, 14L},
	    {0, -1, 2, 0, 2, -7141L, 21L, 8L, 3070L, 0L, 4L},
	    {0, 0, 0, 2, 1, -6302L, -11L, 2L, 3272L, 0L, 4L},
	    {1, 0, 2, -2, 1, 5800L, 10L, 2L, -3045L, 0L, -1L},
	    {2, 0, 2, -2, 2, 6443L, 0L, -7L, -2768L, 0L, -4L},
	    {-2, 0, 0, 2, 1, -5774L, -11L, -15L, 3041L, 0L, -5L},
	    {2, 0, 2, 0, 1, -5350L, 0L, 21L, 2695L, 0L, 12L},
	    {0, -1, 2, -2, 1, -4752L, -11L, -3L, 2719L, 0L, -3L},
	    {0, 0, 0, -2, 1, -4940L, -11L, -21L, 2720L, 0L, -9L},
	    {-1, -1, 0, 2, 0, 7350L, 0L, -8L, -51L, 0L, 4L},
	    {2, 0, 0, -2, 1, 4065L, 0L, 6L, -2206L, 0L, 1L},
	    {1, 0, 0, 2, 0, 6579L, 0L, -24L, -199L, 0L, 2L},
	    {0, 1, 2, -2, 1, 3579L, 0L, 5L, -1900L, 0L, 1L},
	    {1, -1, 0, 0, 0, 4725L, 0L, -6L, -41L, 0L, 3L},
	    {-2, 0, 2, 0, 2, -3075L, 0L, -2L, 1313L, 0L, -1L},
	    {3, 0, 2, 0, 2, -2904L, 0L, 15L, 1233L, 0L, 7L},
	    {0, -1, 0, 2, 0, 4348L, 0L, -10L, -81L, 0L, 2L},
	    {1, -1, 2, 0, 2, -2878L, 0L, 8L, 1232L, 0L, 4L},
	    {0, 0, 0, 1, 0, -4230L, 0L, 5L, -20L, 0L, -2L},
	    {-1, -1, 2, 2, 2, -2819L, 0L, 7L, 1207L, 0L, 3L},
	    {-1, 0, 2, 0, 0, -4056L, 0L, 5L, 40L, 0L, -2L},
	    {0, -1, 2, 2, 2, -2647L, 0L, 11L, 1129L, 0L, 5L},
	    {-2, 0, 0, 0, 1, -2294L, 0L, -10L, 1266L, 0L, -4L},
	    {1, 1, 2, 0, 2, 2481L, 0L, -7L, -1062L, 0L, -3L},
	    {2, 0, 0, 0, 1, 2179L, 0L, -2L, -1129L, 0L, -2L},
	    {-1, 1, 0, 1, 0, 3276L, 0L, 1L, -9L, 0L, 0L},
	    {1, 1, 0, 0, 0, -3389L, 0L, 5L, 35L, 0L, -2L},
	    {1, 0, 2, 0, 0, 3339L, 0L, -13L, -107L, 0L, 1L},
	    {-1, 0, 2, -2, 1, -1987L, 0L, -6L, 1073L, 0L, -2L},
	    {1, 0, 0, 0, 2, -1981L, 0L, 0L, 854L, 0L, 0L},
	    {-1, 0, 0, 1, 0, 4026L, 0L, -353L, -553L, 0L, -139L},
	    {0, 0, 2, 1, 2, 1660L, 0L, -5L, -710L, 0L, -2L},
	    {-1, 0, 2, 4, 2, -1521L, 0L, 9L, 647L, 0L, 4L},
	    {-1, 1, 0, 1, 1, 1314L, 0L, 0L, -700L, 0L, 0L},
	    {0, -2, 2, -2, 1, -1283L, 0L, 0L, 672L, 0L, 0L},
	    {1, 0, 2, 2, 1, -1331L, 0L, 8L, 663L, 0L, 4L},
	    {-2, 0, 2, 2, 2, 1383L, 0L, -2L, -594L, 0L, -2L},
	    {-1, 0, 0, 0, 2, 1405L, 0L, 4L, -610L, 0L, 2L},
	    {1, 1, 2, -2, 2, 1290L, 0L, 0L, -556L, 0L, 0L},
	};

	constexpr int kNTerms = sizeof(kLunisolar) / sizeof(kLunisolar[0]);

	constexpr double kDpsiPlanetary = -0.135e-3 * ASEC_TO_RAD;
	constexpr double kDepsPlanetary = 0.388e-3 * ASEC_TO_RAD;

	inline void mulMatVec(const double M[3][3], Vec3D& v) noexcept
	{
		const double x = M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z;
		const double y = M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z;
		const double z = M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z;
		v.x = x;
		v.y = y;
		v.z = z;
	}

	inline void buildNutationMatrix(double epsA, double dpsi, double deps, double N[3][3]) noexcept
	{
		const double epsTrue = epsA + deps;
		const double cA = std::cos(epsA), sA = std::sin(epsA);
		const double cT = std::cos(epsTrue), sT = std::sin(epsTrue);
		const double cP = std::cos(dpsi), sP = std::sin(dpsi);

		N[0][0] = cP;
		N[0][1] = -sP * cA;
		N[0][2] = -sP * sA;

		N[1][0] = cT * sP;
		N[1][1] = cT * cP * cA + sT * sA;
		N[1][2] = cT * cP * sA - sT * cA;

		N[2][0] = sT * sP;
		N[2][1] = sT * cP * cA - cT * sA;
		N[2][2] = sT * cP * sA + cT * cA;
	}

	inline void transpose3(double M[3][3]) noexcept
	{
		std::swap(M[0][1], M[1][0]);
		std::swap(M[0][2], M[2][0]);
		std::swap(M[1][2], M[2][1]);
	}

	inline void buildPrecessionMatrix(double jd_tt, double P[3][3]) noexcept
	{
		const double T = juliancentury(jd_tt);
		const double zeta =
		    ASEC_TO_RAD *
		    (2.650545 + T * (2306.083227 + T * (0.2988499 + T * (0.01801828 + T * (-5.971e-6 + T * (-3.173e-7))))));
		const double z =
		    ASEC_TO_RAD *
		    (-2.650545 + T * (2306.077181 + T * (1.0927348 + T * (0.01826837 + T * (-28.596e-6 + T * (-2.904e-7))))));
		const double th =
		    ASEC_TO_RAD * T * (2004.191903 + T * (-0.4294934 + T * (-0.04182264 + T * (-7.089e-6 + T * (-1.274e-7)))));

		const double cZ = std::cos(zeta), sZ = std::sin(zeta);
		const double cz = std::cos(z), sz = std::sin(z);
		const double ct = std::cos(th), st = std::sin(th);

		P[0][0] = cz * ct * cZ - sz * sZ;
		P[0][1] = -cz * ct * sZ - sz * cZ;
		P[0][2] = -cz * st;
		P[1][0] = sz * ct * cZ + cz * sZ;
		P[1][1] = -sz * ct * sZ + cz * cZ;
		P[1][2] = -sz * st;
		P[2][0] = st * cZ;
		P[2][1] = -st * sZ;
		P[2][2] = ct;
	}

	inline double meanObliquityIau1976(double jdTt) noexcept
	{
		const double t = juliancentury(jdTt);
		const double arcsec = ((0.00181 * t - 0.0006) * t - 46.815) * t + 84381.448;
		return arcsec * ASEC_TO_RAD;
	}

	void nutationAnglesIau1980(double jdTt, double& dpsi, double& deps) noexcept
	{
		constexpr double rs = 4.84813681109536e-6;
		constexpr double twoPi = 6.2831853071795865;
		constexpr float z1 = 1.0f;
		constexpr float z2 = 2.0f;

		const double t1 = juliancentury(jdTt);
		const float t = static_cast<float>(t1);
		const float t2 = t * t;
		const float t3 = t2 * t;

		const double dl = (485866.733 + 1717915922.633 * t1 + 31.310 * t2 + 0.064 * t3) * rs;
		const double dp = (1287099.804 + 129596581.224 * t1 - 0.577 * t2 - 0.012 * t3) * rs;
		const double df = (335778.877 + 1739527263.137 * t1 - 13.257 * t2 + 0.011 * t3) * rs;
		const double dd = (1072261.307 + 1602961601.328 * t1 - 6.891 * t2 + 0.019 * t3) * rs;
		const double dn = (450160.280 - 6962890.539 * t1 + 7.455 * t2 + 0.008 * t3) * rs;

		const float l = static_cast<float>(std::fmod(dl, twoPi));
		const float p = static_cast<float>(std::fmod(dp, twoPi));
		const float x = static_cast<float>(std::fmod(df, twoPi) * 2.0);
		const float d = static_cast<float>(std::fmod(dd, twoPi));
		const float n = static_cast<float>(std::fmod(dn, twoPi));
		const float cl = std::cos(l), sl = std::sin(l);
		const float cp = std::cos(p), sp = std::sin(p);
		const float cx = std::cos(x), sx = std::sin(x);
		const float cd = std::cos(d), sd = std::sin(d);
		const float cn = std::cos(n), sn = std::sin(n);
		const float cp2 = z2 * cp * cp - z1;
		const float sp2 = z2 * sp * cp;
		const float cd2 = z2 * cd * cd - z1;
		const float sd2 = z2 * sd * cd;
		const float cn2 = z2 * cn * cn - z1;
		const float sn2 = z2 * sn * cn;
		const float cl2 = z2 * cl * cl - z1;
		const float sl2 = z2 * sl * cl;
		const float ca = cx * cd2 + sx * sd2;
		const float sa = sx * cd2 - cx * sd2;
		const float cb = ca * cn - sa * sn;
		const float sb = sa * cn + ca * sn;
		const float cc = cb * cn - sb * sn;
		const float sc = sb * cn + cb * sn;
		const float cv = cx * cd2 - sx * sd2;
		const float sv = sx * cd2 + cx * sd2;
		const float ce = cv * cn - sv * sn;
		const float se = sv * cn + cv * sn;
		const float cf = ce * cn - se * sn;
		const float sf = se * cn + ce * sn;
		const float cg = cl * cd2 + sl * sd2;
		const float sg = sl * cd2 - cl * sd2;
		const float ch = cx * cn2 - sx * sn2;
		const float sh = sx * cn2 + cx * sn2;
		const float cj = ch * cl - sh * sl;
		const float sj = sh * cl + ch * sl;
		const float ck = cj * cl - sj * sl;
		const float sk = sj * cl + cj * sl;
		const float cm = cx * cl2 + sx * sl2;
		const float sm = sx * cl2 - cx * sl2;
		const float cq = cl * cd + sl * sd;
		const float sq = sl * cd - cl * sd;
		const float cr = z2 * cq * cq - z1;
		const float sr = z2 * sq * cq;
		const float cs = cx * cn - sx * sn;
		const float ss = sx * cn + cx * sn;
		const float ct = cs * cl - ss * sl;
		const float st = ss * cl + cs * sl;
		const float cu = cf * cl + sf * sl;
		const float su = sf * cl - cf * sl;
		const float cw = cp * cg - sp * sg;
		const float sw = sp * cg + cp * sg;

		dpsi = -(171996.0f + 174.2f * t) * sn + (2062.0f + 0.2f * t) * sn2 + 46.0f * (sm * cn + cm * sn) - 11.0f * sm -
		       3.0f * (sm * cn2 + cm * sn2) - 3.0f * (sq * cp - cq * sp) - 2.0f * (sb * cp2 - cb * sp2) +
		       (sn * cm - cn * sm) - (13187.0f + 1.6f * t) * sc + (1426.0f - 3.4f * t) * sp -
		       (517.0f - 1.2f * t) * (sc * cp + cc * sp) + (217.0f - 0.5f * t) * (sc * cp - cc * sp) +
		       (129.0f + 0.1f * t) * sb + 48.0f * sr - 22.0f * sa + (17.0f - 0.1f * t) * sp2 -
		       15.0f * (sp * cn + cp * sn) - (16.0f - 0.1f * t) * (sc * cp2 + cc * sp2) - 12.0f * (sn * cp - cn * sp);
		dpsi = dpsi - 6.0f * (sn * cr - cn * sr) - 5.0f * (sb * cp - cb * sp) + 4.0f * (sr * cn + cr * sn) +
		       4.0f * (sb * cp + cb * sp) - 4.0f * sq + (sr * cp + cr * sp) + (sn * ca - cn * sa) -
		       (sp * ca - cp * sa) + (sp * cn2 + cp * sn2) + (sn * cq - cn * sq) - (sp * ca + cp * sa) -
		       (2274.0f + 0.2f * t) * sh + (712.0f + 0.1f * t) * sl - (386.0f + 0.4f * t) * ss - 301.0f * sj -
		       158.0f * sg + 123.0f * (sh * cl - ch * sl) + 63.0f * sd2 + (63.0f + 0.1f * t) * (sl * cn + cl * sn) -
		       (58.0f + 0.1f * t) * (sn * cl - cn * sl) - 59.0f * su - 51.0f * st - 38.0f * sf + 29.0f * sl2;
		dpsi = dpsi + 29.0f * (sc * cl + cc * sl) - 31.0f * sk + 26.0f * sx + 21.0f * (ss * cl - cs * sl) +
		       16.0f * (sn * cg - cn * sg) - 13.0f * (sn * cg + cn * sg) - 10.0f * (se * cl - ce * sl) -
		       7.0f * (sg * cp + cg * sp) + 7.0f * (sh * cp + ch * sp) - 7.0f * (sh * cp - ch * sp) -
		       8.0f * (sf * cl + cf * sl) + 6.0f * (sl * cd2 + cl * sd2) + 6.0f * (sc * cl2 + cc * sl2) -
		       6.0f * (sn * cd2 + cn * sd2) - 7.0f * se + 6.0f * (sb * cl + cb * sl) - 5.0f * (sn * cd2 - cn * sd2) +
		       5.0f * (sl * cp - cl * sp) - 5.0f * (ss * cl2 + cs * sl2) - 4.0f * (sp * cd2 - cp * sd2);
		dpsi = dpsi + 4.0f * (sl * cx - cl * sx) - 4.0f * sd - 3.0f * (sl * cp + cl * sp) + 3.0f * (sl * cx + cl * sx) -
		       3.0f * (sj * cp - cj * sp) - 3.0f * (su * cp - cu * sp) - 2.0f * (sn * cl2 - cn * sl2) -
		       3.0f * (sk * cl + ck * sl) - 3.0f * (sf * cp - cf * sp) + 2.0f * (sj * cp + cj * sp) -
		       2.0f * (sb * cl - cb * sl);
		dpsi = dpsi + 2.0f * (sn * cl2 + cn * sl2) - 2.0f * (sl * cn2 + cl * sn2) + 2.0f * (sl * cl2 + cl * sl2) +
		       2.0f * (sh * cd + ch * sd) + (sn2 * cl - cn2 * sl) - (sg * cd2 - cg * sd2) + (sf * cl2 - cf * sl2) -
		       2.0f * (su * cd2 + cu * sd2) - (sr * cd2 - cr * sd2) + (sw * ch + cw * sh) - (sl * ce + cl * se) -
		       (sf * cr - cf * sr) + (su * ca + cu * sa) + (sg * cp - cg * sp) + (sb * cl2 + cb * sl2) -
		       (sf * cl2 + cf * sl2) - (st * ca - ct * sa) + (sc * cx + cc * sx) + (sj * cr + cj * sr) -
		       (sg * cx + cg * sx);
		dpsi = dpsi + (sp * cs + cp * ss) + (sn * cw - cn * sw) - (sn * cx - cn * sx) - (sh * cd - ch * sd) -
		       (sp * cd2 + cp * sd2) - (sl * cv - cl * sv) - (ss * cp - cs * sp) - (sw * cn + cw * sn) -
		       (sl * ca - cl * sa) + (sl2 * cd2 + cl2 * sd2) - (sf * cd2 + cf * sd2) + (sp * cd + cp * sd);

		deps = (92025.0f + 8.9f * t) * cn - (895.0f - 0.5f * t) * cn2 - 24.0f * (cm * cn - sm * sn) +
		       (cm * cn2 - sm * sn2) + (cb * cp2 + sb * sp2) + (5736.0f - 3.1f * t) * cc + (54.0f - 0.1f * t) * cp +
		       (224.0f - 0.6f * t) * (cc * cp - sc * sp) - (95.0f - 0.3f * t) * (cc * cp + sc * sp) - 70.0f * cb + cr +
		       9.0f * (cp * cn - sp * sn) + 7.0f * (cc * cp2 - sc * sp2) + 6.0f * (cn * cp + sn * sp) +
		       3.0f * (cn * cr + sn * sr) + 3.0f * (cb * cp + sb * sp) - 2.0f * (cr * cn - sr * sn) -
		       2.0f * (cb * cp - sb * sp);
		deps = deps + (977.0f - 0.5f * t) * ch - 7.0f * cl + 200.0f * cs + (129.0f - 0.1f * t) * cj - cg -
		       53.0f * (ch * cl + sh * sl) - 2.0f * cd2 - 33.0f * (cl * cn - sl * sn) + 32.0f * (cn * cl + sn * sl) +
		       26.0f * cu + 27.0f * ct + 16.0f * cf - cl2 - 12.0f * (cc * cl - sc * sl) + 13.0f * ck - cx -
		       10.0f * (cs * cl + ss * sl) - 8.0f * (cn * cg + sn * sg) + 7.0f * (cn * cg - sn * sg) +
		       5.0f * (ce * cl + se * sl) - 3.0f * (ch * cp - sh * sp) + 3.0f * (ch * cp + sh * sp) +
		       3.0f * (cf * cl - sf * sl) - 3.0f * (cc * cl2 - sc * sl2) + 3.0f * (cn * cd2 - sn * sd2) + 3.0f * ce;
		deps = deps - 3.0f * (cb * cl - sb * sl) + 3.0f * (cn * cd2 + sn * sd2) + 3.0f * (cs * cl2 - ss * sl2) +
		       (cj * cp + sj * sp) + (cu * cp + su * sp) + (cn * cl2 + sn * sl2) + (ck * cl - sk * sl) +
		       (cf * cp + sf * sp) - (cj * cp - sj * sp) + (cb * cl + sb * sl);
		deps = deps - (cn * cl2 - sn * sl2) + (cl * cn2 - sl * sn2) - (ch * cd - sh * sd) - (cn2 * cl + sn2 * sl) -
		       (cf * cl2 + sf * sl2) + (cu * cd2 - su * sd2) - (cw * ch - sw * sh) + (cl * ce - sl * se) +
		       (cf * cr + sf * sr) - (cb * cl2 - sb * sl2);

		dpsi *= 1.0e-4;
		deps *= 1.0e-4;
	}

	struct WahrTerm {
		int l, lp, f, d, om;
		double psi, psiT, eps, epsT;
	};

	constexpr WahrTerm kWahr[] = {
	    {3, 0, 0, 0, 0, 2.0f, 0.0f, 0.0f, 0.0f},
	    {2, 1, 0, -2, 0, 1.0f, 0.0f, 0.0f, 0.0f},
	    {2, 0, -2, 0, 0, 11.0f, 0.0f, 0.0f, 0.0f},
	    {2, 0, 0, -2, 0, 48.0f, 0.0f, 1.0f, 0.0f},
	    {2, 0, 0, -4, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {2, 0, 0, 2, 0, 1.0f, 0.0f, 0.0f, 0.0f},
	    {2, 0, 0, 0, 0, 29.0f, 0.0f, -1.0f, 0.0f},
	    {1, -1, 0, -1, 0, -3.0f, 0.0f, 0.0f, 0.0f},
	    {1, -1, 0, -2, 0, 1.0f, 0.0f, 0.0f, 0.0f},
	    {1, -1, 0, 0, 0, 5.0f, 0.0f, 0.0f, 0.0f},
	    {1, 1, 0, -2, 0, -7.0f, 0.0f, 0.0f, 0.0f},
	    {1, 1, 0, 0, 0, -3.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, -2, -2, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, -2, 2, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, -2, 0, 0, 4.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 2, -2, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 2, 0, 0, 3.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 0, -1, 0, -4.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 0, -2, 0, -158.0f, 0.0f, -1.0f, 0.0f},
	    {1, 0, 0, -4, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 0, 2, 0, 6.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 0, 0, 0, 712.0f, 0.1f, -7.0f, 0.0f},
	    {0, 1, -2, 2, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 1, 2, -2, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 1, 0, -2, 0, -4.0f, 0.0f, 0.0f, 0.0f},
	    {0, 1, 0, 2, 0, -1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 1, 0, 1, 0, 1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, 2, -2, 0, -22.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, 2, 0, 0, 26.0f, 0.0f, -1.0f, 0.0f},
	    {0, 0, 0, 2, 0, 63.0f, 0.0f, -2.0f, 0.0f},
	    {0, 0, 0, 1, 0, -4.0f, 0.0f, 0.0f, 0.0f},
	    {-1, -1, 0, 2, 1, 1.0f, 0.0f, 0.0f, 0.0f},
	    {-1, 0, 2, -2, 1, -2.0f, 0.0f, 1.0f, 0.0f},
	    {-1, 0, 2, 2, 1, -10.0f, 0.0f, 5.0f, 0.0f},
	    {-1, 0, 2, 0, 1, 21.0f, 0.0f, -10.0f, 0.0f},
	    {0, -2, 2, -2, 1, -2.0f, 0.0f, 1.0f, 0.0f},
	    {-1, 0, 0, 2, 1, 16.0f, 0.0f, -8.0f, 0.0f},
	    {-1, 0, 0, 1, 1, 1.0f, 0.0f, 0.0f, 0.0f},
	    {-1, 0, 0, 0, 1, -58.0f, -0.1f, 32.0f, 0.0f},
	    {-2, 0, 2, 0, 1, 46.0f, 0.0f, -24.0f, 0.0f},
	    {-2, 0, 0, 2, 1, -6.0f, 0.0f, 3.0f, 0.0f},
	    {-2, 0, 0, 0, 1, -2.0f, 0.0f, 1.0f, 0.0f},
	    {2, 0, -2, 0, 1, 1.0f, 0.0f, 0.0f, 0.0f},
	    {2, 0, 2, -2, 1, 1.0f, 0.0f, -1.0f, 0.0f},
	    {2, 0, 2, 0, 1, -5.0f, 0.0f, 3.0f, 0.0f},
	    {2, 0, 0, -2, 1, 4.0f, 0.0f, -2.0f, 0.0f},
	    {2, 0, 0, 0, 1, 2.0f, 0.0f, -1.0f, 0.0f},
	    {1, 1, 0, -2, 1, -1.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 2, -2, 1, 6.0f, 0.0f, -3.0f, 0.0f},
	    {1, 0, 2, 2, 1, -1.0f, 0.0f, 1.0f, 0.0f},
	    {1, 0, 2, 0, 1, -51.0f, 0.0f, 27.0f, 0.0f},
	    {1, 0, 0, -2, 1, -13.0f, 0.0f, 7.0f, 0.0f},
	    {1, 0, 0, 2, 1, -1.0f, 0.0f, 0.0f, 0.0f},
	    {1, 0, 0, 0, 1, 63.0f, 0.1f, -33.0f, 0.0f},
	    {0, -1, 2, -2, 1, -5.0f, 0.0f, 3.0f, 0.0f},
	    {0, -1, 2, 0, 1, -1.0f, 0.0f, 0.0f, 0.0f},
	    {0, -1, 0, 0, 1, -12.0f, 0.0f, 6.0f, 0.0f},
	    {0, 1, 2, -2, 1, 4.0f, 0.0f, -2.0f, 0.0f},
	    {0, 1, 2, 0, 1, 1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 1, 0, 0, 1, -15.0f, 0.0f, 9.0f, 0.0f},
	    {0, 0, -2, 2, 1, 1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, -2, 0, 1, -1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, 2, -2, 1, 129.0f, 0.1f, -70.0f, 0.0f},
	    {0, 0, 2, 2, 1, -7.0f, 0.0f, 3.0f, 0.0f},
	    {0, 0, 2, 0, 1, -386.0f, -0.4f, 200.0f, 0.0f},
	    {0, 0, 0, -2, 1, -5.0f, 0.0f, 3.0f, 0.0f},
	    {0, 0, 0, 2, 1, -6.0f, 0.0f, 3.0f, 0.0f},
	    {0, 0, 0, 0, 1, 0.0f, 0.0f, 0.0f, 0.0f},
	    {-1, -1, 2, 2, 2, -3.0f, 0.0f, 1.0f, 0.0f},
	    {-1, 0, 4, 0, 2, 1.0f, 0.0f, 0.0f, 0.0f},
	    {-1, 0, 2, 4, 2, -2.0f, 0.0f, 1.0f, 0.0f},
	    {-1, 0, 2, 2, 2, -59.0f, 0.0f, 26.0f, 0.0f},
	    {-1, 0, 2, 0, 2, 123.0f, 0.0f, -53.0f, 0.0f},
	    {-1, 0, 0, 0, 2, 1.0f, 0.0f, -1.0f, 0.0f},
	    {-2, 0, 2, 4, 2, -1.0f, 0.0f, 1.0f, 0.0f},
	    {-2, 0, 2, 2, 2, 1.0f, 0.0f, -1.0f, 0.0f},
	    {-2, 0, 2, 0, 2, -3.0f, 0.0f, 1.0f, 0.0f},
	    {3, 0, 2, -2, 2, 1.0f, 0.0f, 0.0f, 0.0f},
	    {3, 0, 2, 0, 2, -3.0f, 0.0f, 1.0f, 0.0f},
	    {2, 0, 2, -2, 2, 6.0f, 0.0f, -3.0f, 0.0f},
	    {2, 0, 2, 2, 2, -1.0f, 0.0f, 0.0f, 0.0f},
	    {2, 0, 2, 0, 2, -31.0f, 0.0f, 13.0f, 0.0f},
	    {1, -1, 2, 0, 2, -3.0f, 0.0f, 1.0f, 0.0f},
	    {1, 1, 2, -2, 2, 1.0f, 0.0f, -1.0f, 0.0f},
	    {1, 1, 2, 0, 2, 2.0f, 0.0f, -1.0f, 0.0f},
	    {1, 0, 2, -2, 2, 29.0f, 0.0f, -12.0f, 0.0f},
	    {1, 0, 2, 2, 2, -8.0f, 0.0f, 3.0f, 0.0f},
	    {1, 0, 2, 0, 2, -301.0f, 0.0f, 129.0f, -0.1f},
	    {1, 0, 0, 0, 2, -2.0f, 0.0f, 1.0f, 0.0f},
	    {0, -1, 2, 2, 2, -3.0f, 0.0f, 1.0f, 0.0f},
	    {0, -1, 2, 0, 2, -7.0f, 0.0f, 3.0f, 0.0f},
	    {0, 1, 2, 0, 2, 7.0f, 0.0f, -3.0f, 0.0f},
	    {0, 1, 0, 0, 2, 1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, 4, -2, 2, 1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, 2, -1, 2, -1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, 2, 4, 2, -1.0f, 0.0f, 0.0f, 0.0f},
	    {0, 0, 2, 2, 2, -38.0f, 0.0f, 16.0f, 0.0f},
	    {0, 0, 2, 1, 2, 2.0f, 0.0f, -1.0f, 0.0f},
	    {0, 0, 2, 0, 2, -2274.0f, -0.2f, 977.0f, -0.5f},
	    {0, 0, 0, 0, 2, 2062.0f, 0.2f, -895.0f, 0.5f},
	    {0, 2, 0, 0, 0, 17.0f, -0.1f, 0.0f, 0.0f},
	    {0, 1, 0, 0, 0, 1426.0f, -3.4f, 54.0f, -0.1f},
	    {0, -1, 2, -2, 2, 217.0f, -0.5f, -95.0f, 0.3f},
	    {0, 2, 2, -2, 2, -16.0f, 0.1f, 7.0f, 0.0f},
	    {0, 1, 2, -2, 2, -517.0f, 1.2f, 224.0f, -0.6f},
	    {0, 0, 2, -2, 2, -13187.0f, -1.6f, 5736.0f, -3.1f},
	    {0, 0, 0, 0, 1, -171996.0f, -174.2f, 92025.0f, 8.9f},
	};

	void nutationAnglesIau1980Full(double t, double& dpsi, double& deps) noexcept
	{
		const auto argument = [t](double c0, double c1, double c2, double c3) {
			return std::fmod(((c3 * ASEC_TO_RAD * t + c2 * ASEC_TO_RAD) * t + c1 * ASEC_TO_RAD) * t + c0 * ASEC_TO_RAD,
			                 Constants::TWOPI);
		};
		const double l = argument(485866.733, 1717915922.633, 31.310, 0.064);
		const double lp = argument(1287099.804, 129596581.224, -0.577, -0.012);
		const double f = argument(335778.877, 1739527263.137, -13.257, 0.011);
		const double d = argument(1072261.307, 1602961601.328, -6.891, 0.019);
		const double om = argument(450160.280, -6962890.539, 7.455, 0.008);
		dpsi = deps = 0.0;
		for (const auto& term : kWahr) {
			const double arg =
			    std::fmod(term.l * l + term.lp * lp + term.f * f + term.d * d + term.om * om, Constants::TWOPI);
			dpsi += (term.psi + term.psiT * t) * std::sin(arg);
			deps += (term.eps + term.epsT * t) * std::cos(arg);
		}
		dpsi *= 1e-4;
		deps *= 1e-4;
	}

	inline void buildPrecessionMatrixIau1976(double t, double P[3][3], bool dynamical = false) noexcept
	{
		const double zeta = dynamical ? ASEC_TO_RAD * ((0.017998 * t + 0.30188) * t + 2306.2181) * t
		                              : Constants::DEG_TO_RAD * ((0.0000050 * t + 0.0000839) * t + 0.6406161) * t;
		const double z = dynamical ? ASEC_TO_RAD * ((0.018203 * t + 1.09468) * t + 2306.2181) * t
		                           : Constants::DEG_TO_RAD * ((0.0000051 * t + 0.0003041) * t + 0.6406161) * t;
		const double theta = dynamical ? ASEC_TO_RAD * ((-0.041833 * t - 0.42665) * t + 2004.3109) * t
		                               : Constants::DEG_TO_RAD * ((-0.0000116 * t - 0.0001185) * t + 0.5567530) * t;

		const double cZ = std::cos(zeta), sZ = std::sin(zeta);
		const double cz = std::cos(z), sz = std::sin(z);
		const double ct = std::cos(theta), st = std::sin(theta);

		P[0][0] = cz * ct * cZ - sz * sZ;
		P[0][1] = -cz * ct * sZ - sz * cZ;
		P[0][2] = -cz * st;
		P[1][0] = sz * ct * cZ + cz * sZ;
		P[1][1] = -sz * ct * sZ + cz * cZ;
		P[1][2] = -sz * st;
		P[2][0] = st * cZ;
		P[2][1] = -st * sZ;
		P[2][2] = ct;
	}

}  // namespace

namespace EarthRotation {
	namespace Iau2000B {

		double meanObliquity(double jd_tt)
		{
			const double t = juliancentury(jd_tt);
			const double eps0 = 84381.406;
			const double eps =
			    eps0 + t * (-46.836769 + t * (-0.0001831 + t * (0.00200340 + t * (-5.76e-7 + t * (-4.34e-8)))));
			return eps * ASEC_TO_RAD;
		}

		void nutationAngles(double jd_tt, double& dpsi, double& deps)
		{
			const double t = juliancentury(jd_tt);

			const double l = meanAnomalyMoon(t);
			const double lp = meanAnomalySun(t);
			const double F = meanArgLatMoon(t);
			const double D = meanElongMoon(t);
			const double Om = meanLongAscNode(t);

			double sumDpsi = 0.0;
			double sumDeps = 0.0;

			for (int i = kNTerms - 1; i >= 0; --i) {
				const auto& term = kLunisolar[i];
				const double arg = term.nl * l + term.nlp * lp + term.nF * F + term.nD * D + term.nOm * Om;
				const double sa = std::sin(arg);
				const double ca = std::cos(arg);
				sumDpsi += (term.Ap + term.Apt * t) * sa + term.Aps * ca;
				sumDeps += (term.Bp + term.Bpt * t) * ca + term.Bps * sa;
			}

			dpsi = sumDpsi * U2R + kDpsiPlanetary;
			deps = sumDeps * U2R + kDepsPlanetary;
		}

		double equationOfEquinoxes(double jd_tt)
		{
			double dpsi = 0.0, deps = 0.0;
			nutationAngles(jd_tt, dpsi, deps);
			const double epsA = meanObliquity(jd_tt);
			return dpsi * std::cos(epsA);
		}

		void applyNutation(double jd_tt, Vec3D& v, int direction)
		{
			double dpsi = 0.0, deps = 0.0;
			nutationAngles(jd_tt, dpsi, deps);
			const double epsA = meanObliquity(jd_tt);

			double N[3][3];
			buildNutationMatrix(epsA, dpsi, deps, N);

			if (direction < 0)
				transpose3(N);

			mulMatVec(N, v);
		}

		void applyPrecessionNutation(double jd_tt, Vec3D& v, int direction)
		{
			double P[3][3], N[3][3];
			buildPrecessionMatrix(jd_tt, P);

			double dpsi = 0.0, deps = 0.0;
			nutationAngles(jd_tt, dpsi, deps);
			const double epsA = meanObliquity(jd_tt);
			buildNutationMatrix(epsA, dpsi, deps, N);

			if (direction >= 0) {
				mulMatVec(P, v);
				mulMatVec(N, v);
			} else {
				transpose3(N);
				transpose3(P);
				mulMatVec(N, v);
				mulMatVec(P, v);
			}
		}

	}  // namespace Iau2000B

	namespace Iau1976_1980 {

		double gmstRadians(double jd_ut1) { return gmstRadiansMjd(jd_ut1 - 2400000.5); }

		double gmstRadiansMjd(double mjd)
		{
			constexpr double c0 = 24110.54841;
			constexpr double c1 = 8640184.812866;
			constexpr double c2 = 9.3104e-2;
			constexpr double c3 = -6.2e-6;
			constexpr double rap = 1.00273790934;

			const int day = static_cast<int>(mjd);
			const double t = (static_cast<double>(day) - 51544.5) / 36525.0;
			const double gmst0 = ((c3 * t + c2) * t + c1) * t + c0;
			return Helpers::wrapTwoPi(gmst0 * Constants::TWOPI / Constants::SECONDS_PER_DAY +
			                          (mjd - static_cast<double>(day)) * Constants::TWOPI * rap);
		}

		double equationOfEquinoxes(double jd_tt)
		{
			double dpsiArcsec = 0.0;
			double depsArcsec = 0.0;
			nutationAnglesIau1980(jd_tt, dpsiArcsec, depsArcsec);
			return dpsiArcsec * ASEC_TO_RAD * std::cos(meanObliquityIau1976(jd_tt));
		}

		void applyPrecessionNutation(double jd_tt, Vec3D& v, int direction)
		{
			double dpsiArcsec = 0.0;
			double depsArcsec = 0.0;
			nutationAnglesIau1980(jd_tt, dpsiArcsec, depsArcsec);

			double N[3][3];
			buildNutationMatrix(meanObliquityIau1976(jd_tt), dpsiArcsec * ASEC_TO_RAD, depsArcsec * ASEC_TO_RAD, N);
			double P[3][3];
			buildPrecessionMatrixIau1976(juliancentury(jd_tt), P);

			if (direction >= 0) {
				mulMatVec(P, v);
				mulMatVec(N, v);
			} else {
				transpose3(N);
				transpose3(P);
				mulMatVec(N, v);
				mulMatVec(P, v);
			}
		}

		Vec3D dynamicalPoleEquatorial(double jdTt) noexcept
		{
			double dpsiArcsec = 0.0, depsArcsec = 0.0;
			nutationAnglesIau1980Full(juliancentury(jdTt), dpsiArcsec, depsArcsec);
			double N[3][3], P[3][3];
			const double t = juliancentury(jdTt);
			const double eps = (((0.001813 * t - 0.00059) * t - 46.815) * t + 84381.448) * ASEC_TO_RAD;
			buildNutationMatrix(eps, dpsiArcsec * ASEC_TO_RAD, depsArcsec * ASEC_TO_RAD, N);
			buildPrecessionMatrixIau1976(juliancentury(jdTt), P, true);
			Vec3D pole(N[2][0], N[2][1], N[2][2]);
			transpose3(P);
			mulMatVec(P, pole);
			return pole;
		}

	}  // namespace Iau1976_1980

	namespace WithEop {
		namespace {

			Vec3D poleEquatorialIau2006Mjd(double mjdTt, double xpRad, double ypRad, double dxRad, double dyRad,
			                               double ttMinusUt1Seconds) noexcept
			{
				double x, y, s;
				eraXys06a(2400000.5, mjdTt, &x, &y, &s);
				x += dxRad;
				y += dyRad;
				double celestialToIntermediate[3][3], polarMotion[3][3], celestialToTerrestrial[3][3];
				eraC2ixys(x, y, s, celestialToIntermediate);
				const double day = std::floor(mjdTt);
				const double ut1Fraction = (mjdTt - day) - ttMinusUt1Seconds / 86400.0;
				const double earthRotation = eraEra00(2400000.5 + day, ut1Fraction);
				eraPom00(xpRad, ypRad, eraSp00(2400000.5, mjdTt), polarMotion);
				eraC2tcio(celestialToIntermediate, earthRotation, polarMotion, celestialToTerrestrial);
				return {celestialToTerrestrial[2][0], celestialToTerrestrial[2][1], celestialToTerrestrial[2][2]};
			}

			Vec3D poleEquatorialIau1980Mjd(double mjdTt, double xpRad, double ypRad, double dpsiRad, double depsRad,
			                               double ttMinusUt1Seconds) noexcept
			{
				const double t = (mjdTt - 51544.5) / 36525.0;
				double nominalPsi = 0.0, nominalEps = 0.0;
				nutationAnglesIau1980Full(t, nominalPsi, nominalEps);
				const double psi = nominalPsi * ASEC_TO_RAD + dpsiRad;
				const double deps = nominalEps * ASEC_TO_RAD + depsRad;
				const double eps = (((0.001813 * t - 0.00059) * t - 46.815) * t + 84381.448) * ASEC_TO_RAD;
				double N[3][3], P[3][3];
				buildNutationMatrix(eps, psi, deps, N);
				buildPrecessionMatrixIau1976(t, P, true);

				double day = std::floor(mjdTt);
				double seconds = (mjdTt - day) * 86400.0 - ttMinusUt1Seconds;
				const double shift = std::floor(seconds / 86400.0);
				day += shift;
				seconds -= shift * 86400.0;
				constexpr double secondsToRad = Constants::PI / (12.0 * 3600.0);
				constexpr double a = 24110.54841 * secondsToRad;
				constexpr double b = 8640184.812866 * secondsToRad;
				constexpr double c = 0.093104 * secondsToRad;
				constexpr double d = -6.2e-6 * secondsToRad;
				const double t0 = (day - 51544.5) / 36525.0;
				const double midpoint = (day + (seconds / 86400.0) / 2.0 - 51544.5) / 36525.0;
				const double rate =
				    ((3.0 * d * midpoint + 2.0 * c) * midpoint + b) / (36525.0 * 86400.0) + Constants::TWOPI / 86400.0;
				const double gmst = Helpers::wrapTwoPi(((d * t0 + c) * t0 + b) * t0 + a + rate * seconds);
				const double gast = gmst + psi * std::cos(eps);

				const double sx = std::sin(xpRad), cx = std::cos(xpRad);
				const double sy = std::sin(ypRad), cy = std::cos(ypRad);
				const double sg = std::sin(gast), cg = std::cos(gast);
				Vec3D pole(-sx * cg - cx * sy * sg, -sx * sg + cx * sy * cg, cx * cy);
				transpose3(N);
				transpose3(P);
				mulMatVec(N, pole);
				mulMatVec(P, pole);
				return pole;
			}

		}  // namespace

		Vec3D poleEquatorialMjd(double mjdTt, const orbital_core::EarthOrientationSample& eop) noexcept
		{
			if (const auto* offsets = std::get_if<orbital_core::Iau1980PoleOffsets>(&eop.poleOffsets))
				return poleEquatorialIau1980Mjd(mjdTt, eop.xpRad, eop.ypRad, offsets->dpsiRad, offsets->depsRad,
				                                eop.ttMinusUt1Seconds);
			const auto& offsets = std::get<orbital_core::Iau2000PoleOffsets>(eop.poleOffsets);
			return poleEquatorialIau2006Mjd(mjdTt, eop.xpRad, eop.ypRad, offsets.dxRad, offsets.dyRad,
			                                eop.ttMinusUt1Seconds);
		}

		Vec3D terrestrialToEquatorialMjd(double mjdTt, const Vec3D& bodyFixed,
		                                 const orbital_core::EarthOrientationSample& eop) noexcept
		{
			if (const auto* offsets = std::get_if<orbital_core::Iau1980PoleOffsets>(&eop.poleOffsets)) {
				const double t = (mjdTt - 51544.5) / 36525.0;
				double psi = 0.0, deps = 0.0;
				nutationAnglesIau1980Full(t, psi, deps);
				psi = psi * ASEC_TO_RAD + offsets->dpsiRad;
				deps = deps * ASEC_TO_RAD + offsets->depsRad;
				const double eps = (((0.001813 * t - 0.00059) * t - 46.815) * t + 84381.448) * ASEC_TO_RAD;
				double n[3][3], p[3][3];
				buildNutationMatrix(eps, psi, deps, n);
				buildPrecessionMatrixIau1976(t, p, true);
				const double gast =
				    Iau1976_1980::gmstRadiansMjd(mjdTt - eop.ttMinusUt1Seconds / 86400.0) + psi * std::cos(eps);
				const double sx = std::sin(eop.xpRad), cx = std::cos(eop.xpRad);
				const double sy = std::sin(eop.ypRad), cy = std::cos(eop.ypRad);
				const Vec3D polar(cx * bodyFixed.x - sx * bodyFixed.z,
				                  sy * sx * bodyFixed.x + cy * bodyFixed.y + sy * cx * bodyFixed.z,
				                  cy * sx * bodyFixed.x - sy * bodyFixed.y + cy * cx * bodyFixed.z);
				Vec3D out(std::cos(gast) * polar.x - std::sin(gast) * polar.y,
				          std::sin(gast) * polar.x + std::cos(gast) * polar.y, polar.z);
				transpose3(n);
				transpose3(p);
				mulMatVec(n, out);
				mulMatVec(p, out);
				return out;
			}
			const auto& offsets = std::get<orbital_core::Iau2000PoleOffsets>(eop.poleOffsets);
			double x, y, s;
			eraXys06a(2400000.5, mjdTt, &x, &y, &s);
			double c2i[3][3], polar[3][3], c2t[3][3];
			eraC2ixys(x + offsets.dxRad, y + offsets.dyRad, s, c2i);
			const double day = std::floor(mjdTt);
			const double era = eraEra00(2400000.5 + day, mjdTt - day - eop.ttMinusUt1Seconds / 86400.0);
			eraPom00(eop.xpRad, eop.ypRad, eraSp00(2400000.5, mjdTt), polar);
			eraC2tcio(c2i, era, polar, c2t);
			transpose3(c2t);
			Vec3D out = bodyFixed;
			mulMatVec(c2t, out);
			return out;
		}

	}  // namespace WithEop

}  // namespace EarthRotation
