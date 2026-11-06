#include <gtest/gtest.h>
#include "obs/LightDeflection.h"
#include "obs/AstrometricResidual.h"

#include <cmath>

namespace {
	using orbital_core::obs::catalogReferencedDirection;
	using orbital_core::obs::deflectLight;

	constexpr double kMasPerRadian = 180.0 / 3.14159265358979323846 * 3600.0e3;

	struct SkyShift {
		double east{};
		double north{};
	};

	SkyShift skyShiftMas(const Vec3D& from, const Vec3D& to)
	{
		const Vec3D p = from / from.magnitude();
		const Vec3D d = to / to.magnitude() - p;
		const double ra = std::atan2(p.y, p.x);
		const double dec = std::asin(p.z);
		const Vec3D east(-std::sin(ra), std::cos(ra), 0.0);
		const Vec3D north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
		return {d.dot(east) * kMasPerRadian, d.dot(north) * kMasPerRadian};
	}
}  // namespace

TEST(LightDeflection, MatchesSofaReferenceCase)
{
	const Vec3D p(-0.763276255, -0.608633767, -0.216735543);
	const Vec3D e(0.76700421, 0.605629598, 0.211937094);
	const Vec3D deflected = deflectLight(0.00028574, p, p, e, 8.91276983, 3e-10);
	EXPECT_NEAR(deflected.x, -0.7632762548968159627, 1e-12);
	EXPECT_NEAR(deflected.y, -0.6086337670823762701, 1e-12);
	EXPECT_NEAR(deflected.z, -0.2167355431320546947, 1e-12);
}

TEST(LightDeflection, DirectionJacobianMatchesIndependentDifferences)
{
	orbital_core::obs::PreparedObservation observation;
	observation.earthR = Vec3D(0.9, -0.4, 0.05);
	observation.observerOffsetEclAu = Vec3D(0.00002, 0.00001, -0.00001);
	for (const auto convention : {orbital_core::obs::AstrometricConvention::CatalogReferenced,
	                              orbital_core::obs::AstrometricConvention::Absolute}) {
		observation.astrometricConvention = convention;
		for (const Vec3D source : {Vec3D(1.7, 1.0, 0.2), Vec3D(-1.4, 1.2, -0.4), Vec3D(0.91, -0.35, 0.052)}) {
			const Eigen::Matrix3d analytic =
			    orbital_core::obs::opticalObservationDirectionJacobian(source, observation);
			constexpr double h = 1.0e-4;
			for (int col = 0; col < 3; ++col) {
				Vec3D plus = source, minus = source;
				plus[col] += h;
				minus[col] -= h;
				const Vec3D p = orbital_core::obs::opticalObservationDirection(plus, observation);
				const Vec3D m = orbital_core::obs::opticalObservationDirection(minus, observation);
				for (int row = 0; row < 3; ++row)
					EXPECT_NEAR(analytic(row, col), (p[row] - m[row]) / (2 * h), 2.0e-10);
			}
		}
	}
}

TEST(LightDeflection, GaiaPhotocenterDerivativeMatchesIndependentDifferences)
{
	const Vec3D observer(0.9, -0.4, 0.05);
	const Vec3D source(1.7, 1.0, 0.2);
	constexpr double radiusAu = 1.0 / 149597870.7;
	const auto analytic = orbital_core::obs::gaiaPhotocenterCorrection(source, observer, radiusAu);
	ASSERT_GT(analytic.offset.magnitude(), 0.0);
	constexpr double h = 1.0e-5;
	for (int col = 0; col < 3; ++col) {
		Vec3D plus = source, minus = source;
		plus[col] += h;
		minus[col] -= h;
		const auto p = orbital_core::obs::gaiaPhotocenterCorrection(plus, observer, radiusAu).offset;
		const auto m = orbital_core::obs::gaiaPhotocenterCorrection(minus, observer, radiusAu).offset;
		for (int row = 0; row < 3; ++row)
			EXPECT_NEAR(analytic.derivative(row, col), (p[row] - m[row]) / (2 * h), 1.0e-17);
	}
}

TEST(LightDeflection, GaiaPhotocenterUncertaintyPreservesCovariance)
{
	orbital_core::obs::PreparedObservation observation;
	observation.astrometricConvention = orbital_core::obs::AstrometricConvention::GaiaFpr;
	observation.earthR = Vec3D(1.0, 0.0, 0.0);
	observation.sigmaRa = 0.002 * Constants::ARCSEC_RAD;
	observation.sigmaDec = 0.003 * Constants::ARCSEC_RAD;
	observation.sigmaCorr = 0.5;
	observation.photocenterRadiusAu = 1.0 / 149597870.7;
	const Vec3D source(2.0, 0.0, 0.0);
	const auto transform = orbital_core::obs::astrometricResidualTransform(observation, source);
	ASSERT_TRUE(transform.valid);
	const auto reported = orbital_core::obs::residualWeightForReport(observation, 0.0, 0.0, true, transform);
	const double extraArcsec = 0.3 * observation.photocenterRadiusAu / Constants::ARCSEC_RAD;
	EXPECT_NEAR(reported.sigmaRaArcsec, std::hypot(0.002, extraArcsec), 1e-14);
	EXPECT_NEAR(reported.sigmaDecArcsec, std::hypot(0.003, extraArcsec), 1e-14);
}

TEST(LightDeflection, StarAtQuadratureIsDeflectedFourMasAwayFromTheSun)
{
	const Vec3D p(0.0, 1.0, 0.0);
	const Vec3D e(1.0, 0.0, 0.0);
	const Vec3D shift = deflectLight(1.0, p, p, e, 1.0, 1e-6) - p;
	EXPECT_NEAR(shift.magnitude() * kMasPerRadian, 4.0719, 1e-3);
	EXPECT_GT(shift.x, 0.0);
}

TEST(LightDeflection, NearbySourceIsReferencedToTheStarDeflection)
{
	const Vec3D observer(1.0, 0.0, 0.0);
	const Vec3D source = observer + Vec3D(0.0, 0.01, 0.0);
	const Vec3D topocentric = source - observer;
	const Vec3D referenced = catalogReferencedDirection(topocentric, source, observer);
	EXPECT_NEAR(referenced.magnitude(), topocentric.magnitude(), 1e-15);
	const Vec3D shift = referenced / referenced.magnitude() - Vec3D(0.0, 1.0, 0.0);
	EXPECT_LT(shift.x, 0.0);
	EXPECT_NEAR(shift.magnitude() * kMasPerRadian, 4.07, 0.05);
}

TEST(LightDeflection, InfinitelyDistantSourceNeedsNoCorrection)
{
	const Vec3D observer(1.0, 0.0, 0.0);
	const Vec3D direction = Vec3D(0.3, 0.8, 0.52) / Vec3D(0.3, 0.8, 0.52).magnitude();
	const Vec3D source = observer + direction * 1.0e9;
	const Vec3D referenced = catalogReferencedDirection(source - observer, source, observer);
	const SkyShift shift = skyShiftMas(source - observer, referenced);
	EXPECT_NEAR(shift.east, 0.0, 1e-6);
	EXPECT_NEAR(shift.north, 0.0, 1e-6);
}

TEST(LightDeflection, ApophisOccultationOf2022April9)
{
	const Vec3D earth(-9.458715200441830E-01, -3.019114895218548E-01, -1.308738403008057E-01);
	const Vec3D apophis(-5.732508271781787E-01, -6.660725999660179E-01, -2.620123332444549E-01);
	const Vec3D referenced = catalogReferencedDirection(apophis - earth, apophis, earth);
	const SkyShift shift = skyShiftMas(apophis - earth, referenced);
	EXPECT_NEAR(shift.east, 4.996886, 1e-4);
	EXPECT_NEAR(shift.north, 1.370388, 1e-4);
}

TEST(LightDeflection, DegenerateGeometryIsReturnedUnchanged)
{
	const Vec3D zero;
	const Vec3D observer(1.0, 0.0, 0.0);
	const Vec3D source(1.0, 0.0, 0.0);
	const Vec3D result = catalogReferencedDirection(zero, source, observer);
	EXPECT_EQ(result.x, 0.0);
	EXPECT_EQ(result.y, 0.0);
	EXPECT_EQ(result.z, 0.0);
}
