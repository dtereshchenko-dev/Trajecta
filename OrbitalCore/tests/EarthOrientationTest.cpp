#include <gtest/gtest.h>

#include "earth_orientation/EarthOrientation.h"
#include "earth_orientation/EarthPck.h"
#include "earth_orientation/EarthPckPoleCache.h"
#include "earth_orientation/EarthRotation.h"
#include "obs/ObserverGeometry.h"
#include "orbital_core/Constants.h"
#include "orbital_core/RuntimeData.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/BodyFixedFrame.h"

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <future>
#include <limits>
#include <vector>
#include <SpiceUsr.h>

using namespace orbital_core;

TEST(EarthOrientationTest, CachedPckPoleMatchesSpiceAcrossCoverageAndThreads)
{
	const Vec3D axis(0.0, 0.0, 1.0);
	EarthPckPoleCache cache;
	for (const auto* name : {"earth_predict.bpc", "earth_historic.bpc", "earth_latest.bpc"}) {
		if (earthPckLoaded(name))
			ASSERT_TRUE(cache.addFile(earthPckFilePath(name)));
	}
	std::vector<double> epochs;
	std::vector<Vec3D> reference;
	for (double epoch = 39000.0; epoch < 88000.0; epoch += 7.25) {
		Vec3D pole, velocity;
		const bool expected = earthPckState(epoch, axis, pole, velocity);
		Vec3D cached;
		EXPECT_EQ(earthPckPole(epoch, cached), expected) << epoch;
		if (expected) {
			ASSERT_TRUE(cache.pole(epoch, cached)) << epoch;
			EXPECT_LT((cached - pole).magnitude(), 1e-14) << epoch;
			epochs.push_back(epoch);
			reference.push_back(pole);
		}
	}
	if (epochs.empty())
		GTEST_SKIP() << "Earth PCK unavailable";
	std::vector<double> edges;
	for (const auto* name : {"earth_predict.bpc", "earth_historic.bpc", "earth_latest.bpc"}) {
		if (!earthPckLoaded(name))
			continue;
		SpiceInt handle = 0;
		dafopr_c(earthPckFilePath(name).string().c_str(), &handle);
		ASSERT_FALSE(failed_c());
		dafbfs_c(handle);
		SpiceBoolean found = SPICEFALSE;
		while (true) {
			daffna_c(&found);
			ASSERT_FALSE(failed_c());
			if (!found)
				break;
			double descriptor[5]{}, bounds[2]{};
			SpiceInt fields[5]{};
			dafgs_c(descriptor);
			dafus_c(descriptor, 2, 5, bounds, fields);
			ASSERT_FALSE(failed_c());
			if (fields[0] == 3000) {
				edges.push_back(bounds[0]);
				edges.push_back(bounds[1]);
			}
		}
		dafcls_c(handle);
		ASSERT_FALSE(failed_c());
	}
	for (const double edge : edges)
		for (const double offset : {-1.0, -1e-6, 0.0, 1e-6, 1.0}) {
			const double epoch = 51544.5 + (edge + offset) / Constants::SECONDS_PER_DAY;
			Vec3D expected, velocity, actual;
			const bool valid = earthPckState(epoch, axis, expected, velocity);
			EXPECT_EQ(cache.pole(epoch, actual), valid) << epoch;
			if (valid)
				EXPECT_LT((actual - expected).magnitude(), 1e-14) << epoch;
		}
	std::vector<std::future<double>> workers;
	for (int worker = 0; worker < 8; ++worker) {
		workers.push_back(std::async(std::launch::async, [&cache, &epochs, &reference] {
			double maxError = 0.0;
			for (std::size_t i = 0; i < epochs.size(); ++i) {
				Vec3D pole;
				if (!cache.pole(epochs[i], pole))
					return std::numeric_limits<double>::infinity();
				maxError = std::max(maxError, (pole - reference[i]).magnitude());
			}
			return maxError;
		}));
	}
	for (auto& worker : workers)
		EXPECT_LT(worker.get(), 1e-14);
}

TEST(EarthOrientationTest, PckPoleCacheHandlesMissingAndInvalidInput)
{
	EarthPckPoleCache empty;
	Vec3D pole;
	EXPECT_FALSE(empty.pole(60000.0, pole));
	EXPECT_FALSE(empty.pole(std::numeric_limits<double>::quiet_NaN(), pole));
	(void)earthPckLoaded("earth_latest.bpc");
	const auto missing = std::filesystem::temp_directory_path() / "trajecta-missing-pole-cache.bpc";
	ASSERT_FALSE(std::filesystem::exists(missing));
	EXPECT_FALSE(empty.addFile(missing));
	EXPECT_FALSE(empty.pole(60000.0, pole));
}

TEST(EarthOrientationTest, PckPoleCacheRespectsUnsupportedSegmentPriority)
{
	(void)earthPckLoaded("earth_latest.bpc");
	const auto directory = std::filesystem::temp_directory_path();
	const auto supported = directory / "trajecta-test-pole-ecliptic.bpc";
	const auto unsupported = directory / "trajecta-test-pole-equatorial.bpc";
	ASSERT_FALSE(std::filesystem::exists(supported));
	ASSERT_FALSE(std::filesystem::exists(unsupported));
	const auto create = [](const std::filesystem::path& path, const char* frame) {
		SpiceInt handle = 0;
		pckopn_c(path.string().c_str(), "test", 0, &handle);
		if (failed_c()) {
			reset_c();
			return false;
		}
		double coefficients[] = {0.3, 0.0, 0.4, 0.0, 0.5, 0.0};
		pckw02_c(handle, 3000, frame, 0.0, 86400.0, "test", 86400.0, 1, 1, coefficients, 0.0);
		const bool success = !failed_c();
		if (!success)
			reset_c();
		pckcls_c(handle);
		const bool closed = !failed_c();
		if (!closed)
			reset_c();
		return success && closed;
	};
	ASSERT_TRUE(create(supported, "ECLIPJ2000"));
	ASSERT_TRUE(create(unsupported, "J2000"));
	EarthPckPoleCache cache;
	EXPECT_TRUE(cache.addFile(supported));
	Vec3D pole;
	EXPECT_TRUE(cache.pole(51544.75, pole));
	EXPECT_FALSE(cache.pole(51546.0, pole));
	EXPECT_TRUE(cache.addFile(unsupported));
	EXPECT_FALSE(cache.pole(51544.75, pole));
	EXPECT_TRUE(std::filesystem::remove(supported));
	EXPECT_TRUE(std::filesystem::remove(unsupported));
}

TEST(EarthOrientationTest, PckForwardAndInverseAgree)
{
	constexpr double mjdTdb = 53398.0306038786;
	const Vec3D fixed(0.4, -0.6, 0.2);
	Vec3D inertial, velocity, recovered;
	if (!earthPckState(mjdTdb, fixed, inertial, velocity))
		GTEST_SKIP() << "Earth PCK unavailable at test epoch";
	ASSERT_TRUE(earthPckToTerrestrial(mjdTdb, inertial, recovered));
	EXPECT_LT((recovered - fixed).magnitude(), 1e-14);
	Vec3D recoveredVelocity;
	ASSERT_TRUE(earthPckToTerrestrialState(mjdTdb, inertial, velocity * (1.0 / 86400.0), recovered, recoveredVelocity));
	EXPECT_LT((recovered - fixed).magnitude(), 1e-14);
	EXPECT_LT(recoveredVelocity.magnitude(), 1e-18);
	const auto fixedVelocity = hazard::bodyFixedVelocityMjd(hazard::Planet::Earth, Helpers::ecvToEcl(inertial),
	                                                        Helpers::ecvToEcl(velocity) * (1.0 / 86400.0), mjdTdb);
	EXPECT_LT(fixedVelocity.magnitude(), 1e-18);
	const auto ecliptic = Helpers::ecvToEcl(inertial);
	const auto body = hazard::icrfToBodyFixedMjd(hazard::Planet::Earth, ecliptic, mjdTdb);
	EXPECT_LT((body - fixed).magnitude(), 1e-14);
	const auto transform = hazard::bodyFixedTransformMjd(hazard::Planet::Earth, mjdTdb);
	ASSERT_TRUE(transform.valid);
	EXPECT_EQ(transform.frame, "ITRF93");
}

TEST(EarthOrientationTest, PckValidationAndCoverageBoundaries)
{
	const auto directory = getExecutableDir() / "earth_orientation";
	const auto historic = directory / "earth_historic.bpc";
	if (!std::filesystem::exists(historic))
		GTEST_SKIP() << "Earth PCK unavailable";
	EXPECT_TRUE(validateEarthPck(historic));
	EXPECT_FALSE(validateEarthPck(directory / "missing-test.bpc"));
	Vec3D position, velocity;
	const Vec3D fixed(0.5, 0.3, 0.4);
	EXPECT_FALSE(earthPckState(37683.0, fixed, position, velocity));
	EXPECT_TRUE(earthPckState(37685.0, fixed, position, velocity));
	EXPECT_TRUE(earthPckState(61446.0, fixed, position, velocity));
	EXPECT_FALSE(earthPckState(234166.0, fixed, position, velocity));
	const obs::ObserverData station{10.0, 0.8, 0.6, "test"};
	EXPECT_THROW((void)obs::groundObserverStateEclAuMjd(234166.0, 234166.0, station), std::runtime_error);
	EXPECT_FALSE(hazard::bodyFixedTransformMjd(hazard::Planet::Earth, 234166.0).valid);
}

TEST(EarthOrientationTest, C04ClockDoesNotExtrapolateOutsideObservedCoverage)
{
	const auto data = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
	const EarthOrientation table(data);
	if (!table.isAvailable())
		GTEST_SKIP() << "C04 data unavailable";
	constexpr double mjdTt = 59169.5;
	ASSERT_TRUE(table.coversObservedMjdTt(mjdTt));
	double ut1 = 0.0, seconds = 0.0;
	ASSERT_TRUE(table.ttToUt1Mjd(mjdTt, ut1));
	ASSERT_TRUE(table.ut1ToTtSeconds(ut1, seconds));
	EXPECT_NEAR(ut1 + seconds / 86400.0, mjdTt, 1e-11);
	const double dailyTt = 59169.0 + 69.184 / 86400.0;
	ASSERT_TRUE(table.ttToUt1Mjd(dailyTt, ut1));
	EXPECT_NEAR((dailyTt - ut1) * 86400.0, 69.184 + 0.1790137, 1e-6);
	double unchanged = 123.0;
	EXPECT_FALSE(table.ttToUt1Mjd(100000.0, unchanged));
	EXPECT_DOUBLE_EQ(unchanged, 123.0);
}

TEST(EarthOrientationTest, OffsetTypeSelectsTheRotationWithoutImplicitConversion)
{
	constexpr double epoch = 59169.5;
	EarthOrientationSample base;
	base.poleOffsets = Iau2000PoleOffsets{};
	const auto p0 = EarthRotation::WithEop::poleEquatorialMjd(epoch, base);
	auto shifted = base;
	constexpr double dx = 1e-8, dy = -2e-8;
	shifted.poleOffsets = Iau2000PoleOffsets{dx, dy};
	const auto p1 = EarthRotation::WithEop::poleEquatorialMjd(epoch, shifted);
	EXPECT_NEAR(p1.x - p0.x, dx, 1e-15);
	EXPECT_NEAR(p1.y - p0.y, dy, 1e-15);
	EXPECT_NEAR(p1.magnitude(), 1., 1e-14);
	shifted.poleOffsets = Iau1980PoleOffsets{dx, dy};
	const auto other = EarthRotation::WithEop::poleEquatorialMjd(epoch, shifted);
	EXPECT_GT((other - p1).magnitude(), 1e-9);
}

TEST(EarthOrientationTest, MissingFallbackDoesNotInventObserverOrientation)
{
	const auto missing = std::filesystem::temp_directory_path() / "trajecta-missing-orientation-test";
	ASSERT_FALSE(std::filesystem::exists(missing));
	const EarthOrientation orientation(missing);
	const obs::ObserverData station{10.0, 0.8, 0.6, "test"};
	EXPECT_THROW((void)obs::groundObserverStateEclAuMjd(59169.5, 59169.5, station, orientation), std::runtime_error);
}

TEST(EarthOrientationTest, MissingPckAndTimeDataReportsUnavailable)
{
	if (earthPckLoaded("earth_latest.bpc") || earthPckLoaded("earth_historic.bpc") ||
	    earthPckLoaded("earth_predict.bpc"))
		GTEST_SKIP() << "Run in a directory without Earth orientation data";
	const obs::ObserverData station{10.0, 0.8, 0.6, "test"};
	constexpr double epoch = 61446.0;
	EXPECT_THROW((void)obs::groundObserverStateEclAuMjd(epoch, epoch, station), std::runtime_error);
	EXPECT_THROW((void)obs::radarStationStateEclAuMjd(epoch, epoch, station), std::runtime_error);
	EXPECT_FALSE(hazard::bodyFixedTransformMjd(hazard::Planet::Earth, epoch).valid);
}

TEST(EarthOrientationTest, MissingRecentPckUsesRemainingCoverage)
{
	if (earthPckLoaded("earth_latest.bpc") || !earthPckLoaded("earth_historic.bpc") ||
	    !earthPckLoaded("earth_predict.bpc"))
		GTEST_SKIP() << "Requires historic and predictive PCK without recent PCK";
	Vec3D position, velocity;
	const Vec3D fixed(0.5, 0.3, 0.4);
	EXPECT_TRUE(earthPckState(61041.0, fixed, position, velocity));
	EXPECT_TRUE(earthPckState(62502.0, fixed, position, velocity));
}

TEST(EarthOrientationTest, MissingPredictivePckReportsUncoveredFuture)
{
	if (earthPckLoaded("earth_predict.bpc") || !earthPckLoaded("earth_latest.bpc"))
		GTEST_SKIP() << "Requires recent PCK without predictive PCK";
	constexpr double epoch = 62502.0;
	const obs::ObserverData station{10.0, 0.8, 0.6, "test"};
	EXPECT_THROW((void)obs::groundObserverStateEclAuMjd(epoch, epoch, station), std::runtime_error);
	EXPECT_FALSE(hazard::bodyFixedTransformMjd(hazard::Planet::Earth, epoch).valid);
}
