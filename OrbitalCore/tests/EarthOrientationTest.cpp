#include <gtest/gtest.h>

#include "earth_orientation/EarthOrientation.h"
#include "earth_orientation/EarthRotation.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/ObsNode.h"
#include "orbital_core/Types.h"
#include "orbital_core/hazard/BodyFixedFrame.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace {
	using namespace orbital_core;

	double tt(double mjdUtc) { return mjdUtc + Convertors::utcToTdtSeconds(mjdUtc + 2400000.5) / 86400.0; }

	std::string record(double mjd, double xp, double yp, double dut1, double dx = 1.0, double dy = -2.0,
	                   char flag = 'I', bool offsets = true)
	{
		std::string line(185, ' ');
		const auto field = [&](int start, int width, int precision, double value) {
			std::ostringstream out;
			out.imbue(std::locale::classic());
			out << std::fixed << std::setprecision(precision) << std::setw(width) << value;
			line.replace(start, width, out.str());
		};
		field(7, 8, 2, mjd);
		field(18, 9, 6, xp);
		field(37, 9, 6, yp);
		field(58, 10, 7, dut1);
		line[16] = line[57] = flag;
		if (offsets) {
			line[95] = flag;
			field(97, 9, 3, dx);
			field(116, 9, 3, dy);
		}
		return line;
	}

	class EarthOrientationTest : public ::testing::Test {
	   protected:
		void SetUp() override
		{
			directory = std::filesystem::temp_directory_path() /
			            ("orbital-eop-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
			ASSERT_TRUE(std::filesystem::create_directory(directory));
			ASSERT_TRUE(std::filesystem::create_directory(directory / "earth_orientation"));
		}

		void TearDown() override
		{
			std::filesystem::remove(directory / "earth_orientation/finals2000A.all");
			std::filesystem::remove(directory / "earth_orientation");
			std::filesystem::remove(directory);
		}

		void write(std::initializer_list<std::string> lines)
		{
			std::ofstream file(directory / "earth_orientation/finals2000A.all");
			ASSERT_TRUE(file.is_open());
			for (const auto& line : lines)
				file << line << '\n';
		}

		std::filesystem::path directory;
	};

	TEST_F(EarthOrientationTest, FinalsUnitsAndExactEndpoints)
	{
		write({record(59169, .1, -.2, -.175), record(59170, .3, -.4, -.177, 3., -4., 'P')});
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		ASSERT_TRUE(table.isAvailable());
		EXPECT_EQ(table.size(), 2u);
		EXPECT_EQ(table.profile(), EarthOrientationProfile::IAU2000);
		const auto first = table.lookupMjdTt(tt(59169));
		const auto last = table.lookupMjdTt(tt(59170));
		ASSERT_TRUE(first.valid && last.valid);
		ASSERT_TRUE(first.poleOffsetsAvailable && last.poleOffsetsAvailable);
		EXPECT_FALSE(last.extrapolated);
		EXPECT_DOUBLE_EQ(first.xpRad, .1 * Constants::ARCSEC_RAD);
		EXPECT_DOUBLE_EQ(first.ypRad, -.2 * Constants::ARCSEC_RAD);
		EXPECT_DOUBLE_EQ(first.ttMinusUt1Seconds, 69.184 + .175);
		EXPECT_DOUBLE_EQ(std::get<Iau2000PoleOffsets>(first.poleOffsets).dxRad, Constants::ARCSEC_RAD / 1000.);
		EXPECT_DOUBLE_EQ(std::get<Iau2000PoleOffsets>(last.poleOffsets).dyRad, -4. * (Constants::ARCSEC_RAD / 1000.));
	}

	TEST_F(EarthOrientationTest, InterpolatesInTtWithoutLeapSecondRamp)
	{
		write({record(57753, .1, .2, -.4), record(57754, .3, .4, .6)});
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		const auto sample = table.lookupMjdTt((tt(57753) + tt(57754)) / 2.);
		ASSERT_TRUE(sample.valid);
		EXPECT_NEAR(sample.ttMinusUt1Seconds, 68.584, 3e-14);
		EXPECT_NEAR(sample.xpRad / Constants::ARCSEC_RAD, .2, 1e-11);
		EXPECT_NEAR(sample.ypRad / Constants::ARCSEC_RAD, .3, 1e-11);
	}

	TEST_F(EarthOrientationTest, OpticalClockIsIndependentOfTheSelectedPoleProfile)
	{
		write({record(59169, .1, .2, -.175), record(59170, .2, .3, -.177)});
		const double mjdTt = tt(59169.5);
		for (const auto profile :
		     {EarthOrientationProfile::IAU1980, EarthOrientationProfile::IAU2000, EarthOrientationProfile::FullC04}) {
			const EarthOrientation table(directory, profile);
			double mjdUt1 = 59169.5;
			ASSERT_TRUE(table.ttToUt1Mjd(mjdTt, mjdUt1));
			EXPECT_NEAR((mjdUt1 - 59169.5) * 86400.0, -.176, 1e-6);
			EXPECT_EQ(table.lookupMjdTt(mjdTt).valid, profile == EarthOrientationProfile::IAU2000);
		}
	}

	TEST_F(EarthOrientationTest, OpticalClockIteratesAtTheUt1Estimate)
	{
		write({record(59169, .1, .2, -.1), record(59170, .1, .2, -.9)});
		const EarthOrientation table(directory);
		const double mjdTt = tt(59169.5);
		double mjdUt1 = 59169.5;
		ASSERT_TRUE(table.ttToUt1Mjd(mjdTt, mjdUt1));
		const double offset = 69.284 + .8 * (mjdUt1 - 59169.);
		EXPECT_DOUBLE_EQ(mjdUt1, mjdTt - offset / 86400.0);
		const double singlePass = mjdTt - (69.284 + .8 * (mjdTt - 59169.)) / 86400.0;
		EXPECT_GT(std::abs(mjdUt1 - singlePass) * 86400.0, .0006);
	}

	TEST_F(EarthOrientationTest, OpticalClockPreservesContinuousOffsetAcrossLeapSecond)
	{
		write({record(57753, .1, .2, -.4), record(57754, .1, .2, .6), record(57755, .1, .2, .6)});
		const EarthOrientation table(directory);
		for (const double mjdUtc : {57753.5, 57754., 57754.5}) {
			const double mjdTt = tt(mjdUtc);
			double mjdUt1 = mjdUtc;
			ASSERT_TRUE(table.ttToUt1Mjd(mjdTt, mjdUt1));
			EXPECT_DOUBLE_EQ(mjdUt1, mjdTt - 68.584 / 86400.0);
		}
	}

	TEST_F(EarthOrientationTest, OpticalClockHandlesTtJustBeyondTheLastDailyEpoch)
	{
		write({record(59169, .1, .2, -.175), record(59170, .1, .2, -.175)});
		const EarthOrientation table(directory);
		double mjdUt1 = 59170.;
		ASSERT_TRUE(table.ttToUt1Mjd(tt(59170.), mjdUt1));
		EXPECT_DOUBLE_EQ(mjdUt1, tt(59170.) - (69.184 + .175) / 86400.0);
	}

	TEST_F(EarthOrientationTest, UnavailableOpticalClockLeavesCallerFallbackUnchanged)
	{
		double mjdUt1 = 59169.5;
		const EarthOrientation missing(directory);
		EXPECT_FALSE(missing.ttToUt1Mjd(tt(59169.5), mjdUt1));
		EXPECT_DOUBLE_EQ(mjdUt1, 59169.5);

		write({record(59169, .1, .2, -.1), record(59171, .1, .2, -.1)});
		const EarthOrientation table(directory);
		for (const double epoch : {tt(59168.), tt(59170.), tt(59172.), std::numeric_limits<double>::quiet_NaN(),
		                           std::numeric_limits<double>::infinity()}) {
			EXPECT_FALSE(table.ttToUt1Mjd(epoch, mjdUt1));
			EXPECT_DOUBLE_EQ(mjdUt1, 59169.5);
		}
	}

	TEST_F(EarthOrientationTest, PublishedRecordUsesBulletinANotBulletinB)
	{
		write(
		    {"201116 59169.00 I  0.132898 0.000018  0.288443 0.000042  I-0.1790143 0.0000066  0.2295 0.0042  I    "
		     "-0.026   "
		     " "
		     "0.024    -0.006    0.057  0.132881  0.288419 -0.1790137     0.031     0.026  "});
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		const auto sample = table.lookupMjdTt(tt(59169));
		ASSERT_TRUE(sample.valid && sample.poleOffsetsAvailable);
		EXPECT_DOUBLE_EQ(sample.xpRad, .132898 * Constants::ARCSEC_RAD);
		EXPECT_DOUBLE_EQ(sample.ypRad, .288443 * Constants::ARCSEC_RAD);
		EXPECT_DOUBLE_EQ(sample.ttMinusUt1Seconds, 69.184 + .1790143);
		const auto& offsets = std::get<Iau2000PoleOffsets>(sample.poleOffsets);
		EXPECT_DOUBLE_EQ(offsets.dxRad, -.026 * (Constants::ARCSEC_RAD / 1000.));
		EXPECT_DOUBLE_EQ(offsets.dyRad, -.006 * (Constants::ARCSEC_RAD / 1000.));
	}

	TEST_F(EarthOrientationTest, MissingOffsetsDoNotDiscardRotationData)
	{
		write({record(59169, .1, .2, -.1), record(59170, .2, .3, -.2, 0, 0, 'P', false)});
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		const auto sample = table.lookupMjdTt((tt(59169) + tt(59170)) / 2.);
		ASSERT_TRUE(sample.valid);
		EXPECT_FALSE(sample.poleOffsetsAvailable);
		EXPECT_DOUBLE_EQ(std::get<Iau2000PoleOffsets>(sample.poleOffsets).dxRad, 0.);
		EXPECT_DOUBLE_EQ(std::get<Iau2000PoleOffsets>(sample.poleOffsets).dyRad, 0.);
		EXPECT_GT(sample.xpRad, 0.);
		EXPECT_TRUE(table.lookupMjdTt(tt(59169)).poleOffsetsAvailable);
	}

	TEST_F(EarthOrientationTest, MissingMandatoryFieldsAndNonfiniteValuesAreSkipped)
	{
		auto missing = record(59169, .1, .2, -.1);
		missing.replace(58, 10, 10, ' ');
		auto invalid = record(59170, .1, .2, -.1);
		invalid.replace(18, 9, "      nan");
		write({"header", missing, invalid, record(59171, .1, .2, -.1)});
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		EXPECT_EQ(table.size(), 1u);
		EXPECT_TRUE(table.lookupMjdTt(tt(59171)).valid);
		EXPECT_FALSE(table.lookupMjdTt(tt(59170)).valid);
	}

	TEST_F(EarthOrientationTest, DuplicateOrUnsortedEpochsRejectTheTable)
	{
		for (double last : {59169., 59168.}) {
			write({record(59169, .1, .2, -.1), record(last, .1, .2, -.1)});
			const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
			EXPECT_FALSE(table.isAvailable());
			EXPECT_EQ(table.size(), 0u);
		}
	}

	TEST_F(EarthOrientationTest, DoesNotInterpolateAcrossMissingDays)
	{
		write({record(59169, .1, .2, -.1), record(59171, .1, .2, -.1)});
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		EXPECT_FALSE(table.lookupMjdTt(tt(59170)).valid);
		EXPECT_TRUE(table.lookupMjdTt(tt(59169)).valid);
		EXPECT_TRUE(table.lookupMjdTt(tt(59171)).valid);
	}

	TEST_F(EarthOrientationTest, OutsideFinalsCoverageUsesExplicitZeroPoleFallback)
	{
		write({record(59169, .1, .2, -.1), record(59170, .1, .2, -.1)});
		const EarthOrientation table(directory, EarthOrientationProfile::IAU2000);
		for (double epoch : {59168., 59171., 67898.}) {
			const auto sample = table.lookupMjdTt(tt(epoch));
			EXPECT_FALSE(sample.valid);
			EXPECT_TRUE(sample.extrapolated);
			EXPECT_FALSE(sample.poleOffsetsAvailable);
			EXPECT_DOUBLE_EQ(sample.xpRad, 0.);
			EXPECT_DOUBLE_EQ(sample.ypRad, 0.);
			EXPECT_DOUBLE_EQ(sample.ttMinusUt1Seconds, Convertors::utcToTdtSeconds(epoch + 2400000.5));
			EXPECT_DOUBLE_EQ(std::get<Iau2000PoleOffsets>(sample.poleOffsets).dxRad, 0.);
		}
	}

	TEST_F(EarthOrientationTest, UnavailableAndInvalidEpochsKeepTheSelectedOffsetType)
	{
		for (const auto profile : {EarthOrientationProfile::IAU1980, EarthOrientationProfile::IAU2000}) {
			const EarthOrientation table(directory, profile);
			ASSERT_FALSE(table.isAvailable());
			for (const double epoch :
			     {59169., std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
				const auto sample = table.lookupMjdTt(epoch);
				EXPECT_FALSE(sample.valid);
				EXPECT_FALSE(sample.extrapolated);
				EXPECT_TRUE(std::isfinite(sample.ttMinusUt1Seconds));
				EXPECT_EQ(std::holds_alternative<Iau1980PoleOffsets>(sample.poleOffsets),
				          profile == EarthOrientationProfile::IAU1980);
			}
		}
	}

	TEST_F(EarthOrientationTest, ProfilesCannotSilentlyFallBackToAnotherSource)
	{
		write({record(59169, .1, .2, -.1)});
		const EarthOrientation defaultProfile(directory);
		EXPECT_EQ(defaultProfile.profile(), EarthOrientationProfile::IAU1980);
		EXPECT_FALSE(defaultProfile.isAvailable());
		EXPECT_THROW(EarthOrientation(directory, static_cast<EarthOrientationProfile>(-1)), std::invalid_argument);
	}

	TEST_F(EarthOrientationTest, OffsetTypeSelectsTheRotationWithoutImplicitConversion)
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
		const auto legacy = EarthRotation::WithEop::poleEquatorialMjd(epoch, shifted);
		EXPECT_GT((legacy - p1).magnitude(), 1e-9);
	}

	TEST(EarthOrientationObserverTest, ProductionChainsRetainThe2020VvUt1Correction)
	{
		struct RestoreEnvironment {
			std::filesystem::path directory = getExecutableDir();
			IEphemerisSource* source = getDefaultEphemerisSource();

			~RestoreEnvironment()
			{
				setDefaultEphemerisSource(source);
				getExecutableDir() = directory;
			}
		} restore;

		getExecutableDir() = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data";
		DE440Source source;
		ASSERT_TRUE(source.isAvailable());
		setDefaultEphemerisSource(&source);

		constexpr double mjdUtc = 59139.214515046297;
		const double mjdTt = tt(mjdUtc);
		double mjdUt1 = mjdUtc;
		ASSERT_TRUE(earthOrientation().ttToUt1Mjd(mjdTt, mjdUt1));
		EXPECT_DOUBLE_EQ(mjdUt1, 59139.214513060193);
		obs::ParsedObservation parsed;
		parsed.jd = mjdUtc + 2400000.5;
		parsed.hasSplitUtc = true;
		parsed.utcMjdDay = 59139;
		parsed.utcSecondsOfDay = (mjdUtc - parsed.utcMjdDay) * 86400.0;
		parsed.obsCode = "F52";
		parsed.magnitudeText = "20.0";
		const ObsNode node(parsed, ObservatoryRecord{"F52", 203.74409, .936239, .351545});
		EXPECT_NEAR(node.LST, 5.3605187722485672, 2e-12);

		const Vec3D offsetEclAu(2.3981680378233689e-5, -2.3353858328807347e-5, 2.6409633567564912e-5);
		const Vec3D expectedBodyAu(-3.6537912322656392e-5, -1.6072567871481898e-5, 1.4988194117280977e-5);
		const auto bodyKm = hazard::icrfToBodyFixedMjd(hazard::Planet::Earth, offsetEclAu * 149597870.7,
		                                               time::toMjdTdb(mjdTt, time::Scale::TT));
		EXPECT_LT((bodyKm - expectedBodyAu * 149597870.7).magnitude(), 5e-5);
	}

}  // namespace
