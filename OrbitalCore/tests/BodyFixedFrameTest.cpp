#include <gtest/gtest.h>
#include "orbital_core/hazard/BodyFixedFrame.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/Types.h"
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <thread>

namespace {
	namespace hz = orbital_core::hazard;
	const auto data = std::filesystem::path(__FILE__).parent_path() / "data";

	Eigen::Matrix3d matrix(const hz::BodyFixedTransform& rotation)
	{
		Eigen::Matrix3d result;
		for (int i = 0; i < 3; ++i) {
			const auto& r = rotation.equatorialRows[i];
			result.row(i) << r.x, r.y, r.z;
		}
		return result;
	}

	class LibrationStub : public IEphemerisSource {
	   public:
		std::optional<LunarLibration> value;

		EphemerisRecord queryJD(EphemerisPlanet, double jd) override { return {jd, {}, {}}; }

		std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet, double, int) override { return {}; }

		bool isAvailable() const override { return true; }

		std::string name() const override { return "test"; }

		std::optional<LunarLibration> lunarLibrationMjdTdb(double) override { return value; }
	};

	TEST(BodyFixedFrame, MatchesIndependentSpiceMatricesForEveryNonEarthBody)
	{
		DE440Source ephemeris((data.parent_path().parent_path().parent_path() / "data/linux_p1550p2650.440").string());
		ASSERT_TRUE(ephemeris.isAvailable());
		std::ifstream input(data / "body_orientation_expected.json");
		ASSERT_TRUE(input.good());
		const auto fixture = nlohmann::json::parse(input);
		for (const auto& sample : fixture.at("samples")) {
			hz::Planet planet = hz::Planet::Count_;
			for (int i = 0; i < static_cast<int>(hz::Planet::Count_); ++i)
				if (hz::longLabel(static_cast<hz::Planet>(i)) == sample.at("body").get<std::string>())
					planet = static_cast<hz::Planet>(i);
			ASSERT_NE(planet, hz::Planet::Count_);
			const double mjd = sample.at("mjd_tdb");
			SCOPED_TRACE(sample.at("body").get<std::string>() + " " + std::to_string(mjd));
			const auto result = hz::bodyFixedTransformMjd(planet, mjd, &ephemeris);
			ASSERT_TRUE(result.valid);
			EXPECT_EQ(result.frame, sample.at("frame").get<std::string>());
			const Eigen::Matrix3d actual = matrix(result);
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					EXPECT_NEAR(actual(i, j), sample.at("matrix")[i][j].get<double>(), 2e-10);
			EXPECT_NEAR(actual.determinant(), 1.0, 2e-14);
			EXPECT_LT((actual * actual.transpose() - Eigen::Matrix3d::Identity()).norm(), 2e-14);
		}
	}

	TEST(BodyFixedFrame, MissingAndWrongVersionLibrationsAreNotAnApproximateSuccess)
	{
		LibrationStub source;
		EXPECT_FALSE(hz::bodyFixedTransformMjd(hz::Planet::Moon, 60000, &source).valid);
		auto fixed = hz::icrfToBodyFixedMjd(hz::Planet::Moon, {1737.4, 0, 0}, 60000, &source);
		EXPECT_FALSE(hz::latLonFromBodyFixed(hz::Planet::Moon, fixed).valid);
		source.value = LunarLibration{{.1, .4, .7}, 441};
		EXPECT_FALSE(hz::bodyFixedTransformMjd(hz::Planet::Moon, 60000, &source).valid);
		source.value->ephemerisVersion = 440;
		EXPECT_TRUE(hz::bodyFixedTransformMjd(hz::Planet::Moon, 60000, &source).valid);
		source.value->angles.x = std::numeric_limits<double>::quiet_NaN();
		EXPECT_FALSE(hz::bodyFixedTransformMjd(hz::Planet::Moon, 60000, &source).valid);
	}

	TEST(BodyFixedFrame, InvalidEpochAndPlanetDoNotReturnAValidFrame)
	{
		for (int i = 0; i < static_cast<int>(hz::Planet::Count_); ++i) {
			const auto p = static_cast<hz::Planet>(i);
			EXPECT_FALSE(hz::bodyFixedTransformMjd(p, std::numeric_limits<double>::infinity()).valid);
			EXPECT_FALSE(hz::latLonFromBodyFixed(p, {NAN, 1, 0}).valid);
		}
		EXPECT_FALSE(hz::bodyFixedTransformMjd(hz::Planet::Count_, 60000).valid);
		EXPECT_FALSE(hz::latLonFromBodyFixed({0, 0, 0}).valid);
	}

	TEST(BodyFixedFrame, EquatorialAndEclipticPathsAgreeIncludingRetrogradeBodies)
	{
		for (const auto p : {hz::Planet::Earth, hz::Planet::Mercury, hz::Planet::Venus, hz::Planet::Mars,
		                     hz::Planet::Jupiter, hz::Planet::Saturn, hz::Planet::Uranus, hz::Planet::Neptune}) {
			const Vec3D eq{1200, -5600, 4400};
			const auto rotation = hz::bodyFixedTransformMjd(p, 60000);
			ASSERT_TRUE(rotation.valid);
			const auto a = rotation.applyEquatorial(eq);
			const auto b = hz::icrfToBodyFixedMjd(p, Helpers::ecvToEcl(eq), 60000);
			EXPECT_NEAR((a - b).magnitude(), 0, 1e-9);
			EXPECT_NEAR(a.magnitude(), eq.magnitude(), 1e-9);
		}
	}

	TEST(BodyFixedFrame, LunarQueriesAreConcurrentAndDoNotChangePlanetStateCache)
	{
		DE440Source ephemeris((data.parent_path().parent_path().parent_path() / "data/linux_p1550p2650.440").string());
		ASSERT_TRUE(ephemeris.isAvailable());
		const auto before = ephemeris.queryMjd(EphemerisPlanet::Moon, 60000.25);
		const auto expected = hz::bodyFixedTransformMjd(hz::Planet::Moon, 60000.25, &ephemeris);
		ASSERT_TRUE(expected.valid);
		std::atomic<bool> correct = true;
		std::vector<std::thread> workers;
		for (int i = 0; i < 8; ++i)
			workers.emplace_back([&] {
				for (int j = 0; j < 50; ++j) {
					const auto r = hz::bodyFixedTransformMjd(hz::Planet::Moon, 60000.25, &ephemeris);
					const auto body = ephemeris.queryMjd(EphemerisPlanet::Moon, 60000.25);
					if (!r.valid || (matrix(r) - matrix(expected)).norm() != 0 ||
					    (body.R - before.R).magnitude() != 0 || (body.V - before.V).magnitude() != 0)
						correct = false;
				}
			});
		for (auto& worker : workers)
			worker.join();
		EXPECT_TRUE(correct);
		EXPECT_FALSE(ephemeris.lunarLibrationMjdTdb(-1e8).has_value());
	}
}  // namespace
