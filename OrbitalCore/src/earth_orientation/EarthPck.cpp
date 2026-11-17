#include "EarthPck.h"
#include "EarthPckPoleCache.h"

#include "logger/Logger.h"
#include "orbital_core/Constants.h"
#include "orbital_core/RuntimeData.h"
#include "orbital_core/Types.h"

#include <SpiceUsr.h>
#include <array>
#include <cmath>
#include <filesystem>
#include <mutex>

namespace orbital_core {
	namespace {
		class EarthPck {
		   public:
			EarthPck()
			{
				char action[] = "RETURN";
				char reports[] = "NONE";
				erract_c("SET", 0, action);
				errprt_c("SET", 0, reports);
				for (std::size_t index = 0; index < m_names.size(); ++index) {
					const char* name = m_names[index];
					const auto path = earthPckFilePath(name);
					if (!std::filesystem::is_regular_file(path))
						continue;
					furnsh_c(path.string().c_str());
					if (failed_c()) {
						char message[1841]{};
						getmsg_c("LONG", sizeof(message), message);
						reset_c();
						LOG_WARN << "Earth PCK: cannot load " << path.string() << ": " << message;
					} else {
						++m_loaded;
						m_fileLoaded[index] = true;
						if (!m_poleCache.addFile(path))
							LOG_WARN << "Earth PCK pole cache unavailable; using CSPICE";
						LOG_INFO << "Earth PCK loaded: " << path.string();
					}
				}
				if (!m_loaded)
					LOG_WARN << "Earth PCK unavailable; Earth orientation requires C04 or historical UT1";
			}

			bool loaded(std::string_view name) const noexcept
			{
				for (std::size_t i = 0; i < m_names.size(); ++i)
					if (name == m_names[i])
						return m_fileLoaded[i];
				return false;
			}

			bool validate(const std::filesystem::path& path)
			{
				std::lock_guard lock(m_mutex);
				SPICEDOUBLE_CELL(coverage, 20000);
				pckcov_c(path.string().c_str(), 3000, &coverage);
				if (failed_c()) {
					reset_c();
					return false;
				}
				return card_c(&coverage) >= 2;
			}

			bool state(double mjdTdb, const Vec3D& terrestrial, Vec3D& position, Vec3D& velocity)
			{
				std::lock_guard lock(m_mutex);
				if (!m_loaded)
					return false;
				double matrix[6][6];
				const double et = (mjdTdb - 51544.5) * Constants::SECONDS_PER_DAY;
				sxform_c("ITRF93", "J2000", et, matrix);
				if (failed_c()) {
					reset_c();
					if (!m_warnedCoverage) {
						LOG_WARN << "Earth PCK does not cover an observer epoch; checking C04 or historical UT1";
						m_warnedCoverage = true;
					}
					return false;
				}
				const std::array<double, 3> input{terrestrial.x, terrestrial.y, terrestrial.z};
				std::array<double, 6> result{};
				for (int row = 0; row < 6; ++row)
					for (int col = 0; col < 3; ++col)
						result[row] += matrix[row][col] * input[col];
				position = {result[0], result[1], result[2]};
				velocity = {result[3] * Constants::SECONDS_PER_DAY, result[4] * Constants::SECONDS_PER_DAY,
				            result[5] * Constants::SECONDS_PER_DAY};
				return true;
			}

			bool pole(double mjdTdb, Vec3D& result)
			{
				if (m_poleCache.pole(mjdTdb, result))
					return true;
				Vec3D velocity;
				return state(mjdTdb, Vec3D(0.0, 0.0, 1.0), result, velocity);
			}

			bool toTerrestrial(double mjdTdb, const Vec3D& equatorial, Vec3D& terrestrial)
			{
				std::lock_guard lock(m_mutex);
				if (!m_loaded || !std::isfinite(mjdTdb))
					return false;
				double matrix[3][3];
				pxform_c("J2000", "ITRF93", (mjdTdb - 51544.5) * Constants::SECONDS_PER_DAY, matrix);
				if (failed_c()) {
					reset_c();
					return false;
				}
				const std::array<double, 3> input{equatorial.x, equatorial.y, equatorial.z};
				std::array<double, 3> result{};
				for (int row = 0; row < 3; ++row)
					for (int col = 0; col < 3; ++col)
						result[row] += matrix[row][col] * input[col];
				terrestrial = {result[0], result[1], result[2]};
				return true;
			}

			bool toTerrestrialState(double mjdTdb, const Vec3D& position, const Vec3D& velocity,
			                        Vec3D& terrestrialPosition, Vec3D& terrestrialVelocity)
			{
				std::lock_guard lock(m_mutex);
				if (!m_loaded || !std::isfinite(mjdTdb))
					return false;
				double matrix[6][6];
				sxform_c("J2000", "ITRF93", (mjdTdb - 51544.5) * Constants::SECONDS_PER_DAY, matrix);
				if (failed_c()) {
					reset_c();
					return false;
				}
				const std::array<double, 6> input{position.x, position.y, position.z,
				                                  velocity.x, velocity.y, velocity.z};
				std::array<double, 6> result{};
				for (int row = 0; row < 6; ++row)
					for (int col = 0; col < 6; ++col)
						result[row] += matrix[row][col] * input[col];
				terrestrialPosition = {result[0], result[1], result[2]};
				terrestrialVelocity = {result[3], result[4], result[5]};
				return true;
			}

		   private:
			EarthPckPoleCache m_poleCache;
			std::mutex m_mutex;
			const std::array<const char*, 3> m_names{"earth_predict.bpc", "earth_historic.bpc", "earth_latest.bpc"};
			std::array<bool, 3> m_fileLoaded{};
			int m_loaded = 0;
			bool m_warnedCoverage = false;
		};

		EarthPck& earthPck()
		{
			static EarthPck pck;
			return pck;
		}
	}  // namespace

	bool earthPckState(double mjdTdb, const Vec3D& terrestrial, Vec3D& position, Vec3D& velocity)
	{
		return earthPck().state(mjdTdb, terrestrial, position, velocity);
	}

	bool earthPckPole(double mjdTdb, Vec3D& pole) { return earthPck().pole(mjdTdb, pole); }

	bool earthPckLoaded(std::string_view filename) { return earthPck().loaded(filename); }

	bool validateEarthPck(const std::filesystem::path& path) { return earthPck().validate(path); }

	bool earthPckToTerrestrial(double mjdTdb, const Vec3D& equatorial, Vec3D& terrestrial)
	{
		return earthPck().toTerrestrial(mjdTdb, equatorial, terrestrial);
	}

	bool earthPckToTerrestrialState(double mjdTdb, const Vec3D& position, const Vec3D& velocity,
	                                Vec3D& terrestrialPosition, Vec3D& terrestrialVelocity)
	{
		return earthPck().toTerrestrialState(mjdTdb, position, velocity, terrestrialPosition, terrestrialVelocity);
	}
}  // namespace orbital_core
