#include "orbital_core/DE440Source.h"
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include "orbital_core/Constants.h"
#include "orbital_core/Types.h"
#include "io/DataSourceFactory.h"
#include "io/IDataSource.h"
#include "jpleph.h"
#include "jpl_int.h"
#include "ephemerides/SplitTime.h"

#include <array>
#include <atomic>
#include <chrono>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <memory>
#include <span>

namespace {
	std::atomic<std::uint64_t> s_nextSourceId{1};
	constexpr int kJplTargetSun = 11;
	constexpr int kJplCenterSun = 11;
	constexpr int kJplCenterEMB = 3;

	constexpr std::array<int, 11> kPlanetToTarget = {
	    kJplTargetSun, 1, 2, 3, 4, 5, 6, 7, 8, 10, 9,
	};

	constexpr int kJplPositionOnly = 0;
	constexpr int kJplWithVelocities = 1;
	constexpr std::size_t kStateVectorSize = 6;

	const double kCosObliquity = std::cos(Constants::EARTH_INCL);
	const double kSinObliquity = std::sin(Constants::EARTH_INCL);

	[[nodiscard]] inline Vec3D rotateEqToEcl(double x, double y, double z) noexcept
	{
		return Vec3D(x, y * kCosObliquity + z * kSinObliquity, z * kCosObliquity - y * kSinObliquity);
	}

	class ThreadEphemContext {
	   public:
		ThreadEphemContext() noexcept = default;

		ThreadEphemContext(const ThreadEphemContext&) = delete;
		ThreadEphemContext& operator=(const ThreadEphemContext&) = delete;

		~ThreadEphemContext() noexcept { closeContext(); }

		[[nodiscard]] void* acquire(std::uint64_t owner, const std::byte* viewData, std::size_t viewSize) noexcept
		{
			if (m_owner == owner) [[likely]]
				return m_ctx;

			closeContext();

			if (viewData != nullptr && viewSize > 0)
				m_ctx = jpl_init_ephemeris_from_view(viewData, viewSize, nullptr, nullptr);

			m_owner = owner;
			return m_ctx;
		}

	   private:
		void closeContext() noexcept
		{
			if (m_ctx != nullptr) {
				jpl_close_ephemeris(m_ctx);
				m_ctx = nullptr;
			}
		}

		void* m_ctx = nullptr;
		std::uint64_t m_owner = 0;
	};

	thread_local ThreadEphemContext t_ephem;

	struct EphemCacheEntry {
		Vec3D planetPositions[8];
		Vec3D planetVelocities[8];
		Vec3D moonPosition;
		Vec3D moonVelocity;
		EphemerisForceTerms forceTerms;
	};

	struct EphemCacheSlot {
		bool occupied{false};
		orbital_core::detail::SplitTime time{};
		EphemCacheEntry entry;
	};

	constexpr int kEphemCacheSize = 262144;
	constexpr int kEphemCacheMask = kEphemCacheSize - 1;
	constexpr int kEphemProbeDepth = 4;

	using EphemCacheArray = std::array<EphemCacheSlot, kEphemCacheSize>;

	std::atomic<std::uint64_t> s_jdCacheHits{0};
	std::atomic<std::uint64_t> s_jdCacheMisses{0};
	std::atomic<std::uint64_t> s_jdCacheEvictions{0};
	std::atomic<std::uint64_t> s_jplStateNanoseconds{0};

	thread_local std::unique_ptr<EphemCacheArray> t_ephemCachePtr;
	thread_local std::uint64_t t_ephemCacheOwner = 0;

	[[nodiscard]] inline EphemCacheArray& ephemCache() noexcept
	{
		if (!t_ephemCachePtr) [[unlikely]] {
			t_ephemCachePtr = std::make_unique<EphemCacheArray>();
		}

		return *t_ephemCachePtr;
	}

	[[nodiscard]] inline std::uint64_t mixDouble(double value) noexcept
	{
		if (value == 0.0)
			value = 0.0;  // +0 and -0 compare equal as cache keys.
		std::uint64_t bits;
		std::memcpy(&bits, &value, sizeof(bits));
		bits ^= bits >> 33;
		bits *= 0xff51afd7ed558ccdULL;
		bits ^= bits >> 33;
		return bits;
	}

	[[nodiscard]] inline std::size_t hashEpoch(orbital_core::detail::SplitTime time) noexcept
	{
		return static_cast<std::size_t>(mixDouble(time.high) ^ mixDouble(time.low)) & kEphemCacheMask;
	}

	inline const EphemCacheEntry* ephemCacheLookup(orbital_core::detail::SplitTime time) noexcept
	{
		auto& cache = ephemCache();
		const std::size_t h = hashEpoch(time);
		for (int probe = 0; probe < kEphemProbeDepth; ++probe) {
			const auto& s = cache[(h + probe) & kEphemCacheMask];
			if (!s.occupied) {
				s_jdCacheMisses.fetch_add(1, std::memory_order_relaxed);
				return nullptr;
			}
			if (s.time.high == time.high && s.time.low == time.low) {
				s_jdCacheHits.fetch_add(1, std::memory_order_relaxed);
				return &s.entry;
			}
		}
		s_jdCacheMisses.fetch_add(1, std::memory_order_relaxed);
		return nullptr;
	}

	inline void ephemCacheStore(orbital_core::detail::SplitTime time, const Vec3D planetPositions[8],
	                            const Vec3D planetVelocities[8], const Vec3D& moonPosition, const Vec3D& moonVelocity,
	                            const EphemerisForceTerms& forceTerms) noexcept
	{
		auto& cache = ephemCache();
		const std::size_t h = hashEpoch(time);
		for (int probe = 0; probe < kEphemProbeDepth; ++probe) {
			auto& s = cache[(h + probe) & kEphemCacheMask];
			if (!s.occupied || (s.time.high == time.high && s.time.low == time.low)) {
				s.occupied = true;
				s.time = time;

				for (int i = 0; i < 8; ++i) {
					s.entry.planetPositions[i] = planetPositions[i];
					s.entry.planetVelocities[i] = planetVelocities[i];
				}

				s.entry.moonPosition = moonPosition;
				s.entry.moonVelocity = moonVelocity;
				s.entry.forceTerms = forceTerms;
				return;
			}
		}

		s_jdCacheEvictions.fetch_add(1, std::memory_order_relaxed);

		auto& s = cache[h];
		s.occupied = true;
		s.time = time;
		for (int i = 0; i < 8; ++i) {
			s.entry.planetPositions[i] = planetPositions[i];
			s.entry.planetVelocities[i] = planetVelocities[i];
		}

		s.entry.moonPosition = moonPosition;
		s.entry.moonVelocity = moonVelocity;
		s.entry.forceTerms = forceTerms;
	}
}  // namespace

struct DE440Source::Impl {
	const std::uint64_t sourceId = s_nextSourceId.fetch_add(1, std::memory_order_relaxed);
	std::string filePath;
	std::unique_ptr<IDataSource> source;
	std::span<const std::byte> view;
	bool available{false};
	double emrat{0.0};
	double invOnePlusEmrat{0.0};
	std::array<double, 11> gmByPlanet{};
	std::array<double, 10000> gmByAsteroid{};

	bool openView(const std::string& pathStr)
	{
		auto src = DataSourceFactory::open(pathStr);
		if (!src || !src->isOpen())
			return false;

		const auto bytes = src->bytes();
		if (bytes.empty())
			return false;

		void* probe = jpl_init_ephemeris_from_view(bytes.data(), bytes.size(), nullptr, nullptr);
		if (probe == nullptr)
			return false;

		emrat = jpl_get_double(probe, JPL_EPHEM_EARTH_MOON_RATIO);
		invOnePlusEmrat = 1.0 / (1.0 + emrat);

		double gmEarthMoon = 0.0;
		const int nConstants = static_cast<int>(jpl_get_long(probe, JPL_EPHEM_N_CONSTANTS));
		for (int i = 0; i < nConstants; ++i) {
			char rawName[7]{};
			const double value = jpl_get_constant(i, probe, rawName);
			std::string_view name(rawName, 6);
			const auto end = name.find_first_of(" \0");
			if (end != std::string_view::npos)
				name = name.substr(0, end);
			if (name == "GMS")
				gmByPlanet[0] = value;
			else if (name == "GM1")
				gmByPlanet[1] = value;
			else if (name == "GM2")
				gmByPlanet[2] = value;
			else if (name == "GMB")
				gmEarthMoon = value;
			else if (name == "GM4")
				gmByPlanet[4] = value;
			else if (name == "GM5")
				gmByPlanet[5] = value;
			else if (name == "GM6")
				gmByPlanet[6] = value;
			else if (name == "GM7")
				gmByPlanet[7] = value;
			else if (name == "GM8")
				gmByPlanet[8] = value;
			else if (name == "GM9")
				gmByPlanet[10] = value;
			else if (name.size() == 6 && name.starts_with("MA")) {
				unsigned number = 0;
				const auto parsed = std::from_chars(name.data() + 2, name.data() + 6, number);
				if (parsed.ec == std::errc{} && parsed.ptr == name.data() + 6 && number < gmByAsteroid.size() &&
				    std::isfinite(value) && value > 0.0)
					gmByAsteroid[number] = value;
			}
		}
		if (gmEarthMoon > 0.0 && emrat > 0.0) {
			gmByPlanet[9] = gmEarthMoon * invOnePlusEmrat;
			gmByPlanet[3] = gmEarthMoon - gmByPlanet[9];
		}
		jpl_close_ephemeris(probe);

		filePath = pathStr;
		source = std::move(src);
		view = bytes;
		available = true;
		return true;
	}
};

double DE440Source::gravitationalParameter(EphemerisPlanet planet) const noexcept
{
	if (!p_impl || !p_impl->available)
		return 0.0;
	const auto index = static_cast<std::size_t>(planet);
	return index < p_impl->gmByPlanet.size() ? p_impl->gmByPlanet[index] : 0.0;
}

double DE440Source::asteroidGravitationalParameter(int naifId) const noexcept
{
	if (!p_impl || !p_impl->available || naifId < 2000000 || naifId >= 2010000)
		return 0.0;
	return p_impl->gmByAsteroid[static_cast<std::size_t>(naifId - 2000000)];
}

DE440Source::DE440Source() : p_impl(std::make_unique<Impl>())
{
	static constexpr std::array<const char*, 3> kCandidateFilenames = {"linux_p1550p2650.440", "linux_p1550p2650.441",
	                                                                   "de441.bin"};

	for (const char* fn : kCandidateFilenames) {
		const auto path = getExecutableDir() / fn;
		if (!std::filesystem::exists(path))
			continue;

		if (p_impl->openView(path.string()))
			return;
	}
}

DE440Source::DE440Source(std::string_view bspPath) : p_impl(std::make_unique<Impl>())
{
	const std::filesystem::path path(bspPath);
	if (!std::filesystem::exists(path))
		return;

	p_impl->openView(path.string());
}

DE440Source::~DE440Source() = default;
DE440Source::DE440Source(DE440Source&&) noexcept = default;
DE440Source& DE440Source::operator=(DE440Source&&) noexcept = default;

bool DE440Source::isAvailable() const { return p_impl->available; }

std::string DE440Source::name() const { return "DE440"; }

std::string DE440Source::sourcePath() const { return p_impl && p_impl->available ? p_impl->filePath : std::string{}; }

EphemerisRecord DE440Source::queryJD(EphemerisPlanet planet, double jd)
{
	EphemerisRecord record{};
	record.jd = jd;

	if (!p_impl->available)
		return record;

	void* ephem = t_ephem.acquire(p_impl->sourceId, p_impl->view.data(), p_impl->view.size());
	if (ephem == nullptr)
		return record;

	std::array<double, kStateVectorSize> rrd{};

	if (planet == EphemerisPlanet::Sun) {

		if (jpl_pleph(ephem, jd, kJplTargetSun, kJplCenterEMB, rrd.data(), kJplWithVelocities) != 0)
			return record;

		record.R = Vec3D(rrd[0], rrd[1], rrd[2]);
		record.V = Vec3D(rrd[3], rrd[4], rrd[5]);
		return record;
	}

	const int target = kPlanetToTarget[static_cast<std::size_t>(planet)];
	if (jpl_pleph(ephem, jd, target, kJplCenterSun, rrd.data(), kJplWithVelocities) != 0)
		return record;

	record.R = rotateEqToEcl(rrd[0], rrd[1], rrd[2]);
	record.V = rotateEqToEcl(rrd[3], rrd[4], rrd[5]);
	return record;
}

EphemerisRecord DE440Source::queryMjd(EphemerisPlanet planet, double mjd)
{
	EphemerisRecord record{};
	record.jd = 2400000.5 + mjd;
	Vec3D positions[8], velocities[8];
	EphemerisForceTerms forceTerms;
	getAllPlanetStatesMjd(mjd, positions, velocities, planet == EphemerisPlanet::Moon ? &record.R : nullptr,
	                      planet == EphemerisPlanet::Moon ? &record.V : nullptr,
	                      planet == EphemerisPlanet::Pluto ? &forceTerms : nullptr);
	if (planet == EphemerisPlanet::Sun) {
		const double c = std::cos(Constants::EARTH_INCL);
		const double s = std::sin(Constants::EARTH_INCL);
		const auto toSun = [c, s](const Vec3D& earth) {
			return Vec3D(-earth.x, -c * earth.y + s * earth.z, -s * earth.y - c * earth.z);
		};
		record.R = toSun(positions[2]);
		record.V = toSun(velocities[2]);
	} else if (planet == EphemerisPlanet::Pluto) {
		record.R = forceTerms.plutoPosition;
		record.V = forceTerms.plutoVelocity;
	} else if (planet >= EphemerisPlanet::Mercury && planet <= EphemerisPlanet::Neptune) {
		const auto index = static_cast<std::size_t>(planet) - 1;
		record.R = positions[index];
		record.V = velocities[index];
	}
	return record;
}

namespace {
	constexpr std::size_t kJplStateBodies = 13;
	constexpr std::size_t kJplStateCols = 6;
}  // namespace

void DE440Source::getAllPlanetPositions(double jd, Vec3D out[8], Vec3D* outMoon) const
{
	Vec3D unusedVelocities[8];
	Vec3D unusedMoonVelocity;
	getAllPlanetStates(jd, out, unusedVelocities, outMoon, &unusedMoonVelocity);
}

void DE440Source::getAllPlanetStates(double jd, Vec3D outPositions[8], Vec3D outVelocities[8], Vec3D* outMoonPosition,
                                     Vec3D* outMoonVelocity) const
{
	getAllPlanetStatesMjd(jd - 2400000.5, outPositions, outVelocities, outMoonPosition, outMoonVelocity);
}

std::optional<EphemerisRecord> DE440Source::barycentricSunMjd(double mjd)
{
	if (!p_impl->available)
		return std::nullopt;
	void* ephem = t_ephem.acquire(p_impl->sourceId, p_impl->view.data(), p_impl->view.size());
	if (!ephem)
		return std::nullopt;
	int list[14] = {};
	double pv[kJplStateBodies][kJplStateCols] = {}, nut[4] = {};
	if (jpl_state_mjd(ephem, mjd, list, pv, nut, 1) != 0)
		return std::nullopt;
	const auto* sun = static_cast<const jpl_eph_data*>(ephem)->pvsun;
	return EphemerisRecord{2400000.5 + mjd, rotateEqToEcl(sun[0], sun[1], sun[2]),
	                       rotateEqToEcl(sun[3], sun[4], sun[5])};
}

std::optional<LunarLibration> DE440Source::lunarLibrationMjdTdb(double mjdTdb)
{
	if (!p_impl || !p_impl->available || !std::isfinite(mjdTdb))
		return std::nullopt;
	void* ephem = t_ephem.acquire(p_impl->sourceId, p_impl->view.data(), p_impl->view.size());
	if (!ephem)
		return std::nullopt;
	const auto* header = static_cast<const jpl_eph_data*>(ephem);
	if (!header->ipt[12][0] || !header->ipt[12][1] || !header->ipt[12][2])
		return std::nullopt;
	int list[14]{};
	list[11] = 1;
	double pv[kJplStateBodies][kJplStateCols]{}, libration[6]{};
	if (jpl_state_mjd(ephem, mjdTdb, list, pv, libration, 1) != 0)
		return std::nullopt;
	if (!std::isfinite(libration[0]) || !std::isfinite(libration[1]) || !std::isfinite(libration[2]))
		return std::nullopt;
	return LunarLibration{{libration[0], libration[1], libration[2]}, static_cast<int>(header->ephemeris_version)};
}

void DE440Source::getAllPlanetStatesMjd(double mjd, Vec3D outPositions[8], Vec3D outVelocities[8],
                                        Vec3D* outMoonPosition, Vec3D* outMoonVelocity,
                                        EphemerisForceTerms* forceTerms) const
{
	getAllPlanetStatesOffsetMjd(mjd, 0.0, outPositions, outVelocities, outMoonPosition, outMoonVelocity, forceTerms);
}

void DE440Source::getAllPlanetStatesOffsetMjd(double epochMjd, double offsetDays, Vec3D outPositions[8],
                                              Vec3D outVelocities[8], Vec3D* outMoonPosition, Vec3D* outMoonVelocity,
                                              EphemerisForceTerms* forceTerms) const
{
	const auto time = orbital_core::detail::SplitTime::sum(epochMjd, offsetDays);
	if (forceTerms)
		*forceTerms = {};
	if (!p_impl->available) {
		for (int i = 0; i < 8; ++i) {
			outPositions[i] = Vec3D(0.0, 0.0, 0.0);
			outVelocities[i] = Vec3D(0.0, 0.0, 0.0);
		}

		if (outMoonPosition)
			*outMoonPosition = Vec3D(0.0, 0.0, 0.0);
		if (outMoonVelocity)
			*outMoonVelocity = Vec3D(0.0, 0.0, 0.0);

		return;
	}

	if (t_ephemCacheOwner != p_impl->sourceId) {
		t_ephemCachePtr.reset();
		t_ephemCacheOwner = p_impl->sourceId;
	}
	if (const auto* hit = ephemCacheLookup(time)) {
		for (int i = 0; i < 8; ++i) {
			outPositions[i] = hit->planetPositions[i];
			outVelocities[i] = hit->planetVelocities[i];
		}

		if (outMoonPosition)
			*outMoonPosition = hit->moonPosition;
		if (outMoonVelocity)
			*outMoonVelocity = hit->moonVelocity;
		if (forceTerms)
			*forceTerms = hit->forceTerms;

		return;
	}

	void* ephem = t_ephem.acquire(p_impl->sourceId, p_impl->view.data(), p_impl->view.size());
	if (ephem == nullptr) {
		for (int i = 0; i < 8; ++i) {
			outPositions[i] = Vec3D(0.0, 0.0, 0.0);
			outVelocities[i] = Vec3D(0.0, 0.0, 0.0);
		}

		if (outMoonPosition)
			*outMoonPosition = Vec3D(0.0, 0.0, 0.0);
		if (outMoonVelocity)
			*outMoonVelocity = Vec3D(0.0, 0.0, 0.0);

		return;
	}

	int list[14] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0};
	double pv[kJplStateBodies][kJplStateCols] = {};
	double nut[4] = {};

	const auto jplStarted = std::chrono::steady_clock::now();
	const int rval = jpl_state_offset_mjd(ephem, time.high, time.low, list, pv, nut, 1);
	const auto jplElapsed =
	    std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - jplStarted).count();
	s_jplStateNanoseconds.fetch_add(static_cast<std::uint64_t>(jplElapsed), std::memory_order_relaxed);
	if (rval != 0) {
		for (int i = 0; i < 8; ++i) {
			outPositions[i] = Vec3D(0.0, 0.0, 0.0);
			outVelocities[i] = Vec3D(0.0, 0.0, 0.0);
		}

		if (outMoonPosition)
			*outMoonPosition = Vec3D(0.0, 0.0, 0.0);
		if (outMoonVelocity)
			*outMoonVelocity = Vec3D(0.0, 0.0, 0.0);

		return;
	}

	{
		const double inv = p_impl->invOnePlusEmrat;
		pv[2][0] -= pv[9][0] * inv;
		pv[2][1] -= pv[9][1] * inv;
		pv[2][2] -= pv[9][2] * inv;
		pv[2][3] -= pv[9][3] * inv;
		pv[2][4] -= pv[9][4] * inv;
		pv[2][5] -= pv[9][5] * inv;
	}

	const auto* eph = static_cast<const jpl_eph_data*>(ephem);
	const double* pvsun = eph->pvsun;
	EphemerisForceTerms terms;
	terms.sunPosition = rotateEqToEcl(pvsun[0], pvsun[1], pvsun[2]);
	terms.sunVelocity = rotateEqToEcl(pvsun[3], pvsun[4], pvsun[5]);
	terms.sunAcceleration = rotateEqToEcl(pvsun[6], pvsun[7], pvsun[8]);
	terms.plutoPosition = rotateEqToEcl(pv[8][0] - pvsun[0], pv[8][1] - pvsun[1], pv[8][2] - pvsun[2]);
	terms.plutoVelocity = rotateEqToEcl(pv[8][3] - pvsun[3], pv[8][4] - pvsun[4], pv[8][5] - pvsun[5]);
	if (forceTerms)
		*forceTerms = terms;
	for (int i = 0; i < 8; ++i) {
		pv[i][0] -= pvsun[0];
		pv[i][1] -= pvsun[1];
		pv[i][2] -= pvsun[2];
		pv[i][3] -= pvsun[3];
		pv[i][4] -= pvsun[4];
		pv[i][5] -= pvsun[5];
	}

	Vec3D moonPositionEcl;
	Vec3D moonVelocityEcl;
	{
		const double moonEqX = pv[2][0] + pv[9][0];
		const double moonEqY = pv[2][1] + pv[9][1];
		const double moonEqZ = pv[2][2] + pv[9][2];
		moonPositionEcl = rotateEqToEcl(moonEqX, moonEqY, moonEqZ);

		const double moonEqVx = pv[2][3] + pv[9][3];
		const double moonEqVy = pv[2][4] + pv[9][4];
		const double moonEqVz = pv[2][5] + pv[9][5];
		moonVelocityEcl = rotateEqToEcl(moonEqVx, moonEqVy, moonEqVz);
	}

	if (outMoonPosition)
		*outMoonPosition = moonPositionEcl;
	if (outMoonVelocity)
		*outMoonVelocity = moonVelocityEcl;

	for (int i = 0; i < 8; ++i) {
		outPositions[i] = rotateEqToEcl(pv[i][0], pv[i][1], pv[i][2]);
		outVelocities[i] = rotateEqToEcl(pv[i][3], pv[i][4], pv[i][5]);
	}

	ephemCacheStore(time, outPositions, outVelocities, moonPositionEcl, moonVelocityEcl, terms);
}

Vec3D DE440Source::getPlanetPosition(EphemerisPlanet planet, double jd) const
{
	if (!p_impl->available || planet == EphemerisPlanet::Sun)
		return {0.0, 0.0, 0.0};

	void* ephem = t_ephem.acquire(p_impl->sourceId, p_impl->view.data(), p_impl->view.size());
	if (ephem == nullptr)
		return {0.0, 0.0, 0.0};

	std::array<double, kStateVectorSize> rrd{};
	const int target = kPlanetToTarget[static_cast<std::size_t>(planet)];

	if (jpl_pleph(ephem, jd, target, kJplCenterSun, rrd.data(), kJplPositionOnly) != 0)
		return {0.0, 0.0, 0.0};

	return rotateEqToEcl(rrd[0], rrd[1], rrd[2]);
}

Vec3D DE440Source::getPlanetPositionMjd(EphemerisPlanet planet, double mjd) const
{
	return getPlanetPositionOffsetMjd(planet, mjd, 0.0);
}

Vec3D DE440Source::getPlanetPositionOffsetMjd(EphemerisPlanet planet, double epochMjd, double offsetDays) const
{
	if (!p_impl->available || planet == EphemerisPlanet::Sun)
		return {0.0, 0.0, 0.0};

	Vec3D positions[8];
	Vec3D velocities[8];
	Vec3D moonPosition;
	EphemerisForceTerms forceTerms;
	getAllPlanetStatesOffsetMjd(epochMjd, offsetDays, positions, velocities, &moonPosition, nullptr,
	                            planet == EphemerisPlanet::Pluto ? &forceTerms : nullptr);

	if (planet == EphemerisPlanet::Moon)
		return moonPosition;
	if (planet == EphemerisPlanet::Pluto)
		return forceTerms.plutoPosition;

	const int index = static_cast<int>(planet) - static_cast<int>(EphemerisPlanet::Mercury);
	if (index < 0 || index >= 8)
		return {0.0, 0.0, 0.0};
	return positions[index];
}

std::vector<EphemerisRecord> DE440Source::fetchWindow(EphemerisPlanet planet, double targetJD, int windowDays)
{
	std::vector<EphemerisRecord> records;

	if (!p_impl->available)
		return records;

	records.reserve(static_cast<std::size_t>(2 * windowDays + 1));

	for (int day = -windowDays; day <= windowDays; ++day)
		records.push_back(queryJD(planet, targetJD + day));

	return records;
}

DE440Source::CacheStats DE440Source::getCacheStats() noexcept
{
	return CacheStats{
	    s_jdCacheHits.load(std::memory_order_relaxed),
	    s_jdCacheMisses.load(std::memory_order_relaxed),
	    s_jdCacheEvictions.load(std::memory_order_relaxed),
	    s_jplStateNanoseconds.load(std::memory_order_relaxed),
	};
}

void DE440Source::resetCacheStats() noexcept
{
	s_jdCacheHits.store(0, std::memory_order_relaxed);
	s_jdCacheMisses.store(0, std::memory_order_relaxed);
	s_jdCacheEvictions.store(0, std::memory_order_relaxed);
	s_jplStateNanoseconds.store(0, std::memory_order_relaxed);
}
