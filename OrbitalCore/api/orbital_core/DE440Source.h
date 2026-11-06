#pragma once

#include "orbital_core/IEphemerisSource.h"
#include <string>
#include <vector>
#include "orbital_core/OrbitalCoreExport.h"

#include <cstdint>
#include <memory>
#include <string_view>

// Additional force inputs from the same exact-epoch planetary cache.
// Ecliptic J2000; Sun position/velocity/acceleration are barycentric, Pluto is heliocentric.
struct EphemerisForceTerms {
	Vec3D sunPosition;
	Vec3D sunVelocity;
	Vec3D sunAcceleration;
	Vec3D plutoPosition;
	Vec3D plutoVelocity;
};

// ---------------------------------------------------------------------------
// DE440 planetary ephemeris source.
//
// Wraps JPL's jpleph reader for DE440 / DE441 direct-access binary kernels.
// Thread-safe for concurrent queries after
// construction: every worker thread opens its own jpleph context on first
// query and keeps it alive in thread_local storage, so the OMP parallel
// loops in OrbitOptimizer::computeJacobian / NBodyPropagator::propagateAll*
// never serialise on a global jpl_pleph mutex.
// ---------------------------------------------------------------------------
class ORBITAL_CORE_API DE440Source : public IEphemerisSource {
   public:
	DE440Source();
	explicit DE440Source(std::string_view bspPath);

	~DE440Source() override;

	DE440Source(const DE440Source&) = delete;
	DE440Source& operator=(const DE440Source&) = delete;

	DE440Source(DE440Source&&) noexcept;
	DE440Source& operator=(DE440Source&&) noexcept;

	EphemerisRecord queryJD(EphemerisPlanet planet, double jd) override;
	EphemerisRecord queryMjd(EphemerisPlanet planet, double mjd) override;
	std::optional<EphemerisRecord> barycentricSunMjd(double mjd) override;
	std::optional<LunarLibration> lunarLibrationMjdTdb(double mjdTdb) override;
	std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet planet, double targetJD, int windowDays = 10) override;

	[[nodiscard]] bool isAvailable() const override;
	[[nodiscard]] std::string name() const override;
	[[nodiscard]] std::string sourcePath() const;

	// Heliocentric ecliptic position (AU) of a planet at jd.  Returns
	// zero for the Sun (heliocentric Sun is trivially at origin) and for
	// an unavailable source.
	[[nodiscard]] Vec3D getPlanetPosition(EphemerisPlanet planet, double jd) const;

	// MJD-native twin of getPlanetPosition. Routed through the MJD kernel
	// entry so the epoch is not quantized onto the JD grid; see TimeScale.h.
	[[nodiscard]] Vec3D getPlanetPositionMjd(EphemerisPlanet planet, double mjd) const;
	[[nodiscard]] Vec3D getPlanetPositionOffsetMjd(EphemerisPlanet planet, double epochMjd, double offsetDays) const;

	// Gravitational parameter from the active JPL kernel header, in AU^3/day^2.
	// Returns zero when the source is unavailable or the requested constant is
	// absent.  Earth and Moon are split from GMB using the kernel's EMRAT.
	[[nodiscard]] double gravitationalParameter(EphemerisPlanet planet) const noexcept;
	// MAxxxx in the same DE header; NAIF asteroid ID = 2000000 + MPC number.
	// Returns zero if unavailable, as for planetary constants above.
	[[nodiscard]] double asteroidGravitationalParameter(int naifId) const noexcept;

	// Batch heliocentric-ecliptic positions of all 8 perturbing planets
	// (Mercury..Neptune) at `jd`, into `out[0..7]`:
	//   0: Mercury  1: Venus  2: Earth (geocentric)  3: Mars
	//   4: Jupiter  5: Saturn 6: Uranus              7: Neptune
	void getAllPlanetPositions(double jd, Vec3D out[8], Vec3D* outMoon = nullptr) const;

	// Batch heliocentric-ecliptic states of all 8 perturbing planets and,
	// optionally, the Moon. Positions are in AU and velocities in AU/day.
	// This is the state-vector counterpart used by the reference EIH force.
	void getAllPlanetStates(double jd, Vec3D outPositions[8], Vec3D outVelocities[8], Vec3D* outMoonPosition = nullptr,
	                        Vec3D* outMoonVelocity = nullptr) const;

	// Retain sub-JD-ulp stage times during close-encounter integration.
	void getAllPlanetStatesMjd(double mjd, Vec3D outPositions[8], Vec3D outVelocities[8],
	                           Vec3D* outMoonPosition = nullptr, Vec3D* outMoonVelocity = nullptr,
	                           EphemerisForceTerms* forceTerms = nullptr) const;

	void getAllPlanetStatesOffsetMjd(double epochMjd, double offsetDays, Vec3D outPositions[8], Vec3D outVelocities[8],
	                                 Vec3D* outMoonPosition = nullptr, Vec3D* outMoonVelocity = nullptr,
	                                 EphemerisForceTerms* forceTerms = nullptr) const;

	// ---- Diagnostic: per-thread JD-cache statistics ------------------------
	struct CacheStats {
		std::uint64_t hits{0};       // JD found in cache (any probe)
		std::uint64_t misses{0};     // not found — fell through to jpl_state
		std::uint64_t evictions{0};  // probe depth exhausted on store, slot-0 overwritten
		std::uint64_t jplStateNanoseconds{0};

		[[nodiscard]] std::uint64_t totalLookups() const noexcept { return hits + misses; }

		[[nodiscard]] double hitRate() const noexcept
		{
			const auto t = totalLookups();
			return t == 0 ? 0.0 : static_cast<double>(hits) / static_cast<double>(t);
		}
	};

	[[nodiscard]] static CacheStats getCacheStats() noexcept;
	static void resetCacheStats() noexcept;

   private:
	struct Impl;
	std::unique_ptr<Impl> p_impl;
};
