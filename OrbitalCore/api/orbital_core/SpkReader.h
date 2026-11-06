#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <memory>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// SpkReader — minimal DAF/SPK Type 2 Chebyshev parser for JPL ephemeris files.
//
// Designed specifically to consume `sb441-n16.bsp` (the 16 perturbing asteroids
// used in DE441 integration): heliocentric ICRF/J2000 equatorial positions,
// position-only Chebyshev polynomial segments, center=Sun (NAIF id 10).
//
// No dependency on an external kernel toolkit; parses the binary file directly per the
// public NAIF spec (DAF Required Reading + SPK Required Reading).
//
// Public query API returns heliocentric ecliptic J2000 state in AU / AU·day⁻¹,
// matching the convention used by NBodyPropagator and DE440Source.  Times are
// Julian Dates in TDB (seconds past J2000 TDB internally).
//
// Thread safety: the segment index is built once at construction and is
// read-only afterwards.  Each caller thread opens its own FILE* lazily on
// first query (mirrors DE440Source's thread-local context pattern).
// ---------------------------------------------------------------------------
class ORBITAL_CORE_API SpkReader {
   public:
	SpkReader();
	explicit SpkReader(const std::string& bspPath);

	~SpkReader();

	SpkReader(const SpkReader&) = delete;
	SpkReader& operator=(const SpkReader&) = delete;

	bool isAvailable() const noexcept;
	[[nodiscard]] std::string sourcePath() const;

	// NAIF ids of all targets indexed in the file (e.g., 2000001..2000704 for
	// sb441-n16). Empty when !isAvailable().
	std::vector<int> targetIds() const;

	// Fetch heliocentric ecliptic J2000 state at TDB JD.
	// Returns true on success. On out-of-range or unknown naifId, returns
	// false and leaves R/V untouched.
	bool queryState(int naifId, double jdTdb, Vec3D& R, Vec3D& V) const;
	bool queryStateMjd(int naifId, double mjdTdb, Vec3D& R, Vec3D& V) const;
	bool queryStateOffsetMjd(int naifId, double epochMjd, double offsetDays, Vec3D& R, Vec3D& V) const;

	// Position-only fast path (velocity recurrence skipped). Returns zero
	// vector on failure — same convention as DE440Source::getPlanetPosition.
	Vec3D getPosition(int naifId, double jdTdb) const;
	Vec3D getPositionMjd(int naifId, double mjdTdb) const;
	Vec3D getPositionOffsetMjd(int naifId, double epochMjd, double offsetDays) const;

   private:
	struct Impl;
	std::unique_ptr<Impl> d;
};

// ---------------------------------------------------------------------------
// Global default asteroid source (set by application startup, consumed by
// NBodyPropagator). Mirrors getDefaultEphemerisSource() pattern from Types.h.
// ---------------------------------------------------------------------------
ORBITAL_CORE_API SpkReader* getDefaultAsteroidSource();
ORBITAL_CORE_API void setDefaultAsteroidSource(SpkReader* source);
