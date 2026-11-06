#pragma once

#include "orbital_core/OrbitalCoreExport.h"

#include <cstdint>

namespace orbital_core::time {

	// ---------------------------------------------------------------------------
	// Native propagation time is MJD TDB, not JD.
	//
	// Not a style preference: a double holds JD 2 459 170 with an ulp of 4.66e-10 d
	// (40 us), while the same instant as MJD 59 169 has an ulp of 7.28e-12 d
	// (0.6 us) -- 64x finer. Rounding an epoch onto the JD grid displaces every
	// perturber by up to ~2 m along its orbit and biases the planetary
	// acceleration by ~1.8e-10 relative, systematically and with a fixed sign, so
	// it accumulates across a long propagation instead of averaging out.
	//
	// Therefore: everything from NBodyPropagator inwards -- integrators, dynamics,
	// force models, encounter control -- speaks MJD TDB. JD survives only at the
	// outer API, where callers still hand over Orbit::epoch, and is converted
	// exactly once, on entry. Do not reintroduce a JD round trip below that line;
	// jd <-> mjd is exact in binary, but the JD *value itself* is already
	// quantized, and no later conversion can recover what it dropped.
	//
	// Migration: as callers move to MJD, they should reach the MJD-native
	// NBodyPropagator entry points directly and the JD overloads can retire.
	// ---------------------------------------------------------------------------

	// JD of MJD zero. Exactly representable, so the conversions below are exact;
	// what is lossy is holding a modern epoch as JD in the first place.
	inline constexpr double kJdAtMjdZero = 2400000.5;

	[[nodiscard]] inline constexpr double mjdFromJd(double jd) noexcept { return jd - kJdAtMjdZero; }

	[[nodiscard]] inline constexpr double jdFromMjd(double mjd) noexcept { return mjd + kJdAtMjdZero; }

	// ---------------------------------------------------------------------------
	// Time scales at the user boundary.
	//
	// Orbit epochs, prepared observations, integrators, DE/SPK queries and hazard
	// events all use TDB. Convert source timestamps once on input, and convert
	// back only for presentation or a model that explicitly requires TT/UT1.
	//   UTC <-> TT : 32.184 s + TAI-UTC from the same drift / leap-second table the
	//                observation pipeline uses (Convertors::utcToTdtSeconds);
	//   UT1 <-> TT : measured DeltaT, for historical observation formats;
	//   TT  <-> TDB: Calculations::tdbMinusTtSeconds, the periodic series also
	//                used by every input path (~1.7 ms).
	// TT is local to Earth orientation, precession/nutation and clock models;
	// never pass it to a propagation or planetary ephemeris API.
	// ---------------------------------------------------------------------------

	enum class Scale : std::uint8_t {
		UTC,
		TT,
		TDB,
		UT1,
	};

	// MPC80 and RWO use UT1 before 1972, UTC afterwards. ADES records with
	// legacy MPC subFmt retain that convention; native ADES declares UTC.
	[[nodiscard]] inline constexpr Scale legacyObservationScale(double jd) noexcept
	{
		return jd < 2441317.5 ? Scale::UT1 : Scale::UTC;
	}

	// Offset added to a source MJD to obtain TT, in seconds. UTC includes the
	// pre-1972 frequency drift; UT1 uses measured Earth-orientation / DeltaT data.
	// Unavailable epochs (including UTC before 1960) return NaN, never a constant
	// extrapolation. Historical DeltaT is linearly interpolated in UT1.
	[[nodiscard]] ORBITAL_CORE_API double ttOffsetSeconds(double mjd, Scale scale) noexcept;

	// MJD in `scale` -> MJD TT, and back, for terrestrial clock/orientation models.
	// Exact inverses away from leap
	// seconds (to the ~1e-11 d resolution of an MJD double).
	[[nodiscard]] ORBITAL_CORE_API double toMjdTt(double mjd, Scale scale) noexcept;
	[[nodiscard]] ORBITAL_CORE_API double fromMjdTt(double mjdTt, Scale scale) noexcept;

	// Shared input/output boundary of the dynamical model. Already-TDB epochs
	// pass through unchanged, without a TT or JD round trip.
	[[nodiscard]] ORBITAL_CORE_API double toMjdTdb(double mjd, Scale scale) noexcept;
	[[nodiscard]] ORBITAL_CORE_API double fromMjdTdb(double mjdTdb, Scale scale) noexcept;

}  // namespace orbital_core::time
