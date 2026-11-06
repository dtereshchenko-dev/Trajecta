#pragma once

#include "orbital_core/FitConfig.h"
#include "orbital_core/FitRunControl.h"
#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/TimeScale.h"
#include "orbital_core/Vec3D.h"
#include "orbital_core/nongrav/NonGravModel.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace orbital_core {

	// ---------------------------------------------------------------------------
	// Ephemerides of an accepted solution.
	//
	// An ephemeris is a predicted observation at an arbitrary epoch, so it uses the
	// same model as the solution's observations:
	//   - the solution's force model, including its non-gravitational parameters
	//     whenever the solution has them, with close-approach step control; a
	//     solution fitted without N-body uses two-body motion, as in the fit;
	//   - the observer from the MPC station list, with UT1, at the reception epoch;
	//   - the two-pass light-time solution to that observer (reference aber2),
	//     including the Sun's barycentric motion during the delay.
	// Astrometric positions carry light time only (no stellar aberration) in the
	// J2000 equator or ecliptic. The model runs on TDB; UTC and TT are converted at
	// the boundary (TimeScale.h).
	// ---------------------------------------------------------------------------

	enum class EphemerisCenter : std::uint8_t { Heliocentric, Geocentric, Station };

	enum class EphemerisPlane : std::uint8_t { Equatorial, Ecliptic };
	enum class EphemerisMagnitudeType : std::uint8_t { Asteroid, CometTotal, CometNuclear };

	// The accepted solution: heliocentric ECLM J2000 state in AU and AU/day at its
	// fit epoch (MJD TDB), with the effective force model of that fit.
	struct EphemerisSolution {
		Vec3D position;
		Vec3D velocity;
		double epochMjdTdb{std::numeric_limits<double>::quiet_NaN()};
		NBodySettings forces;
		NonGravModel nonGrav;
		std::optional<double> absoluteMagnitude;
		double slopeG{0.15};
		EphemerisMagnitudeType magnitudeType{EphemerisMagnitudeType::Asteroid};
		double cometSlope{10.0};
		double cometPhaseCoefficient{0.0};
	};

	struct EphemerisRequest {
		EphemerisSolution solution;
		EphemerisCenter center{EphemerisCenter::Geocentric};
		std::string station;
		EphemerisPlane plane{EphemerisPlane::Equatorial};
		time::Scale scale{time::Scale::UTC};
		double startMjd{0.0};
		double endMjd{0.0};
		double stepDays{1.0};
		bool lightTime{true};
	};

	// Angles in radians, rates in radians per day, distances in AU. Position and
	// velocity are relative to the centre in the requested plane: astrometric
	// (emission position minus observer at reception) with light time, geometric
	// otherwise.
	struct EphemerisRow {
		double mjd{};  // grid epoch in the requested scale
		double mjdUtc{};
		double mjdTt{};
		double mjdTdb{};
		double longitude{};      // RA or ecliptic longitude, [0, 2 pi)
		double latitude{};       // Dec or ecliptic latitude
		double longitudeRate{};  // d(longitude)/dt * cos(latitude)
		double latitudeRate{};
		double distanceAu{};                                          // Delta, from the centre
		double heliocentricAu{};                                      // r
		double elongation{std::numeric_limits<double>::quiet_NaN()};  // heliocentric: NaN
		double phase{std::numeric_limits<double>::quiet_NaN()};       // heliocentric: NaN
		std::optional<double> visualMagnitude;
		Vec3D position;
		Vec3D velocity;
	};

	struct EphemerisResult {
		std::vector<EphemerisRow> rows;
		bool cancelled{false};
		bool lightTimeApplied{false};
		bool nBody{false};
		bool nonGravitational{false};
	};

	inline constexpr std::size_t kMaxEphemerisEpochs = 100000;

	// Counts start + i*step <= end using the same double arithmetic as generation.
	// Equal endpoints yield one epoch. For an interval, steps below its MJD spacing
	// are rejected instead of emitting duplicate dates. No fixed time tolerance is
	// added to the requested interval; exceeding the output cap returns cap + 1.
	[[nodiscard]] ORBITAL_CORE_API std::size_t ephemerisEpochCount(double startMjd, double endMjd, double stepDays);

	// Cooperative cancellation reaches the integrator corrector stages, grid and
	// observation loops. A cancelled result never contains partial ephemeris rows.
	[[nodiscard]] ORBITAL_CORE_API EphemerisResult computeEphemeris(const EphemerisRequest& request,
	                                                                const FitRunControl& control = {});

}  // namespace orbital_core
