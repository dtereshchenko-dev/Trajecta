#pragma once

#include "orbital_core/Vec3D.h"
#include <vector>
#include <string>
#include <optional>

struct EphemerisRecord {
	double jd;  // JD TDB
	Vec3D R;    // Position (AU, ecliptic heliocentric)
	Vec3D V;    // Velocity (AU/day)
	double dobs = 0;
	double phaseAngle = 0;
	double mv = 0;
	double vr = 0;
};

enum class EphemerisPlanet { Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Moon, Pluto };

// JPL principal-axis lunar Euler angles (phi, theta, psi), radians.
// Their reference frame is tied to the DE version; do not mix PA and ME.
struct LunarLibration {
	Vec3D angles;
	int ephemerisVersion{};
};

class IEphemerisSource {
   public:
	virtual ~IEphemerisSource() = default;

	// All epochs in this interface are TDB; velocities are derivatives per TDB day.
	virtual EphemerisRecord queryJD(EphemerisPlanet planet, double jd) = 0;

	virtual std::vector<EphemerisRecord> fetchWindow(EphemerisPlanet planet, double targetJD, int windowDays = 10) = 0;

	virtual bool isAvailable() const = 0;
	virtual std::string name() const = 0;

	// Barycentric Sun in J2000 ecliptic AU and AU/day, for radar relativity.
	// This is not the historical queryMjd(Sun) optical observer convention.
	virtual std::optional<EphemerisRecord> barycentricSunMjd(double) { return std::nullopt; }

	// Lunar orientation uses the same MJD TDB as position queries.
	virtual std::optional<LunarLibration> lunarLibrationMjdTdb(double) { return std::nullopt; }

	// Backends with a split-date kernel override this to retain sub-JD precision.
	virtual EphemerisRecord queryMjd(EphemerisPlanet planet, double mjd) { return queryJD(planet, 2400000.5 + mjd); }
};
