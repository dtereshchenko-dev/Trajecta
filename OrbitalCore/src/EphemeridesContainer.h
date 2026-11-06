#pragma once
#include "orbital_core/Vec3D.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <memory>
#include <string>

class ORBITAL_CORE_API EphemeridesContainer {
   public:
	struct CelestianEphemerisData {
		double dobs{}, mv{}, vr{};
		double phaseAngle{};
		double jd{};
		Vec3D R;
		Vec3D V;
	};

	enum class Planet { Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Moon };
	enum class Error { NoError = 0, InvalidData, InsufficientData, IOError };

	explicit EphemeridesContainer(const Planet planetId);
	~EphemeridesContainer();

	CelestianEphemerisData get(const double& targetJD);
	CelestianEphemerisData getMjd(double targetMjd);

	Error error() const noexcept;
	std::string errorString() const noexcept;

   private:
	EphemeridesContainer() = delete;
	EphemeridesContainer(const EphemeridesContainer&) = delete;
	EphemeridesContainer& operator=(const EphemeridesContainer&) = delete;

	struct Impl;
	std::unique_ptr<Impl> d;
};
