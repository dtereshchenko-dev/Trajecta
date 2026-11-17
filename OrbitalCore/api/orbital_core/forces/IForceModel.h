#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/DE440Source.h"
#include "orbital_core/ForceModelSettings.h"

#include <Eigen/Core>
#include <cstdint>
#include <string_view>
#include <concepts>

class DE440Source;
class SpkReader;
class NonGravModel;

namespace orbital_core::forces {
	struct ForceContext {
		Vec3D R;
		Vec3D V;
		// Keep the stage epoch and offset separate through DE/SPK lookup.
		double epochMjd = 0.0;
		double offsetDays = 0.0;
		// Rounded views for diagnostics and models without fine time dependence.
		double mjd;
		// Same instant as a Julian Date, for the few consumers still on that
		// convention. Derived, never a source: forming it costs ~40 us of
		// resolution that cannot be recovered by converting back.
		double jd;
		double centralNewtonianGm = 0.0;
		double sunGm = 0.0;
		double planetGms[8]{};
		double moonGm = 0.0;
		double earthGm = 0.0;
		double plutoGm = 0.0;
		EphemerisForceTerms ephemerisTerms;

		// Filled by PlanetaryPerturbations and reused by EIH/EarthJ2 so a
		// force evaluation performs only one batched DE440 lookup.
		Vec3D planetPositions[8]{};
		Vec3D planetVelocities[8]{};
		Vec3D moonPosition{};
		Vec3D moonVelocity{};
		bool havePlanets = false;
		bool havePlanetVelocities = false;

		// Per-call config carried from the propagator into the pipeline.
		DE440Source* ephemerisSource = nullptr;
		SpkReader* asteroidSource = nullptr;
		uint16_t bodyMask = 0;
		int targetNaifId = 0;
		bool useSb441 = false;
		ForceModelSettings forceModel;
		bool computePartials = true;
		const NonGravModel* nonGrav = nullptr;

		// Sets the epoch from its one authoritative form. Use this rather than
		// assigning the fields: jd must stay a derived view of mjd, never an
		// independent value someone set from a coarser source.
		void setEpochMjd(double epoch, double offset = 0.0) noexcept
		{
			epochMjd = epoch;
			offsetDays = offset;
			mjd = epoch + offset;
			jd = mjd + 2400000.5;
		}
	};

	template <typename T>
	concept ForceModel =
	    requires(const T& m, ForceContext& ctx, Vec3D& a, Eigen::Matrix3d& dadr, Eigen::Matrix3d& dadv) {
		    { m.evaluate(ctx, a, dadr, dadv) } noexcept;
		    { m.name() } noexcept -> std::convertible_to<std::string_view>;
	    };
}  // namespace orbital_core::forces
