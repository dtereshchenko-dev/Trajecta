#pragma once

#include "orbital_core/OrbitalCoreExport.h"
#include "orbital_core/StateVectors.h"

namespace Calculations {

	// Universal-variable propagation. State uses AU and V/k; dtDays is in days.
	[[nodiscard]] ORBITAL_CORE_API StateVectors propagateKeplerState(const StateVectors& initial, double dtDays);

	// initialRate differentiates [R, V/k] and dtRate differentiates dtDays.
	[[nodiscard]] ORBITAL_CORE_API StateVectors propagateKeplerStateRate(const StateVectors& initial, double dtDays,
	                                                                     const StateVectors& initialRate,
	                                                                     double dtRate);

	namespace detail {
		// tau and tauRate use canonical time (k times days).
		[[nodiscard]] bool keplerPropUV(const Vec3D& R0, const Vec3D& V0, double tau, Vec3D& Rout, Vec3D& Vout,
		                                const StateVectors* initialRate = nullptr, double tauRate = 0.0,
		                                Vec3D* positionRate = nullptr, Vec3D* velocityRate = nullptr) noexcept;
	}  // namespace detail

}  // namespace Calculations
