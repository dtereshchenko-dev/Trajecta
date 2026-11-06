#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include <Eigen/Dense>

// ---------------------------------------------------------------------------
// VariationalEqs
//
// Header-only utilities for the variational equations of motion that ride
// alongside the state propagation in `IIntegrator::integrateStateWithStm`
// and `integrateStateWithStmAndSens`.
//
// Two related quantities are evolved next to (r, v):
//
//   * Φ ∈ ℝ⁶ˣ⁶ — the State Transition Matrix.  Initial value Φ(t₀, t₀) = I.
//                Satisfies dΦ/dt = F·Φ where F is the **state Jacobian**
//                of the equations of motion (see `buildStateJacobian`).
//                Maps initial-condition perturbations forward in time:
//                δx(t) = Φ(t, t₀) · δx(t₀).
//
//   * Ψ ∈ ℝ⁶ˣ⁴ — the **non-grav sensitivity matrix**.  Columns track how
//                the state responds to each fittable non-grav parameter
//                p_k ∈ {A₁, A₂, A₃, DT}.  Initial value Ψ(t₀, t₀) = 0.
//                Satisfies dΨ/dt = F·Ψ + G where G ∈ ℝ⁶ˣ⁴ stacks the
//                non-grav forcing partials (see `buildNongravSensForcing`)
//                onto the lower three rows (upper three are 0 because
//                d(dr/dt)/dp_k = 0 — non-grav touches acceleration only).
//
// `buildStateJacobian` and `buildNongravSensForcing` are kept here as
// `inline` helpers (no object file, no virtual dispatch) so the variational
// structure has a single source of truth and any future integrator can pull
// it in via `#include "propagator/VariationalEqs.h"` without copy-pasting.
// ---------------------------------------------------------------------------

namespace VariationalEqs {

	// 6×6 state Jacobian for r̈ = a(r, v, t):
	//
	//   F = [[  0  ,   I  ],
	//        [ ∂a/∂r, ∂a/∂v]]
	//
	// `dAdR` and `dAdV` are the 3×3 partials of acceleration w.r.t. position
	// and velocity respectively, supplied by the force model (gravity blocks
	// + Earth-J₂ + plugin forces all add into them in `accelerationAndPartials`).
	// For purely position-dependent forces `dAdV` is zero, but we let the caller
	// pass whatever Eigen-derived 3×3 it has — keeping the function ignorant of
	// which physical effects contribute.
	//
	// Returned by-value for use by an augmented integrator RHS.
	[[nodiscard]] inline Eigen::Matrix<double, 6, 6> buildStateJacobian(const Eigen::Matrix3d& dAdR,
	                                                                    const Eigen::Matrix3d& dAdV) noexcept
	{
		Eigen::Matrix<double, 6, 6> F;
		F.block<3, 3>(0, 0).setZero();
		F.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity();
		F.block<3, 3>(3, 0) = dAdR;
		F.block<3, 3>(3, 3) = dAdV;
		return F;
	}

	// 3×4 non-grav sensitivity forcing G:
	//
	//   G = [ ∂a_NG/∂A₁ | ∂a_NG/∂A₂ | ∂a_NG/∂A₃ | ∂a_NG/∂DT ]
	//
	// Three column vectors of analytic A-partials come from the RTN-direction
	// helper `NonGravModel::sensitivityPartialsA`; the DT column comes from
	// `sensitivityPartialDT`, which closes the Kepler-back-propagation
	// derivative analytically (no FD on |r_past|).
	//
	// Lives on the lower 3 rows of dΨ/dt = F·Ψ + Ĝ where Ĝ ∈ ℝ⁶ˣ⁴ has
	// `Ĝ.topRows<3>() = 0` (non-grav touches r̈, not ṙ) and
	// `Ĝ.bottomRows<3>() = G`.
	//
	// `t` is needed only for the DT path (Kepler back-propagation evaluates
	// at t − DT).  When the model has no DT (`Yarkovsky`,
	// `None`) the dDT column comes back as zero and the parameter is harmless.
	[[nodiscard]] inline Eigen::Matrix<double, 3, 4> buildNongravSensForcing(const NonGravModel& ng, const Vec3D& R,
	                                                                         const Vec3D& V, double t) noexcept
	{
		Vec3D dA1, dA2, dA3;
		ng.sensitivityPartialsA(R, V, dA1, dA2, dA3);
		const Vec3D dDT = ng.sensitivityPartialDT(R, V, t);

		Eigen::Matrix<double, 3, 4> G;
		G(0, 0) = dA1.x;
		G(1, 0) = dA1.y;
		G(2, 0) = dA1.z;
		G(0, 1) = dA2.x;
		G(1, 1) = dA2.y;
		G(2, 1) = dA2.z;
		G(0, 2) = dA3.x;
		G(1, 2) = dA3.y;
		G(2, 2) = dA3.z;
		G(0, 3) = dDT.x;
		G(1, 3) = dDT.y;
		G(2, 3) = dDT.z;
		return G;
	}
}  // namespace VariationalEqs
