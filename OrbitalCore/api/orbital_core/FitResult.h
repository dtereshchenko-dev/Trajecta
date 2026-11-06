#pragma once

#include "orbital_core/Orbit.h"
#include "orbital_core/Residuals.h"
#include "orbital_core/nongrav/NonGravModel.h"
#include "orbital_core/uncertainty/DerivedSigmas.h"

#include <Eigen/Core>
#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

namespace orbital_core {

	// 6×6 cartesian (R [AU], V [AU/day]) state covariance, heliocentric
	// ecliptic.  Same `DontAlign` storage convention as ElementCov6x6
	// (no AVX-32 alignment requirements; safe on heap / DLL boundary).
	// Used by impact-risk / B-plane / LOV pipelines that need uncertainty
	// in physical state space rather than orbital-element space.
	using StateCov6x6 = Eigen::Matrix<double, 6, 6, Eigen::DontAlign>;
	using FitParameterCovariance = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::DontAlign>;

	// COT = (q [AU], e, i, Omega, omega, nu [rad]), heliocentric ecliptic.
	// Distinct from ElementCov6x6: its first coordinate is NOT semi-major axis.
	// Both matrices include all fitted non-grav contributions and cross terms.
	struct CometaryElementCovariance {
		using Matrix6 = Eigen::Matrix<double, 6, 6, Eigen::DontAlign>;
		Matrix6 covariance{Matrix6::Zero()};
		Matrix6 formalCovariance{Matrix6::Zero()};
		double epochJdTdb{std::numeric_limits<double>::quiet_NaN()};
	};

	// Post-fit uncertainty info — covariance matrices + 1-σ vectors.  Populated
	// by the orchestrator after LM convergence.  Empty when not computed (e.g.
	// failed fit, or path that skips covariance).  Future consumers:
	//   * MOID propagation — element covariance projected through close-
	//     encounter geometry to bound impact probability.
	//   * impact line-of-variations sampling.
	//   * Per-parameter σ reporting for downstream analysis.
	struct Uncertainty {
		// RMS-scaled orbital-element covariance (6×6: a, e, i, Ω, ω, ν) —
		// heliocentric ecliptic, same frame as Orbit's getters.  6th index is
		// true anomaly (ν), not mean anomaly — σ_M is propagated separately
		// into `derivedSigmas` below.  This drives reported orbit uncertainty;
		// Element-space visualization may use formal covariance separately.
		uncertainty::ElementCov6x6 elementCovariance{uncertainty::ElementCov6x6::Zero()};

		// Available on the near-parabolic chart at the actual display epoch.
		// At e == 1 the a-based elementCovariance is NaN/unavailable, while
		// this regular COT block and the fit-epoch Cartesian covariance remain valid.
		std::optional<CometaryElementCovariance> cometaryElementCovariance;

		// 6×6 cartesian state covariance in (R [AU], V [AU/day]) at the
		// **natural LM fit epoch** (== fitEpochOrbit->getEpoch()) — no STM
		// propagation. `cartesianStateCovariance` includes the reference-style
		// objective scale, while `cartesianStateCovarianceFormal` is the
		// unscaled (JᵀWJ)⁻¹ covariance. Impact MC and LOV use the same
		// final objective-scaled covariance; samplers never rescale it internally.
		//
		// Rationale: state-space covariance is only meaningful at the
		// epoch where the fit lives.  STM-
		// propagating it to the user's display epoch via linear cov
		// propagation amplifies the V·Δt cross-coupling by up to ×25-50
		// for 30-60 day extrapolations, breaking downstream impact-risk /
		// b-plane consumers. Downstream pipelines that need uncertainty at a later
		// epoch should either (a) draw samples here and integrate them
		// non-linearly through N-body to the target time, or
		// (b) chain-rule through to orbital
		// elements at the display epoch (already provided in
		// `elementCovariance`).
		StateCov6x6 cartesianStateCovariance{StateCov6x6::Zero()};
		StateCov6x6 cartesianStateCovarianceFormal{StateCov6x6::Zero()};
		double covarianceObjectiveScale{1.0};

		// Full marginal covariance at fitEpochOrbit->getEpoch() in physical
		// parameter layout:
		//   [R_x,R_y,R_z, V_x,V_y,V_z, non-grav...]
		// with R in AU, V in AU/day, non-grav in the model's physical units.
		// Non-grav order follows finalNonGrav.fitParameterMask():
		//   Yarkovsky: [A2]
		//   comet/direct models: [A1,A2], [A1,A2,A3], or [A1,A2,A3,DT]
		FitParameterCovariance fitParameterCovariance;
		FitParameterCovariance fitParameterCovarianceFormal;

		// Single source of truth for all 1-σ values:
		//   * direct: σ_a, σ_e, σ_i, σ_Ω, σ_ω, σ_ν  (= sqrt(diag(elementCovariance)))
		//   * derived: σ_M, σ_Tp, σ_q, σ_Q, σ_n, σ_P, σ_Py, σ_1/a
		// Populated from the available conventional or cometary chart.
		// Undefined parabolic a/M/n/Q/P sigmas are NaN, never zero placeholders.
		uncertainty::DerivedSigmas derivedSigmas;

		// Orbit condition/quality parameter computed from the same runoff
		// terms as `derivedSigmas` (NaN when undefined).
		double orbitQualityU{std::numeric_limits<double>::quiet_NaN()};

		// True only when the conventional (a,e,i,Omega,omega,nu) block is available.
		bool hasElementCovariance = false;

		// Cartesian availability is independent of angular-chart singularities.
		[[nodiscard]] bool hasStateCovariance() const noexcept
		{
			return cartesianStateCovariance.allFinite() && cartesianStateCovarianceFormal.allFinite() &&
			       !cartesianStateCovariance.isZero(0.0) && !cartesianStateCovarianceFormal.isZero(0.0);
		}

		[[nodiscard]] bool empty() const noexcept
		{
			return !hasElementCovariance && !cometaryElementCovariance.has_value() && !hasStateCovariance();
		}
	};

	enum class FitSelectionStatus {
		NotApplied,
		Stable,
		CorrectionFailed,
		PassLimit,
		InsufficientObservations,
		SafeguardRollback,
	};

	struct FitResult {
		// Orbit (R, V) at the user-facing display epoch.  R, V are
		// non-linearly N-body-propagated from the natural LM fit
		// epoch via the same force model the LM uses.  Orbital
		// elements / σ_a, σ_e ... are reported at this epoch (paired
		// with `uncertainty.elementCovariance`, which is chain-ruled
		// through to here from the fit-epoch state cov).
		Orbit_p orbit;

		// Orbit (R, V) at the **natural LM fit epoch** — the same
		// epoch where `uncertainty.cartesianStateCovariance` is
		// defined.  Impact-risk / obs-MC / b-plane pipelines start
		// propagation from this state. State and covariance live together
		// at the fit epoch; non-linear propagation handles transport to
		// any future close-approach time.
		// Identical to `orbit` when outputEpoch == fit epoch (e.g.
		// after a σ-clip refit that warm-starts at the propagated
		// orbit).
		Orbit_p fitEpochOrbit;

		// Residuals for the observations that actively participated in the
		// final fit. UI RMS/cost and covariance diagnostics use this set.
		Residuals residuals;

		// Optional residual field over the original immutable observation
		// list. This is diagnostic/export-only: rejected observations keep
		// their residuals without becoming part of the fitted objective.
		Residuals allResiduals;

		NonGravModel finalNonGrav;
		Uncertainty uncertainty;
		// Unweighted angular RMS in arcsec, independent of residual row output.
		double cost{0.0};
		// Dimensionless sqrt(r^T W r / m) for the final fitted observations,
		// using the objective's optical covariance and radar channel sigmas.
		// m counts scalar residuals (RA/Dec, delay, Doppler), not observations.
		// Independent of report rows and covariance scaling; NaN without
		// uncertainty weighting or measurements.
		double normalizedRms{std::numeric_limits<double>::quiet_NaN()};
		int iterations{0};
		// FitSession additionally requires automatic selection to stabilize;
		// a valid nominal state can be returned when this is false.
		bool converged{false};
		bool cancelled{false};
		FitSelectionStatus selectionStatus{FitSelectionStatus::NotApplied};
		int sigmaClipPassesRun{0};
		int sigmaClipRejected{0};
		std::vector<std::size_t> keptObsIndices;
	};
}  // namespace orbital_core
