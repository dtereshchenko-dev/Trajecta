#pragma once

#include "orbital_core/Vec3D.h"
#include "orbital_core/OrbitalCoreExport.h"
#include <Eigen/Core>
#include <cstdint>

enum class NgModel : int {
	None = 0,
	MarsdenH2O = 1,      // Marsden, Sekanina & Yeomans 1973 — water-ice sublimation
	MarsdenCO = 2,       // CO sublimation (r0 ≈ 5 AU) — distant-active comets
	InverseSquare = 3,   // g(r) = (1 AU / r)², optionally with DT
	Yarkovsky = 4,       // g(r) = (1 AU / r)² — asteroid Yarkovsky (A2 only)
	InverseSquareDT = 5  // Legacy identifier; normalized to InverseSquare on input.
};

enum class NgParameter : std::uint8_t {
	A1 = 0,
	A2 = 1,
	A3 = 2,
	DT = 3,
};

enum class NgParameterMask : std::uint8_t {
	None = 0,
	A1 = 1u << 0,
	A2 = 1u << 1,
	A3 = 1u << 2,
	DT = 1u << 3,
};

[[nodiscard]] constexpr NgParameterMask operator|(NgParameterMask a, NgParameterMask b) noexcept
{
	return static_cast<NgParameterMask>(static_cast<std::uint8_t>(a) | static_cast<std::uint8_t>(b));
}

[[nodiscard]] constexpr NgParameterMask operator&(NgParameterMask a, NgParameterMask b) noexcept
{
	return static_cast<NgParameterMask>(static_cast<std::uint8_t>(a) & static_cast<std::uint8_t>(b));
}

[[nodiscard]] constexpr bool hasNgParameter(NgParameterMask mask, NgParameterMask bit) noexcept
{
	return static_cast<std::uint8_t>(mask & bit) != 0;
}

// Non-gravitational acceleration model (Marsden radial/transverse/normal terms
// plus optional activity-time-offset DT).  PHYSICS-ONLY: owns enabled flag,
// A1/A2/A3, DT, the g(r) law, the analytic acceleration partials, and the
// two-body Kepler back-propagation used when DT ≠ 0.
//
// `nongrav/NonGravProfile.h` adapts the model to the optimizer parameter
// vector. The force model remains independent of least-squares fitting.
class ORBITAL_CORE_API NonGravModel {
   public:
	// ---- Physics ----
	void setModel(NgModel m) noexcept;

	[[nodiscard]] NgModel model() const noexcept { return m_model; }

	void setEnabled(bool e) noexcept { m_enabled = e; }

	[[nodiscard]] bool enabled() const noexcept { return m_enabled && m_model != NgModel::None; }

	void setParams(double a1, double a2, double a3) noexcept
	{
		m_A1 = a1;
		m_A2 = a2;
		m_A3 = a3;
	}

	[[nodiscard]] double A1() const noexcept { return m_A1; }

	[[nodiscard]] double A2() const noexcept { return m_A2; }

	[[nodiscard]] double A3() const noexcept { return m_A3; }

	// Activity-time offset DT (days): a_ng(t) = A × g(r(t − DT)).
	// Supported by the Marsden and inverse-square laws; ignored by Yarkovsky/None.
	void setDT(double dtDays) noexcept { m_DT = dtDays; }

	[[nodiscard]] double DT() const noexcept { return m_DT; }

	[[nodiscard]] double gr(double rAU) const noexcept;  // heliocentric-distance factor
	[[nodiscard]] Vec3D acceleration(const Vec3D& r, const Vec3D& v, double jd) const noexcept;  // AU/day², inertial

	// State derivatives at fixed A/DT, in AU and AU/day. Both RTN directions
	// and g(r(t-DT)) depend on the state. These belong in the STM generator,
	// independently of the direct fit-parameter sensitivity forcing.
	void accelerationPartials(const Vec3D& r, const Vec3D& v, Eigen::Matrix3d& dadr,
	                          Eigen::Matrix3d& dadv) const noexcept;

	// Partial derivatives of a_NG wrt A1/A2/A3 at state (r, v).  These are
	// the three RTN direction vectors scaled by g(rForG) — independent of
	// current A values, so they are well-defined even when A=0 (which is
	// where the LM fit starts before the non-grav pass).  Used by the STM
	// sensitivity equations to avoid finite differences on A columns.
	//
	// Partials are zero when the model is disabled or NgModel::None.
	void sensitivityPartialsA(const Vec3D& r, const Vec3D& v, Vec3D& dA1, Vec3D& dA2, Vec3D& dA3) const noexcept;

	// Partial ∂a_NG/∂DT at state (r, v).  Zero when the model has no DT
	// parameter (Yarkovsky / None). Derived by
	// differentiating the Kepler back-propagation in rPastMagnitude —
	// no numerical differentiation of |r_past| needed.
	[[nodiscard]] Vec3D sensitivityPartialDT(const Vec3D& r, const Vec3D& v, double jd) const noexcept;

	// Number of parameters this model would contribute to an LM fit, IF it
	// were registered with a `NonGravProfile`.  Useful for sizing STM
	// sensitivity columns inside the propagator (which doesn't otherwise
	// touch fit-time concerns):
	//   - None:                    0
	//   - Yarkovsky:               1 (A2 only)
	//   - InverseSquare:           3 by default; DT may also be fitted
	//   - MarsdenH2O / MarsdenCO:  4 (A1, A2, A3, DT)
	// Returns 0 when the model is disabled.
	void setFitParameterMask(NgParameterMask mask) noexcept;
	[[nodiscard]] NgParameterMask fitParameterMask() const noexcept;
	[[nodiscard]] int fitParameterCount() const noexcept;
	[[nodiscard]] NgParameter fitParameterAt(int fitIndex) const noexcept;
	[[nodiscard]] int sensitivityColumnForFitParameter(int fitIndex) const noexcept;
	[[nodiscard]] double fitParameterValue(int fitIndex) const noexcept;
	void setFitParameterValue(int fitIndex, double value) noexcept;

	[[nodiscard]] static NgParameterMask defaultFitParameterMaskFor(NgModel model) noexcept;
	[[nodiscard]] static NgParameterMask sanitizeFitParameterMaskFor(NgModel model, NgParameterMask mask) noexcept;

	[[nodiscard]] static constexpr NgModel canonicalModel(NgModel model) noexcept
	{
		return model == NgModel::InverseSquareDT ? NgModel::InverseSquare : model;
	}

	[[nodiscard]] int numFitParams() const noexcept;

	// Max possible contribution (used for compile-time sizing where needed).
	static constexpr int maxFitParams = 4;

   private:
	// Return |r(t − dtDays)| via two-body Kepler back-propagation from (r, v).
	// Handles both elliptical (ε < 0) and hyperbolic (ε > 0) regimes.  Returns
	// |r| unchanged when dtDays == 0.
	[[nodiscard]] double rPastMagnitude(const Vec3D& r, const Vec3D& v, double dtDays) const noexcept;

	// Analytic gradient of the same radius, including implicit differentiation
	// of the elliptic/hyperbolic Kepler equation; columns are (r, v).
	[[nodiscard]] Eigen::Matrix<double, 1, 6> rPastGradient(const Vec3D& r, const Vec3D& v,
	                                                        double dtDays) const noexcept;
	[[nodiscard]] double grDerivative(double rAU) const noexcept;

	// Analytic ∂|r(t − dtDays)|/∂dtDays via differentiation of the same
	// Kepler back-propagation used in rPastMagnitude.  At dtDays=0 the
	// closed-form limit is −(r·v)/|r| (radial velocity, negated).  Returns
	// 0 for degenerate geometries (circular, near-parabolic).
	[[nodiscard]] double rPastDerivDT(const Vec3D& r, const Vec3D& v, double dtDays) const noexcept;

	// A supported DT also affects propagation when fixed at a nonzero value.
	[[nodiscard]] bool usesDT() const noexcept
	{
		return m_model == NgModel::MarsdenH2O || m_model == NgModel::MarsdenCO || m_model == NgModel::InverseSquare;
	}

	NgModel m_model{NgModel::MarsdenH2O};
	bool m_enabled{false};
	NgParameterMask m_fitMask{NgParameterMask::A1 | NgParameterMask::A2 | NgParameterMask::A3 | NgParameterMask::DT};
	double m_A1 = 0.0, m_A2 = 0.0, m_A3 = 0.0;  // AU/day²
	double m_DT = 0.0;                          // days
};
