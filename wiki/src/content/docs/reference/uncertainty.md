---
title: Uncertainty and interpretation
description: Fit covariance, uncertainty propagation and interpretation of results.
---

An orbit fit produces a nominal solution and a covariance matrix. The matrix
describes how uncertainties in the fitted parameters vary together, based on
the measurements, weights and model used in the fit. Trajecta uses it to
estimate uncertainties in orbital elements, ephemerides and close approaches,
and to sample alternative orbits.

## Fit covariance

At the fit epoch, the parameter vector contains the heliocentric Cartesian
position and velocity. If non-gravitational coefficients are estimated, the
matrix also includes those fitted parameters and their correlations with the
orbital state. Fixed coefficients are not extra fitted dimensions.

The fit forms a Jacobian of the residuals with respect to the parameters,
applies the observation weights and excludes residuals outside the accepted
fit mask. The formal covariance is obtained from the weighted Jacobian using
singular-value decomposition. If the Jacobian is rank-deficient or invalid,
the fit does not publish a usable covariance.

Trajecta also stores an objective-scaled covariance. For `n` active scalar
residuals, `p` fitted parameters and weighted residual sum of squares `χ²`,
the variance factor is:

```text
factor = max(n, χ²) / (n − p)
scaled covariance = factor × formal covariance
```

When the inputs are invalid or `n ≤ p`, the implementation uses a factor of
1. The factor is applied to the covariance, not directly to its standard
deviations. Both formal and scaled matrices are available in `covariance.txt`;
the scaled covariance is used for the solution diagnostics and calculations
described below. **COV scale** in the interface displays this factor.

## Reading the matrix

For parameter `i`, the standard deviation shown as `σ` is the square root of
the diagonal covariance term `Cᵢᵢ`. A covariance between two parameters has
the product of their units. The corresponding dimensionless correlation is:

```text
ρᵢⱼ = Cᵢⱼ / sqrt(Cᵢᵢ Cⱼⱼ)
```

Correlations near +1 or −1 mean that the parameters vary strongly together
or in opposite directions in the local fit. They are useful for identifying
poorly separated parameters, but do not on their own describe how uncertainty
will evolve over time.

On **Solution → Covariance**, switch between **COV** and **Correlation** and
select a cell to see its value and the individual sigmas. The Cartesian view
uses AU and AU/day. When available, the equinoctial view shows `a`, `h`, `k`,
`p`, `q` and `λ`; the matrix displays `λ` in degrees. Fitted
non-gravitational parameters appear alongside the orbital parameters.

The matrix view is associated with the **fit epoch** shown above it. The
orbital elements displayed elsewhere can instead use a different output epoch.
Compare covariances only after checking their epoch, coordinate basis and units.

**Weak direction** is the eigenvector associated with the largest eigenvalue
of the covariance in the selected basis; **Weak σ** is its square root. Its
components depend on the basis and units, so they should not be compared as
unit-independent measures of the contribution of each parameter. In the
equinoctial weak-direction calculation, `λ` uses radians even though the
displayed covariance matrix uses degrees.

## Propagation to another epoch

For a small change in the fitted parameters `δx₀`, the variational equations
provide a local linear approximation to the change at a later epoch:

```text
δx(t) ≈ A(t) δx₀
C(t)  ≈ A(t) C₀ A(t)ᵀ
```

For a six-parameter fit, `A` contains the state transition matrix. When the fit
includes non-gravitational parameters, it also contains the corresponding
sensitivities, so their cross-covariances contribute to the propagated state.
Transformations from Cartesian state to orbital elements use the appropriate
local Jacobian as well.

This describes the local behaviour around the nominal trajectory. A strongly
nonlinear encounter can distort an uncertainty region beyond what one linear
covariance describes; sampled trajectories can help examine that behaviour.

## Ephemeris uncertainties

When uncertainty output is requested and a suitable covariance is available,
Trajecta propagates it to each ephemeris epoch and projects the position
uncertainty onto the sky. The resulting two-dimensional covariance provides
the **3σ error ellipse**: its major and minor semiaxes, position angle and
area. Component uncertainties and distance uncertainty are also available in
the output, as is uncertainty in the calculated observer range rate.

These are local, covariance-based estimates for the selected observer and
ephemeris settings. The 3σ ellipse is not the same statement as a 3σ interval
for one coordinate. See [Ephemerides and close approaches](../../ephemerides-and-close-approaches/)
for the output and light-time conventions.

## Close-approach uncertainties

The close-approach catalogue follows the nominal trajectory and reports its
time and distance from the centre of each encountered body. Where a usable
covariance is available, Trajecta propagates it to the encounter and estimates
a local **3σ time uncertainty**.

The minimum and maximum distances use a linearized 3σ B-plane ellipse and
nominal-encounter-energy focusing. They are not simply the nominal distance
plus or minus three times a distance sigma. An invalid covariance or degenerate
encounter can leave uncertainty values unavailable. The catalogue does not
calculate impact probability.

## Sampling and impact probability

**Covariance MC** samples a multivariate normal distribution centred on the
fitted state at the fit epoch. It uses the full scaled parameter covariance,
including fitted non-gravitational terms and their correlations. Each sample is
then propagated to the selected search interval. The covariance must have the
expected dimensions and pass the positive-definiteness check.

**Observation MC** instead perturbs the accepted measurements and refits an
orbit for each trial. **MOV** explores a grid of fitted nodes with weights.
These methods describe uncertainty differently, so their estimates need not
match. Their settings and interpretation are described in
[Orbital Monte Carlo](../../guide/monte-carlo/) and
[Manifold of Variations](../../guide/mov/).

For an orbital MC event, Trajecta divides the number of completed trials that
crossed the boundary in that event by the number of completed trials. The
displayed 95% Wilson interval describes finite-sample uncertainty under the
adopted model. Fit failures, propagation failures, cancelled and unstarted
trials are not counted as misses; review their counts when interpreting an
incomplete run. No crossings in a finite sample does not establish zero impact
probability.

## Boundary and atmosphere

A boundary crossing refers to the selected target, altitude and search
interval. Atmospheric MC is a separate calculation, conditional on its supplied
Earth entry states. Its distribution of outcomes is not the probability of
the object entering Earth's atmosphere in the first place.

All of these estimates depend on the adopted observations, weights and
physical model. Formal covariance and finite-sample intervals do not by
themselves include every systematic error or unmodelled force.
