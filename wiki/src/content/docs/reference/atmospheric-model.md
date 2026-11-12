---
title: Atmospheric model
description: Atmospheric flight, fragmentation, environment and ground-impact distributions.
---

Atmospheric Monte Carlo uses StrewnField to propagate an ensemble of Earth
entry states through the atmosphere. The [user guide](../../guide/atmospheric-mc/)
explains the controls and report; this page describes the model behind them.

## Model and scope

StrewnField is a physics-based, semi-empirical model of atmospheric entry and
fragment fall. It combines equations of motion and ablation with parameterized
fragmentation, material properties and aerodynamic assumptions. It is designed
to explore possible outcomes under those assumptions, not to serve as a
high-fidelity simulation of a body's internal structure or its full breakup
history.

### Model assumptions

- **Fragmentation:** breakup uses a strength threshold and a prescribed mass
  distribution. Crack propagation, gradual erosion and aerodynamic interaction
  between fragments are not explicitly simulated.
- **Aerodynamics:** drag uses an effective area and a Mach-dependent
  coefficient; individual fragment orientation and rotation are not modelled.
- **Ablation:** mass loss follows an effective ablation law with a sampled
  coefficient, rather than a detailed thermal model of each fragment.
- **Atmosphere:** automatic GFS weather uses one snapshot and zero vertical
  wind. The sampled wind multiplier scales the entire wind vector uniformly. It
  does not change wind direction or independently perturb the wind profile at
  different altitudes.
- **Uncertainty:** the ensemble samples the specified input distributions. Its
  results depend on the selected entry states, material and environment.

The resulting field depends on the entry states, material and environment
chosen for the calculation. Increasing the number of realizations samples those
assumptions more thoroughly; it does not remove uncertainty in the assumptions
themselves. Retrospective comparisons with recovered meteorites are available
for [2023 CX1](../../verification/2023-cx1/) and
[2008 TC3](../../verification/2008-tc3/).

## State and equations of motion

Each tracked fragment has position **r**, velocity **v** and mass `m`. Position
and velocity are expressed in Earth-centred, Earth-fixed (ECEF) coordinates,
in metres and metres per second; time is in seconds. The position is converted
to WGS84 latitude, longitude and ellipsoidal height to sample the environment.

Writing `u = v - wind` for velocity relative to the air, the equations are:

```text
dr/dt = v
dv/dt = gravity + Coriolis + centrifugal - ρ Cd A |u| u / (2m)
dm/dt = -σa ρ A |u|³ / 2
```

Here `ρ` is atmospheric density, `A` is projected area, `Cd` is the drag
coefficient and `σa` is the ablation coefficient. The sampled wind and drag
multipliers scale the wind and `Cd` respectively. The rotating-frame terms are
`-2 Ω × v` and `-Ω × (Ω × r)`. Gravity includes a central Earth term and
Earth J2 in this ECEF model. This atmospheric gravity model is separate from
the N-body force model used to produce the orbital entry states.

For a spherical fragment of density `ρm`, `A = π (3m / (4πρm))^(2/3)`.
The core also has area-only brick and cone shapes. Drag coefficient is linearly
interpolated from a Mach-number table, with a constant value beyond Mach 4.
The speed of sound used for Mach number is `sqrt(401.874 T)` in m/s, with
temperature `T` in kelvin.

In entry mode, ablation and further breakup are evaluated while air-relative
speed is at least the configured luminous-speed threshold, 3 km/s by default.
Below that threshold mass remains constant and the model follows dark flight.
This speed threshold is a model switch, not a prediction of observed luminosity.

## Breakup and fragment masses

Breakup occurs when dynamic pressure reaches a fragment's strength. The UI
uses `q = ρ |u|²`; the core also supports the `q = ½ρ |u|²` convention. These
strength values differ by a factor of two for the same state.

For each breakup, a configured dust fraction can be removed from the parent.
Let `M` be the mass remaining after dust, `f` the largest-fragment fraction and
`β` the mass exponent. The largest fragment has mass `f M`. Subsequent
individual masses are generated in descending order, with `0 < f < 1` and
`0 < β < 1`:

```text
m(k) = M [f^(-β) + k β f^(1-β) / (1-β)]^(-1/β),  k = 1, 2, ...
```

Generation stops at the minimum fragment mass or when the remaining mass is
exhausted. The configured remainder policy determines whether an eligible
remaining piece is retained. Mass assigned to dust, fragments below the cutoff
and other residuals is recorded separately in the mass budget.

To avoid integrating every small fragment individually, the descending sequence
is kept individually up to a configured limit. The tail is grouped into
logarithmic mass bins; each bin is represented by a limited number of tracks
with multiplicities that preserve the count and mass of the group. First and
subsequent breakups can use different binning limits.

The strength of a child of mass `mc` from a parent of mass `mp` is
`min(strength_cap, parent_strength × (mp / mc)^alpha)`. Children receive
random lateral directions around the parent's velocity. Kick magnitudes depend
on air density, material density, relative air speed and the ratio of parent to
child diameter. The mass-weighted mean kick is removed so the represented
children retain the parent's mean velocity.

## Numerical integration and events

Flight uses an adaptive Dormand–Prince 5(4) Runge–Kutta method, not the orbital
RA15 integrator. The state has seven components. The embedded difference is
scaled using component-specific absolute tolerances and a relative tolerance;
the maximum normalized component controls step acceptance. The step is adjusted
with a safety factor and bounded by the configured limits.

Ground contact, breakup and the mass cutoff are checked on accepted steps.
When a sign change is detected, the solver uses cubic Hermite interpolation and
bisection to locate the earliest event within that step. Integration also has
flight-time and step-attempt limits, and a realization can time out or fail.
These outcomes are recorded separately from ground impacts.

## Atmosphere and terrain

The prepared environment is a regular latitude, longitude and height grid.
During flight, density, temperature and east/north/up wind components are
interpolated linearly in height and bilinearly in latitude and longitude. Terrain
height is bilinearly interpolated separately; ground contact uses clearance
above WGS84 terrain height. A trajectory outside the prepared horizontal grid
or above its atmospheric coverage produces an error.

For automatic weather, Trajecta selects the nearest three-hour NOAA GFS valid
snapshot, with no temporal interpolation; the valid time must be within 90
minutes of the midpoint of the entry-time range. At each GFS column, geopotential height is converted to
geometric height and adjusted using EGM2008 geoid undulation. Density at pressure
levels is computed as `p / (287.0531 T)` using the dry ideal-gas approximation.
Horizontal wind comes from GFS; vertical wind is set to zero.

When preparing the regular height grid, temperature and wind are interpolated
linearly between pressure levels, while density is interpolated logarithmically.
Below the lowest level the values are held constant. Above the highest level,
the profile transitions over 10 km to the standard atmosphere, with horizontal
wind fading to zero. The prepared grid is then interpolated during flight as
described above. Terrain combines Copernicus DEM and EGM2008 conversion to
ellipsoidal heights. The saved environment records source and timing metadata.

## Ground-impact density and contours

For a landed representative `i`, the count-map weight is
`wi = run_weight × multiplicity`; the mass-map weight is additionally
multiplied by landed fragment mass. The two maps are estimated and normalized
separately. Positions are projected to a local east/north plane and binned on
a regular grid before Gaussian smoothing.

Fixed KDE uses an isotropic Gaussian with the configured sigma. Adaptive KDE
starts with a Gaussian pilot whose covariance is the weighted sample covariance
scaled by `Neff^(-1/3)`, where `Neff = (Σ wi)² / Σ wi²`. At each sample the
pilot density `fi` determines a local scale:

```text
g = exp(Σ wi ln(fi) / Σ wi)
λi = clamp((fi / g)^(-alpha), lower_limit, upper_limit)
```

The implementation groups these scales into logarithmic classes and convolves
each weighted histogram with the corresponding Gaussian using FFT. If the
sample covariance is degenerate, it falls back to the configured fixed kernel.
The final grid is normalized to sum to one. HPD thresholds are found by sorting
cells by density and accumulating their normalized values to 50%, 80% or 99%.
The contours describe the chosen number- or mass-weighted distribution, not
the chance that at least one fragment lands in a region.

## Jensen–Shannon convergence

The stopping monitor is separate from the report KDE. It bins cumulative
impacts by fragment multiplicity and by landed mass, applies a fixed Gaussian
kernel, and normalizes each map. It compares successive cumulative maps with
base-2 Jensen–Shannon divergence:

```text
M = (P + Q) / 2
JS(P, Q) = 1/2 Σ Pi log2(Pi / Mi) + 1/2 Σ Qi log2(Qi / Mi)
```

Zero terms contribute zero. Both count and mass divergences must remain below
the configured threshold for the required consecutive checks. An incomplete
realization or an impact outside the monitor grid blocks early convergence.
The monitor's histograms use multiplicity and fragment mass; unlike the report
maps, this implementation does not multiply them by run weights.

## Code and comparisons

The implementation is in `StrewnField/src/TrajectoryIntegrator.cpp`,
`FragmentationModel.cpp`, `Environment.cpp`, `FieldDensity.cpp` and
`Convergence.cpp`; automatic weather preparation is in
`StrewnField/app/GfsSource.cpp`. See the event comparisons for
[2023 CX1](../../verification/2023-cx1/) and [2008 TC3](../../verification/2008-tc3/).

## Further development

Useful next steps are to compare additional events and study sensitivity to
material and weather assumptions, then assess improvements such as varying wind
direction and vertical structure, alternative aerodynamic representations and
more flexible fragmentation or erosion models. These are directions for future
work, not features included in the current model.
