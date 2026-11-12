---
title: Photometry
description: Magnitude models, fitted absolute magnitudes and size estimates.
---

Trajecta estimates brightness parameters from the photometry attached to
observations and uses them for displayed photometric quantities and magnitude
predictions. These calculations are separate from fitting the orbit.

## Asteroids: H–G model

For heliocentric distance `r` and observer distance `Δ` in AU, phase angle
`α` in radians, absolute magnitude `H` and the chosen slope parameter `G`:

```text
V = H + 5 log10(r Δ) - 2.5 log10((1-G) Φ1 + G Φ2)
Φ1 = exp(-3.33 tan(α/2)^0.63)
Φ2 = exp(-1.87 tan(α/2)^1.22)
```

`G` is supplied in **Settings → Photometry** and held fixed; its default is
0.15. For each usable magnitude, the code applies a fixed band-to-V correction
and removes distance and phase terms to obtain an individual H estimate. It
combines these with weights `1 / σmag`. This is not inverse-variance weighting.
An observation with an unsupported band or invalid geometry is skipped.

Where a magnitude uncertainty is available it is used; otherwise the parser
can infer an uncertainty from the recorded decimal precision. The photometry
routine uses 0.7 mag as a fallback when no positive uncertainty is present.
The reported photometric RMS uses squared normalized residuals divided by the
sum of squared weights. The fit reads the full residual collection when
available, so inclusion in the orbital fit and inclusion in this photometric
estimate are not equivalent.

## Comets

Total and nuclear magnitudes are fitted separately when their type is known.
For either component, the predicted magnitude is:

```text
m = H + 5 log10(Δ) + K log10(r) + PHCOF × αdeg
```

`K1` is used for total magnitude and `K2` for nuclear magnitude. The nuclear
component uses the configured phase coefficient `PHCOF` in mag/degree; the
total component uses zero. Defaults are `K1 = 10`, `K2 = 10` and `PHCOF = 0`.
The slopes and phase coefficient are fixed settings, not fitted parameters.
Measurements without a known total/nuclear type do not enter these component
estimates. Usable magnitudes are combined with `1 / σmag` weights.
Unlike the asteroid H–G estimator, the comet estimator requires a positive,
finite magnitude uncertainty after input processing. Measurements without one
are skipped.

## Size and atmospheric calculations

For an asteroid, a diameter can be estimated from H and assumed geometric
albedo `p` using `D(km) = 1329 × 10^(-H/5) / sqrt(p)`. A size range therefore
reflects the assumed albedo range, rather than a direct diameter measurement.
In atmospheric MC, **H and albedo** keeps H fixed and samples the specified
albedo and material-density distributions to obtain initial mass.

For output frames, light-time settings and other ephemeris quantities, see
[Ephemerides and close approaches](../../ephemerides-and-close-approaches/).
The implementation is in `OrbitalCore/src/Photometry.cpp` and
`OrbitalCore/api/orbital_core/Types.h`.
