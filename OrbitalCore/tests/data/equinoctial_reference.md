# Equinoctial Reference

`equinoctial_reference.json` contains five double-precision reference calls to
an independent Fortran coordinate converter, with analytic Jacobians.
CAR units are AU and AU/day; EQU is `[a,h,k,p,q,lambda]` in AU/radians.
The central mass is solar (`rhs=1`, `center=0`, `gms=0.01720209895^2`).

The reference was generated on 2026-09-26 in an Ubuntu VM from unmodified
libraries. No fit or ephemeris propagation is needed.

The sixth, nearly circular diagnostic case (`h=1e-8,k=-2e-8`) is deliberately
not an analytic-derivative reference: the reference's pericentre-based derivatives
lose precision there. Circular and nearly circular states are covered by
independent five-point derivative, inverse and round-trip tests instead.
