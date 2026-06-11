# Numerical schemes

Chaa is a method-of-lines Godunov code: reconstruct primitive variables
to faces, solve a Riemann problem per face, accumulate the
finite-volume divergence and sources, advance with a strong-stability-
preserving Runge–Kutta scheme.

## Reconstruction (`--recon=`)

| option | order | stencil | notes |
|---|---|---|---|
| `constant` | 1 | 1 cell | donor cell |
| `linear` | 2 | 3 cells | MUSCL/PLM with `--limiter=minmod\|vanleer\|mc` |
| `limo3` | 3 | 3 cells | Čada & Torrilhon (2009) limiter with the radius-of-curvature switch, as in Idefix |
| `ppm` | up to 4 (smooth) | 5 cells | parabolic face values with the extremum-preserving limiter of Colella & Sekora (2008) / Peterson & Hammett (2013); needs `NG ≥ 3` ghost layers (the default) |
| `wenoz` | 5 (smooth) | 5 cells | WENO-Z (Borges et al. 2008), as in AthenaPK; needs `NG ≥ 3` |

Measured on the isentropic-vortex accuracy test (64², one advection
period, L1 density error):

| scheme | L1(ρ) |
|---|---|
| `linear` (van Leer) + rk2 | 2.1 × 10⁻³ |
| `limo3` + rk2 | 9.5 × 10⁻⁴ |
| `ppm` + rk3 | 2.2 × 10⁻⁴ |

`limo3` and `ppm` are direct ports of the Idefix implementations
(`src/fluid/RiemannSolver/slopeLimiter.hpp`), including the positivity
fallbacks for density and pressure.

## Riemann solvers (`--riemann=`)

| option | description |
|---|---|
| `llf` | Rusanov / local Lax–Friedrichs — most diffusive, most robust |
| `hll` | two-wave HLL with Davis wave-speed estimates |
| `hllc` | HLL + restored contact wave (default) |
| `exact` | iterative Godunov solver (Toro ch. 4): Newton iteration on the star pressure, self-similar sampling at the face |

For the isothermal EOS the states carry \(p = \rho c_s^2\) and the
sound speed is \(\sqrt{p/\rho}\); the approximate solvers reduce to
their isothermal counterparts.

## Time integration (`--integrator=`)

| option | scheme | matches Idefix |
|---|---|---|
| `euler` | forward Euler | `nstages=1` |
| `rk2` | SSP RK2 (Heun) | `nstages=2` |
| `rk3` | SSP RK3 (Shu–Osher) | `nstages=3` |
| `vl2` | van Leer predictor-corrector (midpoint) | AthenaPK's `vl2` |

On the acoustic linear-wave test (128 cells, one period), `wenoz`+`vl2`
preserves the eigenmode to 3×10⁻⁴ of its amplitude.

All reconstructions evaluate their stencils with the local cell
spacing, so they run unchanged on the logarithmic and stretched grid
laws (formally the limiters keep their uniform-grid coefficients —
second-order accuracy is retained, the higher-order schemes degrade
gracefully on strongly stretched meshes).

The time step combines the advective CFL condition over all active
directions (using physical cell sizes — \(r\,\Delta\theta\),
\(R\,\Delta\phi\), …) with the explicit diffusion limits of viscosity
and conduction.

## Positivity

Density and pressure floors after every stage, plus per-face fallbacks
in the reconstruction (PLM/LimO3/PPM revert toward first order if a
face state would go negative).
