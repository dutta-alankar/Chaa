# Physics modules

Every module below is enabled by a single runtime parameter and
composes freely with the others (and with any geometry/grid law unless
stated). Each one carries a quantitative CI test, quoted with its
measured result.

## Equation of state

- **Ideal gas** (default): γ-law, `--gam`.
- **Isothermal / locally isothermal**: `--eos=isothermal --csIso=…
  [--csSlope=…]`; p = ρcs² is enforced after every stage, with
  cs ∝ R_cyl^slope for disks (csSlope = −½ ⇒ constant aspect ratio).
  *Validated against the exact isothermal Riemann solution
  (`sod-1d-iso`).*

## Gravity

- **Constant**: `--grav1/2/3` — acceleration along the coordinate axes
  (`rt` test).
- **Central point mass**: `--gravCentral=GM [--gravEps]` — Keplerian
  disks in any geometry (`disk-cavity`: equilibrium drift 0.12 % over
  10 t.u.).
- **Self-gravity**: `--sgFourPiG=4πG [--sgTol --sgMaxIter]` — solves
  ∇²Φ = 4πGρ on the **evolving density** once per step with a
  conjugate-gradient iteration built on the same finite-volume metric
  operators as the flux divergence (any geometry/grid law). Fully
  periodic domains subtract the mean density (Jeans swindle);
  non-periodic sides impose Φ=0 ghosts. The potential warm-starts from
  the previous step. *`selfgrav-kick` verifies the acceleration field
  against the analytic Poisson solution of a sinusoidal density to
  0.08 %.*
- **Custom forces/potentials**: the `problemBodyForce(i,j,k,t)` hook —
  see [Set up your own problem](../custom-problem.md).

## Rotation: shearing box and FARGO

- **Shearing box**: `--omegaRot=Ω --shearQ=q` adds the Coriolis and
  tidal sources of the local rotating frame (Cartesian; x1 = radial,
  x2 = azimuthal), with **shear-periodic** radial boundaries
  (`--bcX1min/max=shear-periodic`) that wrap the periodic image with
  the time-dependent azimuthal offset and velocity jump.
  *`epicycle-shearbox`: a uniform radial kick oscillates at the
  epicyclic frequency κ=√(2(2−q))Ω — measured ⟨vx⟩(π/κ) = −waveAmp to
  five digits.*
- **FARGO orbital advection**: `--fargo=on` splits the azimuthal
  velocity into a background w(R) (Keplerian √(GM/R) in polar, −qΩx in
  the shearing box) plus a residual: the azimuthal Riemann solve runs
  in the comoving frame (fluxes transformed back exactly), the CFL
  condition sees only the *residual* velocity, and a conservative
  slope-limited remap shifts each radial row by w·dt after every step.
  Cold Keplerian disks step ~v_K/cs times faster. *`epicycle-fargo`
  reproduces the epicyclic frequency exactly; `disk-cavity-fargo` holds
  the same rotational equilibrium as the non-FARGO run.*

## Dissipation and cooling

- **Viscosity** `--mu`: full stress tensor in Cartesian, the dominant
  τ_Rφ term in cylindrical/polar (*Taylor–Couette steady profile to
  0.35 %*), with viscous heating.
- **Thermal conduction** `--kappa`: q = −κ∇T along every active
  direction (*entropy-mode decay rate to 3 %*, cross-checked against
  Idefix to the same level).
- **Passive-scalar diffusion** `--scDiff`: conservative ρD∇s flux.
- **Optically thin cooling** `--coolLambda0 --coolAlpha --coolTfloor`:
  Λ(T) = Λ₀Tᵅ integrated **exactly** over each step (Townsend 2009),
  operator split, unconditionally stable (*matches the analytic
  power-law solution to round-off in `cooling-box`*).
- All explicit diffusion terms share the time-step limit
  Δt ≤ ½·cflVisc·Δℓ²/(χ n_dim).

## Turbulence driving

`--forceAmp --forceTcorr --forceKmin/Kmax --forceSeed`: solenoidal
spectral forcing with Ornstein-Uhlenbeck temporal correlation
(`src/Forcing.chpl`), AthenaPK-style. *`turbulence-2d` drives an
isothermal box to v_rms ≈ 1 with a mixing tracer dye.*

## Tracers

- **Tracer fields**: `-DCHAA_NSCAL=n` passive scalars in the state
  vector, advected with mass-flux-consistent upwinding (provably
  bounded — the Sod dye rides the contact exactly in [0,1]).
- **Tracer particles**: see [Particles](particles.md).

## Composition matrix

| | cart | cyl | polar | sph | log/stretch grids |
|---|---|---|---|---|---|
| isothermal EOS | ✓ | ✓ | ✓ | ✓ | ✓ |
| central gravity | ✓ | ✓ | ✓ | ✓ | ✓ |
| self-gravity | ✓ | ✓ | ✓ | ✓ | ✓ |
| viscosity | ✓ | τ_Rφ | τ_Rφ | — | ✓ |
| conduction / cooling / scalars | ✓ | ✓ | ✓ | ✓ | ✓ |
| shearing box + shear-periodic | ✓ | — | — | — | x2 uniform |
| FARGO | ✓ (SB) | — | ✓ | — | x2 uniform |
| turbulence driving | ✓ | — | — | — | ✓ |
