# Idefix HD tests in Chaa

[Idefix](https://github.com/idefix-code/idefix) (Lesur et al. 2023) is
a Kokkos-based finite-volume code whose HD test suite (`test/HD/`) is a
useful cross-code reference. This page maps every Idefix HD test onto
Chaa, and documents what was ported.

## Schemes ported from Idefix

| Idefix | Chaa | notes |
|---|---|---|
| `ORDER=1` (donor cell) | `--recon=constant` | |
| `ORDER=2` PLM (`PLM_LIM`: VanLeer/MinMod/McLim) | `--recon=linear --limiter=vanleer\|minmod\|mc` | same limiter formulas |
| `ORDER=3` LimO3 | `--recon=limo3` | direct port (Čada & Torrilhon 2009, radius-of-curvature switch, positivity fallback) |
| `ORDER=4` PPM | `--recon=ppm` | direct port of the PH13/CS08/CD11 extremum-preserving limiter chain |
| `nstages=1/2/3` | `--integrator=euler\|rk2\|rk3` | same SSP RK schemes |
| `solver tvdlf/hll/hllc` | `--riemann=llf\|hll\|hllc` | (+ Chaa adds `exact`) |
| `csiso constant/userdef` | `--eos=isothermal --csIso=… [--csSlope=…]` | locally isothermal cs ∝ R_cyl^slope covers the disk tests' cs profiles |
| `[Gravity] potential central` | `--gravCentral=GM [--gravEps=…]` | point mass at the origin |
| `TDiffusion explicit constant` | `--kappa=…` | explicit thermal conduction |
| `viscosity explicit constant` | `--mu=…` | full tensor in Cartesian, τ_Rφ in cylindrical/polar |

## Test-by-test mapping

| Idefix test | status | Chaa equivalent |
|---|---|---|
| `sod` | ✅ ported | `sod-1d-cart` (CI) |
| `sod-iso` | ✅ ported | `sod-1d-iso` (CI) — validated against the exact isothermal Riemann solution |
| `MachReflection` | ✅ ported | `dmr` (CI) — identical setup: Mach-10 shock at 60°, post-shock state (8, 8.25 sin 60°, −8.25 cos 60°, 116.5), time-dependent top boundary |
| `SedovBlastWave` | ✅ ported | `sedov-3d-idefix` (CI) — γ=5/3 on [−0.5,0.5]³, radius vs similarity solution |
| `KHI` | ✅ ported | `khi-2d-iso` (CI) — same initial condition (isothermal cs=10, ±1 shear across the wavy interface) |
| `thermalDiffusion` | ✅ ported | `thermal-diffusion` (CI) — same setup (ρ=1−A sin 2πx, p=1); decay rate vs the analytic entropy-mode rate |
| `RWI-cavity` | ✅ profile ported | `disk-cavity` (CI) — same density-cavity profile, locally isothermal cs=h/√R, central gravity, exact rotational-equilibrium init; CI runs a short equilibrium-maintenance check rather than the t=1000 RWI growth (no Fargo orbital advection in Chaa) |
| `ViscousFlowPastCylinder` | ◑ equivalent physics | `cylinder-flow` (CI) — Chaa solves it on a Cartesian grid with an immersed cylinder rather than Idefix's log-radial polar mesh, since Chaa's curvilinear viscosity currently implements only τ_Rφ |
| `ViscousDisk` | ✖ not ported | needs the full viscous stress tensor in spherical coordinates |
| `FargoPlanet` | ✖ not ported | needs Fargo orbital advection and an embedded-planet potential |
| `ShearingBox` | ✖ not ported | needs shearing-periodic boundaries and tidal/Coriolis source terms |
| `VSI` | ✖ not ported | long-time spherical viscous disk run (same missing pieces as ViscousDisk) |

The unported tests all require infrastructure beyond a general-purpose
HD solver (orbital advection, shearing boxes, full curvilinear
Navier–Stokes); they are natural future work.

## Cross-code accuracy snapshot

On the isentropic-vortex accuracy test (64², one period, L1(ρ)):
`linear` 2.1×10⁻³ → `limo3` 9.5×10⁻⁴ → `ppm` 2.2×10⁻⁴ — the expected
ordering of the Idefix scheme family.
