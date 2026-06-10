# Test suite

Every case below runs in CI on each push (and locally via
`ctest --test-dir build` or `tests/run_case.sh <case>|all`), and is
validated **quantitatively** — against exact solutions, similarity
scalings, analytic steady states or strict symmetry requirements.

## Shock tubes & 1D

| case | problem | geometry | validation |
|---|---|---|---|
| `sod-1d-cart` | Sod shock tube | 1D Cartesian | L1(ρ) vs **exact Riemann solution** (1.6×10⁻³ at 400 cells) |
| `sod-1d-exact` | Sod with `--riemann=exact` | 1D Cartesian | L1(ρ) vs exact solution (1.5×10⁻³) |
| `sod-1d-iso` | isothermal Sod (Idefix `sod-iso`) | 1D Cartesian | L1(ρ) vs **exact isothermal Riemann solution**; p=ρcs² exact |
| `sod-1d-cyl` / `sod-1d-sph` | radial shock tube | 1D cyl/sph | bounds, shock presence |
| `twoblast-1d` | Woodward–Colella blast waves | 1D Cartesian | peak ρ≈6.07 at x≈0.78 |
| `thermal-diffusion` | decaying entropy mode (Idefix `thermalDiffusion`) | 1D Cartesian | decay rate vs Γ=κ(γ−1)k²/γ within 6 % |

## Sedov–Taylor blasts

| case | geometry | validation |
|---|---|---|
| `sedov-1d-sph` | 1D spherical | radius vs similarity solution (α from Kamm & Timmes), 0.4 % |
| `sedov-2d-cyl` | 2D (R,z) | sphericity R↔z < 3 % (measured 0.01 %), radius < 5 % |
| `sedov-2d-sph` | 2D (r,θ) | radius independent of θ (< 1 %) |
| `sedov-3d-cart` | 3D (x,y,z) | radius < 8 %, octant symmetry |
| `sedov-3d-sph` | 3D (r,θ,φ) | radius angle-independent to round-off |
| `sedov-3d-idefix` | 3D, γ=5/3 (Idefix `SedovBlastWave`) | radius vs similarity solution |
| `blast-2d-polar` | 2D (R,φ), off-centre | mirror symmetry < 10⁻¹⁰ |
| `blast-3d-polar` | 3D (R,φ,z) | φ & z symmetry, valid XDMF |

## Multidimensional & instabilities

| case | problem | validation |
|---|---|---|
| `riemann2d` | Lax–Liu config 3 | diagonal symmetry < 10⁻¹⁰, bounds |
| `dmr` | double Mach reflection, M=10 (≡ Idefix `MachReflection`) | Mach-stem foot position, ρmax |
| `kh` | Kelvin–Helmholtz (adiabatic) | instability growth |
| `khi-2d-iso` | Kelvin–Helmholtz (Idefix `KHI`, isothermal cs=10) | growth from the interface seed |
| `rt` | Rayleigh–Taylor (uniform gravity) | mode growth, bounds |
| `vortex` | isentropic vortex, one period | L1(ρ)=2.1×10⁻³ (64², PLM), **exact** mass conservation |
| `vortex-limo3` | same with `--recon=limo3` | L1(ρ)=9.5×10⁻⁴ |
| `vortex-ppm` | same with `--recon=ppm --integrator=rk3` | L1(ρ)=2.2×10⁻⁴ |

## Diffusion, rotation, gravity

| case | problem | validation |
|---|---|---|
| `taylor-couette` | viscous flow between rotating cylinders | steady v_φ(R) vs **analytic Couette profile**, 0.35 % |
| `cylinder-flow` | viscous flow past an immersed cylinder | exact no-slip solid, wake deficit |
| `disk-cavity` | locally isothermal Keplerian disk with a cavity (Idefix `RWI-cavity` profile) | rotational equilibrium drift < 2 % over 10 t.u., v_φ Keplerian to 3 % |

## Running locally

```sh
ctest --test-dir build -j 4                 # everything
tests/run_case.sh sedov-2d-cyl              # one case
PY=/path/to/python tests/run_case.sh all    # custom python (needs numpy, h5py)
```

Case definitions are one line each in `tests/cases.conf`; validators
live in `tests/validate/validate.py`.
