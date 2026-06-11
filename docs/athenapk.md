# AthenaPK features in Chaa

[AthenaPK](https://github.com/parthenon-hpc-lab/athenapk) is the
Athena-flavoured (M)HD code built on the Parthenon AMR framework. This
page maps its **hydro** feature set onto Chaa: what was brought in, how
to use it, and what is explicitly out of scope.

## Features ported

| AthenaPK | Chaa | notes |
|---|---|---|
| passive scalars (`nscalars`) | `-DCHAA_NSCAL=n` tracer fields | ride in the state vector; upwinded on the mass flux; bounded (validated on the Sod contact and cloud tests) |
| tracer particles | `--nParticles=…` | RK2 advection in the interpolated velocity field, any grid law, periodic wrapping; written beside every dump |
| (tabulated) cooling, Townsend integration | `--coolLambda0, --coolAlpha, --coolTfloor` | power-law Λ(T)=Λ₀Tᵅ with the **exact** Townsend (2009) integration, operator split — matches the analytic solution to round-off in the `cooling-box` test |
| turbulence driver (spectral OU forcing) | `--forceAmp, --forceTcorr, --forceKmin/Kmax, --forceSeed` | solenoidal modes with Ornstein-Uhlenbeck temporal correlation (`src/Forcing.chpl`) |
| WENO-Z reconstruction | `--recon=wenoz` | Borges et al. (2008), 5-point stencil |
| PPM, PLM, donor cell | `--recon=ppm / linear / constant` | already present |
| VL2 (predictor-corrector) integrator | `--integrator=vl2` | plus RK1/2/3 |
| viscosity / thermal conduction | `--mu`, `--kappa` | already present (explicit) |
| passive-scalar diffusion | `--scDiff` | conservative ρD∇s flux |
| diode ("outflow, no inflow") boundaries | `--bcX*=outflow-diode` (and `inflow-diode`) | |

## Test problems recreated (hydro)

| AthenaPK problem | Chaa | CI case |
|---|---|---|
| `sod` | `--problem=sod` | `sod-1d-cart` (+ exact-solution L1) |
| `linear_wave` | `--problem=linearWave` | `linear-wave` — one period with `wenoz`+`vl2` preserves the eigenmode to 0.03 % of its amplitude |
| `blast` / `sedov` | `--problem=blast/sedov` | many (all geometries) |
| `kelvin_helmholtz` | `--problem=kh` (with tracer dye) | `kh` |
| `cloud` (wind-cloud) | `--problem=cloud` | `cloud-wind` — χ=10 cloud, dyed with a tracer, diode boundaries |
| `turbulence` | `--problem=turbulence` | `turbulence-2d` — driven to v_rms ≈ 1 with a mixing dye |
| `rand_blast`, `field_loop`, `cpaw`, `precipitator`, … | not ported | MHD-dependent or covered by existing blasts |

## Static and adaptive mesh refinement

**Not ported — by design, stated plainly.** SMR/AMR is the raison
d'être of the Parthenon framework underneath AthenaPK: a forest of mesh
blocks, prolongation/restriction, flux correction at refinement
boundaries, and dynamic load balancing. Grafting that onto Chaa's
single distributed-domain architecture would be a rewrite of the mesh
layer, and a half-working AMR would compromise an otherwise carefully
validated solver. What Chaa *does* offer today for resolution control
are the **non-uniform grid laws** (logarithmic and geometric stretching
per direction — effective static refinement toward a region of
interest, validated by the `sedov-1d-log` test with a 100× spacing
range), and the architecture notes in
[Architecture & code flow](architecture.md) sketch the block-tree
design a future AMR layer would need.

## Quick recipes

```sh
# driven isothermal turbulence with a mixing dye:
./build/bin/chaa --paramsFile=src/problems/turbulence_runtime_params.ini

# wind-cloud with tracers and diode boundaries:
./build/bin/chaa --paramsFile=src/problems/cloud_runtime_params.ini \
                 --nParticles=256

# radiative blast: Sedov with power-law cooling:
./build/bin/chaa --problem=sedov --geometry=spherical --nx1=512 \
   --x1min=0 --x1max=1.2 --bcX1min=reflect --tstop=0.5 \
   --coolLambda0=0.05 --coolAlpha=0.5
```
