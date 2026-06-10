# Chaa ☕

[![CI](https://github.com/dutta-alankar/Chaa/actions/workflows/ci.yml/badge.svg)](https://github.com/dutta-alankar/Chaa/actions/workflows/ci.yml)

**Chaa** (চা — *tea*, served hot) is a finite-volume hydrodynamics solver
written in [Chapel](https://chapel-lang.org). It solves the compressible
Euler / Navier–Stokes equations (continuity, momentum, energy) on uniform
structured grids in **1D, 2D and 3D**, in **Cartesian, cylindrical, polar
and spherical** coordinates.

It grew out of (and takes its design inspiration from) the
[`advection` branch of `1d-fluid-finite-volume`](https://github.com/dutta-alankar/1d-fluid-finite-volume/tree/advection),
generalising that 1D Chapel advection solver to the full Euler system in
curvilinear geometry.

## Why Chapel?

The entire solver is written as plain data-parallel `forall` loops over
block-distributed arrays (`StencilDist`). Chapel's distributed-array
abstraction keeps every halo exchange and remote access **under the hood**:
there is exactly one explicit communication call in the code
(`updateFluff()`, which refreshes the stencil halo caches), and *zero*
message-passing logic. The same binary runs

- multithreaded on a laptop (`CHPL_COMM=none`),
- distributed across nodes (`CHPL_COMM=gasnet`, launched with `-nl N`),

with no change to the numerics — that is the performance portability the
code is built around.

## Numerics

| ingredient        | options                                                          |
|-------------------|------------------------------------------------------------------|
| reconstruction    | piecewise constant, piecewise linear (`minmod`, `vanleer`, `mc`) |
| Riemann solver    | Rusanov (`llf`), `hll`, `hllc`, `exact` (iterative Godunov)      |
| time integration  | `euler`, SSP `rk2`, SSP `rk3`                                    |
| source terms      | well-balanced curvilinear geometry, uniform gravity              |
| diffusion         | explicit viscosity: full stress tensor (Cartesian), τ<sub>Rφ</sub> (cylindrical/polar) |
| positivity        | density/pressure floors, slope-limited face states               |

The geometric source terms are discretised with the same area/volume
factors as the flux divergence, so a uniform-pressure fluid is balanced to
machine precision in every geometry (no spurious accelerations near axes).

### Coordinates (PLUTO conventions)

| `--geometry=`  | x1 | x2 | x3 | notes                                       |
|----------------|----|----|----|---------------------------------------------|
| `cartesian`    | x  | y  | z  |                                             |
| `cylindrical`  | R  | z  | φ  | axisymmetric: `nx3 = 1`, v<sub>φ</sub> evolves passively |
| `polar`        | R  | φ  | z  |                                             |
| `spherical`    | r  | θ  | φ  |                                             |

Any subset of dimensions may be active (`nx2 = nx3 = 1` gives 1D x/R/r,
etc.), covering 1D (x / r), 2D (x,y / R,z / R,φ / r,θ) and 3D
(x,y,z / R,φ,z / r,θ,φ).

## Building

Requires CMake ≥ 3.16, Chapel ≥ 2.8 (`brew install chapel`, or the
[`chapel/chapel` docker image](https://hub.docker.com/r/chapel/chapel)) and,
for HDF5 output, libhdf5 (`brew install hdf5` / `apt install libhdf5-dev`).

```sh
cmake -B build                  # -DCHAA_HDF5=OFF to drop the HDF5 dependency
cmake --build build             # -> build/bin/chaa
ctest --test-dir build          # run the validated test suite
```

Compile-time parameters (ghost-layer count, HDF5 support) live in
[`src/compile_params.chpl`](src/compile_params.chpl) and map to CMake
options (`-DCHAA_NG=…`, `-DCHAA_HDF5=…`).

## Running

Runtime options are resolved with the precedence

```
command line  (--key=value)   >   runtime_params.ini   >   built-in default
```

[`runtime_params.ini`](runtime_params.ini) (looked up in the working
directory, or anywhere via `--paramsFile=…`) holds the run configuration;
any key can still be overridden per-run on the command line:

```sh
./build/bin/chaa --problem=sod --nx1=400 --tstop=0.2 --outFormats=txt,vtk,hdf5
./build/bin/chaa --problem=sedov --geometry=spherical --nx1=512 \
                 --x1min=0 --x1max=1.2 --bcX1min=reflect --tstop=0.5
./build/bin/chaa --paramsFile=my_run.ini

# multi-locale (with CHPL_COMM=gasnet build):
./build/bin/chaa -nl 4 --problem=sedov --geometry=cartesian \
                 --nx1=256 --nx2=256 --nx3=256
```

Boundary conditions per side (`--bcX1min=…` etc.): `periodic`, `outflow`,
`reflect`, `axis`, `inflow`, `userdef`.

### Output formats (`--outFormats=`)

- `txt` — column ASCII (1D)
- `vtk` — legacy VTK; rectilinear grid for Cartesian, structured grid
  mapped to physical space for curvilinear meshes (wedges, annuli, shells
  render correctly in ParaView/VisIt)
- `hdf5` — HDF5 datasets plus an `.xmf` XDMF companion that loads
  directly in ParaView and VisIt (2D/3D)

## Test problems

All of the classic HD test problems (largely following the PLUTO test
suite) are built in, validated quantitatively in CI on every push:

| case | problem | geometry | validation |
|------|---------|----------|------------|
| `sod-1d-cart` | Sod shock tube | 1D Cartesian | L1(ρ) vs **exact Riemann solution** < 1.2 % |
| `sod-1d-exact` | Sod with the `exact` Godunov solver | 1D Cartesian | L1(ρ) vs exact solution |
| `sod-1d-cyl`, `sod-1d-sph` | radial shock tube | 1D cyl/sph | bounds, shock presence |
| `twoblast-1d` | Woodward–Colella interacting blasts | 1D Cartesian | peak ρ ≈ 6 at x ≈ 0.78 |
| `sedov-1d-sph` | Sedov–Taylor | 1D spherical | shock radius vs similarity solution (α from Kamm & Timmes), err < 4 % |
| `sedov-2d-cyl` | Sedov–Taylor | 2D (R,z) | sphericity R↔z < 3 %, radius < 5 % |
| `sedov-2d-sph` | Sedov–Taylor | 2D (r,θ) | radius θ-independent (< 1 %) |
| `sedov-3d-cart` | Sedov–Taylor | 3D (x,y,z) | radius < 8 %, octant symmetry |
| `sedov-3d-sph` | Sedov–Taylor | 3D (r,θ,φ) | radius angle-independent |
| `blast-2d-polar` | blast wave | 2D (R,φ) | mirror symmetry < 1e-10 |
| `blast-3d-polar` | blast wave | 3D (R,φ,z) | φ & z symmetry, XDMF validity |
| `riemann2d` | Lax–Liu config 3 | 2D Cartesian | diagonal symmetry < 1e-10, bounds |
| `dmr` | double Mach reflection (Mach 10) | 2D Cartesian | Mach-stem foot position, ρ<sub>max</sub> |
| `kh` | Kelvin–Helmholtz | 2D Cartesian | instability growth |
| `rt` | Rayleigh–Taylor (gravity) | 2D Cartesian | mode growth, bounds |
| `vortex` | isentropic vortex | 2D Cartesian | L1(ρ) after a full period < 8e-3, exact mass conservation |
| `taylor-couette` | Taylor–Couette (viscous) | 1D cylindrical | steady v<sub>φ</sub>(R) vs **analytic Couette profile** < 2 % |
| `cylinder-flow` | viscous flow past a cylinder | 2D Cartesian | no-slip solid, wake deficit |

Run them locally:

```sh
ctest --test-dir build -j 4        # or: tests/run_case.sh all | <case-name>
```

(validation needs `python3` with `numpy` and `h5py`).

Representative measured accuracies (Apple Silicon, Chapel 2.8):
Sod L1(ρ) = 1.6e-3 at 400 cells; Sedov shock radius within 0.5 % of the
similarity solution in 1D spherical and 0.01 % aspherical in (R,z);
Taylor–Couette steady profile within 0.35 % of analytic; 3D spherical
Sedov shock radius independent of angle to machine precision.

## Layout

```
runtime_params.ini        runtime configuration file (overridden by CLI)
src/compile_params.chpl   compile-time parameters (config params)
src/Cli.chpl              command-line layer (config consts + sentinels)
src/IniReader.chpl        runtime_params.ini parser
src/Params.chpl           effective parameters (CLI > ini > default)
src/Grid.chpl             closed-form mesh + curvilinear metric factors
src/State.chpl            StencilDist-distributed field arrays
src/Eos.chpl              gamma-law EOS, prim<->cons
src/Recon.chpl            slope-limited reconstruction
src/Riemann.chpl          LLF / HLL / HLLC / exact
src/Hydro.chpl            sweeps, geometric sources, viscosity, gravity, CFL
src/Boundary.chpl         ghost-cell boundary conditions
src/Problems.chpl         problem registry/dispatcher
src/problems/*.chpl       one file per test problem (IC + user BCs)
src/Evolve.chpl           SSP Runge-Kutta stepping
src/Output.chpl           txt / VTK / XDMF writers
src/Hdf5IO.chpl           minimal HDF5 bindings + writer
src/Chaa.chpl             driver
```

Adding a problem = one new file in `src/problems/` with a `setup()` proc
plus a one-line registration in `src/Problems.chpl`.

## License

MIT — see [LICENSE](LICENSE).
