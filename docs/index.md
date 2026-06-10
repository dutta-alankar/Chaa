# Chaa ☕

**Chaa** stands for **Ch**apel-based **H**ydrodynamics for
**A**strophysical **A**pplications — with an obvious playful reference
to Bengali চা (*tea*, served hot). It is a finite-volume solver for the
compressible Euler / Navier–Stokes equations, written entirely in
[Chapel](https://chapel-lang.org). It runs in **1D, 2D and 3D** on
uniform structured grids in **Cartesian, cylindrical, polar and
spherical** coordinates, and the same binary scales from a laptop to a
cluster — Chapel's distributed arrays handle all communication under
the hood.

[Get started :material-rocket-launch:](getting-started.md){ .md-button .md-button--primary }
[Set up your own problem :material-flask:](custom-problem.md){ .md-button }

## Highlights

- **Performance portable by construction** — the numerics are plain
  data-parallel `forall` loops over `StencilDist` block-distributed
  arrays. One explicit communication call (`updateFluff()`) in the whole
  code; zero message-passing logic.
- **Full geometry support** — well-balanced curvilinear source terms:
  a uniform-pressure fluid is exactly static in every geometry, and a 2D
  (R,z) Sedov blast stays spherical to 0.01 %.
- **Modern shock-capturing schemes** — donor-cell, slope-limited PLM,
  third-order LimO3 and PPM reconstruction; Rusanov, HLL, HLLC and exact
  (Godunov) Riemann solvers; SSP RK1/RK2/RK3 time stepping.
- **Physics modules** — ideal-gas and (locally) isothermal equations of
  state, explicit viscosity, explicit thermal conduction, constant and
  central point-mass gravity.
- **Output for real workflows** — ASCII tables (1D), legacy VTK, and
  HDF5 with XDMF companions that load directly into ParaView and VisIt.
- **Validated, continuously** — 26 test problems (Sod, Sedov–Taylor,
  blast waves, double Mach reflection, Kelvin–Helmholtz,
  Rayleigh–Taylor, isentropic vortex, Taylor–Couette, flow past a
  cylinder, thermal diffusion, rotating disks, …) are checked
  quantitatively against exact and similarity solutions in CI on every
  push.

## A taste

```sh
cmake -B build && cmake --build build
./build/bin/chaa --problem=sedov --geometry=spherical --nx1=512 \
                 --x1min=0 --x1max=1.2 --bcX1min=reflect --tstop=0.5
```

gives a Sedov–Taylor blast whose shock lands within 0.5 % of the
analytic similarity solution.

## Provenance

Chaa grew out of the
[`advection` branch of `1d-fluid-finite-volume`](https://github.com/dutta-alankar/1d-fluid-finite-volume/tree/advection),
a 1D Chapel finite-volume experiment, and borrows test problems and
scheme choices from the [PLUTO](https://plutocode.ph.unito.it) and
[Idefix](https://github.com/idefix-code/idefix) astrophysical codes.
MIT licensed.
