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
  third-order LimO3 and PPM, and WENO-Z reconstruction; Rusanov, HLL,
  HLLC and exact (Godunov) Riemann solvers; SSP RK1/RK2/RK3 and VL2
  time stepping.
- **Physics modules** — ideal-gas and (locally) isothermal equations of
  state, explicit viscosity, thermal conduction, optically thin cooling
  (exact Townsend integration), constant, central point-mass and
  Poisson **self-gravity**, a custom body-force hook, **shearing box**
  with shear-periodic boundaries, **FARGO orbital advection**, OU
  turbulence driving, passive tracer fields and fully distributed
  Lagrangian tracer particles.
- **Flexible meshes** — uniform, logarithmic and geometrically
  stretched (with uniform anchor blocks) grid laws per direction, all
  with closed-form metrics.
- **Output for real workflows** — ASCII tables (1D), legacy VTK, and
  HDF5 with XDMF companions that load directly into ParaView and VisIt;
  bundled python tools for
  [field plots and analytic comparisons](user-guide/plotting.md).
- **Validated, continuously** — 45 test cases (Sod, Sedov–Taylor,
  blast waves, double Mach reflection, Kelvin–Helmholtz,
  Rayleigh–Taylor, isentropic vortex, Taylor–Couette, flow past a
  cylinder, thermal diffusion, rotating disks, shearing boxes,
  self-gravity, …) are checked quantitatively against exact and
  similarity solutions in CI on every push, and matched-configuration
  runs are [cross-validated against Idefix and
  AthenaPK](cross-validation.md).
- **Scales, measurably** — near-ideal thread scaling on a node and
  ~90 % strong-scaling efficiency of the distributed code path at 4
  locales; see [Benchmarks & scaling](benchmarks.md).

## A taste

```sh
cmake -B build && cmake --build build
./build/bin/chaa --problem=sedov --geometry=spherical --nx1=512 \
                 --x1min=0 --x1max=1.2 --bcX1min=reflect --tstop=0.5
```

gives a Sedov–Taylor blast whose shock lands within 0.5 % of the
analytic similarity solution:

![Sedov density profile on the analytic similarity shock radius](assets/plots/sedov1d-radius.png)

(plotted with the bundled `tools/plot_compare.py` — every test problem
has [a figure like this](test-problems.md).)

## Provenance

Chaa grew out of the
[`advection` branch of `1d-fluid-finite-volume`](https://github.com/dutta-alankar/1d-fluid-finite-volume/tree/advection),
a 1D Chapel finite-volume experiment, and borrows test problems and
scheme choices from the [PLUTO](https://plutocode.ph.unito.it) and
[Idefix](https://github.com/idefix-code/idefix) astrophysical codes.
MIT licensed.
