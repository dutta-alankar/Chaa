# Configuration

Chaa resolves every runtime parameter with the precedence

```
command line  (--key=value)   >   runtime_params.ini   >   built-in default
```

- the **command line** uses Chapel `config const` flags: `--nx1=256`,
  `--riemann=hllc`, …
- the **parameter file** is an INI-style file looked up as
  `runtime_params.ini` in the working directory, or anywhere via
  `--paramsFile=path`. Keys use exactly the same names as the flags;
  `[sections]` are cosmetic; `#` and `;` start comments.
- **compile-time** parameters live in
  [`src/compile_params.chpl`](https://github.com/dutta-alankar/Chaa/blob/main/src/compile_params.chpl)
  and map to CMake options.

## Run control

| key | default | meaning |
|---|---|---|
| `problem` | `sod` | which problem setup to run (see [test suite](../tests.md)) |
| `geometry` | `cartesian` | `cartesian` \| `cylindrical` \| `polar` \| `spherical` |
| `tstop` | `0.2` | stop time |
| `maxSteps` | 10⁹ | hard step limit |
| `cfl` | `0.4` | advective CFL number |
| `cflVisc` | `0.3` | safety factor of the explicit-diffusion limit |
| `dtMax` | 10³⁰ | optional dt ceiling |
| `logEvery` | `100` | console log cadence (steps) |
| `paramsFile` | `runtime_params.ini` | parameter file location |

## Grid

| key | default | meaning |
|---|---|---|
| `nx1, nx2, nx3` | `128, 1, 1` | cells per direction (1 = inactive dimension) |
| `x1min, x1max` … | `0, 1` | domain extent per direction |
| `gridX1, gridX2, gridX3` | `uniform` | grid law per direction: `uniform` \| `log` (spacing grows ∝x, needs xmin>0) \| `log-dec` (spacing shrinks with x) \| `stretch` (geometric progression) |
| `stretchX1..3` | `1.05` | spacing ratio for `stretch` (>1 grows, <1 shrinks) |

Coordinate meanings per geometry are listed in
[equations & geometry](equations.md).

## Physics

| key | default | meaning |
|---|---|---|
| `gam` | `1.4` | adiabatic index (ideal EOS) |
| `eos` | `ideal` | `ideal` \| `isothermal` |
| `csIso` | `1.0` | isothermal sound speed |
| `csSlope` | `0.0` | locally isothermal: cs = csIso·R_cyl^csSlope |
| `mu` | `0` | dynamic viscosity (explicit) |
| `kappa` | `0` | thermal conductivity (explicit, ideal EOS only) |
| `grav1..3` | `0` | constant gravity along the coordinate axes |
| `gravCentral` | `0` | GM of a central point mass at the origin |
| `gravEps` | `0` | gravitational softening length |
| `rhoFloor`, `prsFloor` | 10⁻¹², 10⁻¹⁴ | positivity floors |
| `coolLambda0, coolAlpha, coolTfloor` | `0, 0.5, 10⁻⁶` | optically thin cooling Λ(T)=Λ₀Tᵅ (exact Townsend integration) |
| `scDiff` | `0` | passive-scalar diffusivity |
| `forceAmp, forceTcorr, forceKmin, forceKmax, forceSeed` | off | Ornstein-Uhlenbeck turbulence driving |
| `nParticles, partSeed` | `0` | Lagrangian tracer particles |

## Numerics

| key | default | options |
|---|---|---|
| `recon` | `linear` | `constant`, `linear`, `limo3`, `ppm`, `wenoz` |
| `limiter` | `vanleer` | `minmod`, `vanleer`, `mc` (PLM only) |
| `riemann` | `hllc` | `llf`, `hll`, `hllc`, `exact` |
| `integrator` | `rk2` | `euler`, `rk2`, `rk3`, `vl2` |

## Boundary conditions

One per side: `bcX1min`, `bcX1max`, `bcX2min`, `bcX2max`, `bcX3min`,
`bcX3max`, each one of

| value | meaning |
|---|---|
| `zero-gradient` | copy the nearest interior cell (default; `outflow` is a legacy alias) |
| `outflow-diode` | zero-gradient, but the normal velocity is clamped so nothing can flow *in* |
| `inflow-diode` | zero-gradient, but the normal velocity is clamped so nothing can flow *out* |
| `periodic` | periodic wrap |
| `reflect` | mirror, normal velocity flipped |
| `axis` | reflect + azimuthal velocity flipped (use at r=0, θ=0, θ=π) |
| `inflow` | fixed state from `inRho, inVx1..3, inPrs` |
| `userdef` | problem-defined hook (see [custom problems](../custom-problem.md)) |

## Output

| key | default | meaning |
|---|---|---|
| `outFormats` | `txt` | comma list among `txt`, `vtk`, `hdf5` |
| `outDt` | `0` | dump interval (≤0: only initial and final states) |
| `outDir` | `output` | output directory (created if missing) |
| `restartDt` | `0` | restart-file cadence in simulation time (0: only at the end of the run and on a graceful stop) — see [Stopping & restarting](restart.md) |
| `restart` | `false` | *(command line only)* resume from `<outDir>/restart.chaa` |

## Problem-specific parameters

Each bundled problem reads its own knobs (all also settable from the
ini file): `sodX0`, `sodRhoL/R`, `sodVxL/R`, `sodPrsL/R`;
`sedovE0`, `sedovR0`, `sedovRhoAmb`, `sedovPrsAmb`; `blastPin/Pout`,
`blastRhoIn/Out`, `blastR0`; `cen1..3` (explosion / object centre);
`inRho`, `inVx1..3`, `inPrs` (inflow state); `tcOmegaIn/Out`
(Taylor–Couette wall rotation); `cylRad` (cylinder radius);
`vortexBeta`; `khRhoIn/Out`, `khV0`, `khPert`, `khPrs`; `rtRhoTop/Bot`,
`rtPrs0`, `rtPert`; `diskH0`, `diskJumpR`, `diskJumpW`; `twAmp`; `cloudChi`, `cloudRad`;
`waveAmp`.

Every bundled problem also ships its canonical configuration as
`src/problems/<problem>_runtime_params.ini`, runnable directly:

```sh
./build/bin/chaa --paramsFile=src/problems/cloud_runtime_params.ini
```

## Compile-time parameters

| CMake option | Chapel param | default | meaning |
|---|---|---|---|
| `-DCHAA_NG=n` | `NG` | `3` | ghost layers (2 suffices below `ppm`) |
| `-DCHAA_HDF5=ON/OFF` | `hdf5Enabled` | `ON` | HDF5 writer compiled in |
| `-DCHAA_NSCAL=n` | `NSCAL` | `1` | passive tracer fields carried in the state vector |
| `-DCHAA_CHPL_FLAGS=…` | — | `--fast` | extra `chpl` flags |
