# Architecture & code flow

This page walks through Chaa's code logic in execution order: what
happens at startup, how an integration step flows through the source
files, and which module owns which responsibility.

## Source map

| file | owns | depends on |
|---|---|---|
| `src/compile_params.chpl` | compile-time knobs: `NG` (ghost layers), `hdf5Enabled`, `NSCAL` (tracer fields) | — |
| `src/Cli.chpl` | every runtime flag as a `config const` with an "unset" sentinel | — |
| `src/IniReader.chpl` | parses `runtime_params.ini` / `--paramsFile` into a key→value map | Cli |
| `src/Params.chpl` | *effective* parameters (CLI > ini > default), option parsing into integer codes, state-vector layout (`NVAR`, `NSCAL`, `NTOT`, slot indices) | Cli, IniReader |
| `src/Grid.chpl` | closed-form coordinates under the four grid laws, metric factors (areas, volumes, centroids), physical cell sizes, the inverse coordinate maps | Params |
| `src/State.chpl` | the distributed arrays: `V` (primitives), `U` (conservatives), `U0`, `RHS`, `FLX`, `solveMask`, all on one `StencilDist` domain | Params, Grid |
| `src/Eos.chpl` | prim↔cons conversion, sound speed (ideal/isothermal), physical flux, `mkPrim` | Params |
| `src/Recon.chpl` | face reconstruction: constant/PLM/LimO3/PPM/WENO-Z | Params |
| `src/Riemann.chpl` | LLF/HLL/HLLC/exact solvers + tracer upwinding | Params, Eos |
| `src/Hydro.chpl` | `computeRHS`: sweeps, divergence, sources, diffusion, CFL | Grid, State, Eos, Recon, Riemann, Forcing |
| `src/Evolve.chpl` | SSP RK / VL2 staging, floors, operator-split cooling | Hydro, Boundary |
| `src/Boundary.chpl` | ghost-cell boundary conditions, halo refresh | State, Eos, Problems |
| `src/Problems.chpl` + `src/problems/*` | one module per test problem: `setup()`, optional `userBC`, `internalBC`; one `*_runtime_params.ini` per problem | Params, Grid, State, Eos |
| `src/Forcing.chpl` | Ornstein-Uhlenbeck spectral turbulence driving | Params, Grid |
| `src/Particles.chpl` | Lagrangian tracer particles | Params, Grid, State |
| `src/Output.chpl` | txt/VTK writers, XDMF, parallel piece orchestration | State, Hdf5IO, Particles |
| `src/Hdf5IO.chpl` | minimal HDF5 bindings + block writer | Grid, State |
| `src/Chaa.chpl` | the driver: banner, sanity checks, main loop | everything |

## Startup sequence

Chapel initialises modules before `main()` runs, in dependency order —
this is how the three-level parameter system resolves with no explicit
orchestration:

```
Cli            config consts get command-line values (or sentinels)
  └─ IniReader  module-init statement loadIni() parses the parameter file
       └─ Params   each `const` resolves: CLI value if set, else ini, else default
            └─ Grid/State  domains and distributed arrays are created
```

`main()` (`src/Chaa.chpl`) then:

1. prints the logo and configuration banner; runs `sanityChecks()`
   (geometry/EOS/scheme compatibility);
2. `problemInit()` dispatches to the selected problem's `setup()`,
   which fills the primitive array `V` over the interior `DInt`;
3. for the isothermal EOS, `p = rho cs^2(x)` is enforced;
4. `U = prim2cons(V)` over the interior; `problemInternalBC` applies
   immersed solids; `applyBCs(0)` fills ghosts and refreshes halos;
5. `initForcing()` / `initParticles()` if enabled;
6. writes dump 0 and enters the main loop.

## The main loop

```
while t < tstop:
    dt = computeDt()                      # Hydro: CFL + diffusion limits
    updateForcing(dt)                     # Forcing: OU amplitude update
    advance(dt, t)                        # Evolve: the integrator
    advanceParticles(dt)                  # Particles: RK2 in frozen field
    [writeOutputs at the outDt cadence]   # Output
```

### Inside `advance` (src/Evolve.chpl)

Each SSP Runge-Kutta stage has the canonical form
`U <- cA*U0 + cB*(U + dt*RHS)` and runs:

```
stage:
    computeRHS(tStage)          (Hydro)
    update U over DInt          (skipping solveMask-ed cells)
    applyFloorsAndPrims()       cons->prim with floors; isothermal p reset
    problemInternalBC(tStage)   re-impose immersed solids
    applyBCs(tStage)            (Boundary) ghosts + updateFluff()
```

`euler` is one stage, `rk2` two, `rk3` three (Shu–Osher), `vl2` is a
half-step predictor + full-step corrector. If cooling is enabled, an
operator-split exact (Townsend) integration of `dT/dt = -(γ-1)ρΛ(T)`
runs after the stages.

### Inside `computeRHS` (src/Hydro.chpl)

```
zero RHS
for each active direction d:                       # sweepDir(d)
    forall faces:                                  #   distributed
        reconstruct wL,wR        (Recon: 3- or 5-cell stencil)
        F = riemannFlux(wL,wR,d) (Riemann; tracers upwinded on F_rho)
        F -= viscous, conduction, scalar-diffusion fluxes (if enabled)
        FLX[face] = F
    FLX.updateFluff()                              #   halo refresh
    forall cells:
        RHS += (A_left*F_left - A_right*F_right) * invVol   # FV divergence
addSources(t):                                     # cell-centred, forall
    well-balanced curvilinear terms, gravity (constant + central),
    OU forcing, residual curvilinear viscous term
```

The geometric area/volume factors (`fA1`, `invV1`, `g2`, …) come from
`Grid` as closed-form functions of the index under the active grid law,
so the same sweep code runs on uniform, log and stretched grids in
every geometry.

### Data-parallelism and communication

All fields live on one block-distributed, halo-padded domain
(`State.DAll`, a `StencilDist`). Every loop above is a `forall` that
Chapel runs across locales; reads like `V[idx-e]` near a partition edge
are served from the locale-local halo cache. The **only** explicit
communication calls in the code are `updateFluff()` — after boundary
application (`Boundary.syncHalos`) and after the face-flux loop. Ghost
cells at physical boundaries are *owned* elements of the padded domain,
so boundary conditions are ordinary distributed foralls.

## Boundary conditions (src/Boundary.chpl)

Sides are applied in x1, x2, x3 order, each pass spanning the full
extent of the other directions (corners end up consistent). Each side
is a thin-slab `forall` doing one of: zero-gradient copy, diode
(zero-gradient + one-way normal-velocity clamp), periodic copy, mirror
(±sign flips for reflect/axis), fixed inflow state, or a problem hook
(`userdef` → `Problems.problemUserBC`). Afterwards ghost conservatives
are rebuilt from the ghost primitives and halos refresh.

## Output paths (src/Output.chpl, src/Hdf5IO.chpl)

- **single locale**: one txt/VTK/HDF5(+XDMF) file per dump over the
  full interior;
- **multiple locales**: each locale concurrently writes its own block
  (`DInt.localSubdomain()`) as an independent HDF5/VTK piece file, and
  locale 0 writes one `.xmf` *spatial collection* that stitches the
  pieces seamlessly in ParaView/VisIt.

## What is deliberately *not* here

Chaa is a single-mesh code. **Static and adaptive mesh refinement are
out of scope** for the current architecture: block-structured AMR needs
a tree of mesh blocks with prolongation/restriction operators, flux
correction at fine-coarse faces and dynamic load balancing — the very
machinery the [Parthenon](https://github.com/parthenon-hpc-lab/parthenon)
framework (on which AthenaPK is built) provides. A Chapel-native
equivalent would replace `State`'s single `StencilDist` domain with a
collection of per-block domains plus an oct-tree index — a
rearchitecture, not a feature; it is tracked as future work rather than
shipped half-done. The same applies to orbital advection (Fargo) and
shearing boxes.
