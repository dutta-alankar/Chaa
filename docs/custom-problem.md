# Set up your own problem

Adding a physics setup to Chaa is one new file plus a one-line
registration. This page builds a complete example — a dense blob hit
by a wind ("blob test") — and shows every optional hook along the way.

## Step 1 — create the problem file

Create `src/problems/Blob.chpl`. A problem module needs one proc,
`setup()`, which fills the primitive array `V` over the interior domain
`DInt`:

```chapel
/* Blob.chpl — dense cloud in a supersonic wind. */
module Blob {
  use Params, Grid, State, Eos;
  use Math;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);              // cell-centre coordinates
      const r = sqrt((x - cen1)**2 + (y - cen2)**2);
      const inside = r < blastR0;                // reuse an existing knob
      V[i,j,k] = mkPrim(if inside then 10.0 else 1.0,   // rho
                        if inside then 0.0 else inVx1,  // vx1
                        0.0, 0.0,                       // vx2, vx3
                        inPrs);                         // p
    }
  }
}
```

What you have access to:

- `DInt` — the (distributed) interior domain; iterate it with `forall`.
- `V[i,j,k]` — the primitive tuple `(rho, vx1, vx2, vx3, prs)`
  (slots `IRHO, IVX1, IVX2, IVX3, IPRS`).
- coordinates from `Grid`: cell centres `x1c(i), x2c(j), x3c(k)`, faces
  `x1f(i)…`, the physical position `physPos(i,j,k)` (mapped to
  Cartesian, useful for distances in curvilinear geometry), and
  `cellVol(i,j,k)` (see how `Sedov.chpl` measures its deposit volume).
- every parameter from `Params` — either reuse generic knobs
  (`inRho`, `cen1..3`, `blastR0`, …) or add your own (step 4).

!!! tip
    Don't worry about ghost cells — boundary conditions fill them right
    after `setup()` runs. With `--eos=isothermal` your pressure is also
    overwritten by ρ·cs² automatically.

## Step 2 — register it

Add three lines to `src/Problems.chpl`:

```chapel
  import ..., CylinderFlow, Blob;          // 1: import

  proc problemInit() {
    select problem {
      ...
      when "blob" do Blob.setup();         // 2: dispatch
```

(and, only if you define the hooks below, register them in
`problemUserBC` / `problemInternalBC` — line 3.)

## Step 3 — build and run

```sh
cmake --build build
./build/bin/chaa --problem=blob --nx1=256 --nx2=128 \
   --x1min=0 --x1max=4 --x2min=-1 --x2max=1 --cen1=1 --cen2=0 \
   --blastR0=0.15 --inVx1=2.7 --inPrs=0.714 --bcX1min=inflow \
   --tstop=2 --outDt=0.2 --outFormats=hdf5
```

CMake globs `src/problems/*.chpl`, so the new file is compiled in
automatically.

## Step 4 — ship a problem parameter file

Every bundled problem carries its canonical configuration as
`src/problems/<problem>_runtime_params.ini` — give yours one too:

```ini
# src/problems/blob_runtime_params.ini
problem  = blob
geometry = cartesian
nx1 = 256
nx2 = 128
x1min = 0
x1max = 4
x2min = -1
x2max = 1
bcX1min = inflow
bcX1max = outflow-diode
inVx1 = 2.7
inPrs = 0.714
cen1 = 1
blastR0 = 0.15
tstop = 2
outFormats = hdf5
```

Then the whole run is reproducible from one line (and any key can still
be overridden on the command line, which always wins):

```sh
./build/bin/chaa --paramsFile=src/problems/blob_runtime_params.ini --nx1=512
```

## Step 5 — add simulation-specific parameters

To make, say, the blob density ratio configurable as `--blobChi=…` and
as `blobChi = …` in parameter files, add it to the **three-layer
parameter system** (command line > ini file > built-in default):

1. `src/Cli.chpl` — declare the command-line flag with an "unset"
   sentinel (`UNSET_R` for reals, `UNSET_I` for ints, `UNSET_S` for
   strings):
   ```chapel
   config const blobChi = UNSET_R;
   ```
2. `src/Params.chpl` — resolve the effective value; the second argument
   is the ini-file key (keep it identical to the flag name), the third
   is the built-in default:
   ```chapel
   const blobChi = valR(Cli.blobChi, "blobChi", 10.0);
   ```
   (`valI` / `valS` for integer and string parameters.)
3. use `blobChi` anywhere — your `setup()`, boundary hooks, validators.
4. document the default in your problem's `*_runtime_params.ini`.

Before adding a new name, check whether a generic knob already fits:
`inRho/inVx1..3/inPrs` (ambient/inflow state), `cen1..3` (a centre),
`blastR0`/`cloudRad` (a radius), `mu`, `kappa`, `gravCentral`, … —
reusing them keeps the namespace small.

## Step 6 — tracers and particles (optional, free)

- The build carries `NSCAL` passive tracer fields (default 1; set
  `-DCHAA_NSCAL=…`). Dye any region in `setup()`:
  ```chapel
  var w = mkPrim(rho, vx, vy, 0.0, p);
  if ISC < NTOT then w(ISC) = if inside then 1.0 else 0.0;
  V[i,j,k] = w;
  ```
  Tracers are advected (bounded, mass-flux-consistent), optionally
  diffused with `--scDiff`, and written as `sc0…` in every output.
- Lagrangian tracer particles need no problem code at all:
  `--nParticles=N` scatters them uniformly and writes
  `<problem>.particles.NNNN.txt` beside every dump. To seed them
  yourself, give your module a `particleInit` hook and register it in
  `Problems.problemParticleInit` (return `false` to fall back to the
  random scatter):
  ```chapel
  // in src/problems/Blob.chpl — a line of tracers across the wind
  proc particleInit(ref pos: [?D] 3*real): bool {
    for p in D do
      pos[p] = (cen1, -1.0 + 2.0*(p + 0.5)/D.size, 0.0);
    return true;
  }
  ```
  ```chapel
  // in src/Problems.chpl
  proc problemParticleInit(ref pos: [] 3*real): bool {
    select problem {
      when "blob" do return Blob.particleInit(pos);
      otherwise do return false;
    }
  }
  ```
  (`IsentropicVortex.particleInit` is a working example: with
  `--partRingR > 0` it puts the tracers on a ring around the vortex —
  validated by the `vortex-particles-ring` CI case.) Particles are
  automatically distributed across locales by position and migrate as
  they move; see [Tracer particles](user-guide/particles.md).

## Step 7 — custom forces and potentials (optional)

For an external body force or a static potential, implement an
acceleration proc in your problem module and register it in the
`problemBodyForce` dispatcher (`src/Problems.chpl`):

```chapel
// in src/problems/Blob.chpl
inline proc accel(i: int, j: int, k: int, t: real): 3*real {
  // e.g. a softened point-mass potential at (cen1, cen2):
  const dx = x1c(i) - cen1, dy = x2c(j) - cen2;
  const ir3 = 1.0/(dx*dx + dy*dy + gravEps**2)**1.5;
  return (-gravCentral*dx*ir3, -gravCentral*dy*ir3, 0.0);
}
```

```chapel
// in src/Problems.chpl
const problemHasBodyForce = problem == "blob";
proc problemBodyForce(i, j, k, t: real): 3*real {
  select problem {
    when "blob" do return Blob.accel(i, j, k, t);
    otherwise do return (0.0, 0.0, 0.0);
  }
}
```

The acceleration is applied with the matching energy source
ρ v·a in every stage. (For gravity from the *gas itself*, just enable
`--sgFourPiG`; for a fixed central mass, `--gravCentral` — no code
needed.)

## Step 8 — custom boundary conditions (optional)

For anything beyond the built-in BC types, set a side to
`--bcX2min=userdef` and provide a hook. The ghost slab is just a slice
of the padded domain; see `DoubleMach.chpl` (time-dependent inflow) and
`TaylorCouette.chpl` (no-slip rotating walls) for working templates:

```chapel
  proc userBC(side: int, t: real) {
    if side == 2 {                              // x2min
      const Dg = DAll[1-ng1..nx1+ng1, 1-ng2..0, 1..nx3];
      forall (i, j, k) in Dg {
        var w = V[i, 1-j, k];                   // mirror cell
        w(IVX2) = -w(IVX2);                     // ... and customise
        V[i,j,k] = w;
      }
    }
  }
```

Register it in `Problems.problemUserBC`. Side ids: 0/1 = x1min/max,
2/3 = x2min/max, 4/5 = x3min/max.

## Step 9 — internal (immersed) boundaries (optional)

To carve solid objects out of the flow, set `solveMask` in `setup()`
and re-impose the solid state after every stage with an `internalBC(t)`
hook — `CylinderFlow.chpl` is a complete example. Masked cells are
skipped by the update and the time-step computation.

## Step 10 — add a validated test (recommended)

1. add a line to `tests/cases.conf`:
   ```
   blob-2d | --problem=blob --nx1=128 --nx2=64 ... --tstop=1 --outFormats=hdf5
   ```
2. add a `blob_2d(outdir)` checker in `tests/validate/validate.py` and
   register it in the `CASES` dict (`load_h5`, `check(...)`,
   `finish()` are provided);
3. run `tests/run_case.sh blob-2d`, and add `blob-2d` to the matrix in
   `.github/workflows/ci.yml` so it runs on every push.

## Checklist

- [x] `src/problems/MyProblem.chpl` with `proc setup()`
- [x] import + `when "myproblem" do MyProblem.setup();` in `src/Problems.chpl`
- [x] `src/problems/myproblem_runtime_params.ini` with the canonical run
- [ ] simulation-specific knobs in `Cli.chpl` + `Params.chpl` (if needed)
- [ ] tracer dye / particles (if useful)
- [ ] `userBC` / `internalBC` / `problemBodyForce` hooks (if needed)
- [ ] a case line + validator, wired into CI
