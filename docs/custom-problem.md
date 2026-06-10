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
  use Params, Grid, State;
  use Math;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);             // cell-centre coordinates
      const r = sqrt((x - cen1)**2 + (y - cen2)**2);
      const inside = r < blastR0;               // reuse an existing knob
      V[i,j,k] = (if inside then 10.0 else 1.0, // rho
                  if inside then 0.0 else inVx1, // vx1
                  0.0, 0.0,                      // vx2, vx3
                  inPrs);                        // p
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

## Step 4 — add your own parameters (optional)

To make, say, the blob density ratio configurable as `--blobChi=…` and
from `runtime_params.ini`:

1. `src/Cli.chpl` — declare the sentinel flag:
   ```chapel
   config const blobChi = UNSET_R;
   ```
2. `src/Params.chpl` — resolve it (command line > ini > default):
   ```chapel
   const blobChi = valR(Cli.blobChi, "blobChi", 10.0);
   ```
3. use `blobChi` in your `setup()`.

## Step 5 — custom boundary conditions (optional)

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

## Step 6 — internal (immersed) boundaries (optional)

To carve solid objects out of the flow, set `solveMask` in `setup()`
and re-impose the solid state after every stage with an `internalBC(t)`
hook — `CylinderFlow.chpl` is a complete example. Masked cells are
skipped by the update and the time-step computation.

## Step 7 — add a validated test (recommended)

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
- [ ] new knobs in `Cli.chpl` + `Params.chpl` (if needed)
- [ ] `userBC` / `internalBC` hooks (if needed)
- [ ] a case line + validator, wired into CI
