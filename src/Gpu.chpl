/* Gpu.chpl — optional GPU execution engine (compile-time opt-in).
 *
 * Build with a Chapel configured for GPUs (CHPL_LOCALE_MODEL=gpu) and
 * -sgpuEnabled=true (CMake: -DCHAA_GPU=ON); the CPU build is untouched
 * when the flag is off.
 *
 * Design: the hot loops in Hydro/Evolve are written once as generic
 * kernel procs.  In a GPU build each locale carves its x1 slab of the
 * grid into one padded block per visible GPU (with co-locales each
 * locale sees exactly one GPU, so multi-node runs keep the x1-only
 * locale split).  A block holds device-resident copies of V/U/U0/RHS/
 * FLX and the whole RK stage runs as device kernels; only the thin
 * ghost/halo shell is staged through the host StencilDist arrays each
 * stage, so every boundary condition (including user-defined ones)
 * runs unchanged on the host.  Output, restart, particles, FARGO and
 * the self-gravity CG solve keep their host implementations and pull
 * the interior off the devices when needed.
 */
module Gpu {
  use Params, Grid, State, Eos;
  use List, RangeChunk, BlockDist, GpuDiagnostics;
  import CompileParams.gpuEnabled;
  import Forcing;
  import SelfGravity;
  import Boundary;
  import Problems;

  /* one x1 slab of the domain, resident on one GPU */
  class GpuBlock {
    const dev: locale;             // the GPU sublocale
    const lo, hi: int;             // interior x1 range of this block
    const DA:  domain(3);          // padded (ghost/halo shell included)
    const DI:  domain(3);          // interior
    const DF1: domain(3);          // face domains (face = left face)
    const DF2: domain(3);
    const DF3: domain(3);
    const DSG: domain(3);          // PHI domain (empty without self-gravity)

    var V, U, U0, RHS, FLX: [DA] StateVec;
    var WL, WR: [DA] StateVec;     // face states (two-pass flux kernels)
    var mask: [DA] bool;
    var DTB: [DI] real;            // per-cell dt buffer for the reduction
    var PHI: [DSG] real;           // device copy of the potential

    /* pack buffers for the strided x2/x3 ghost slabs: a slab is packed
       into contiguous device memory by a kernel and moved with one DMA
       (copying the strided slices directly degenerates into thousands
       of tiny transfers per stage) */
    const DPk2: domain(1);
    const DPk3: domain(1);
    var pk2: [DPk2] StateVec;
    var pk3: [DPk3] StateVec;

    // forcing mode tables (device copies; fa1/fa2 refreshed every step)
    var fkv, fe1, fe2: [0..#Forcing.maxModes] 3*real;
    var fph, fa1, fa2: [0..#Forcing.maxModes] real;
  }

  const LocD = {0..#numLocales} dmapped new blockDist({0..#numLocales});
  var locBlocks: [LocD] list(unmanaged GpuBlock);
  var nGpuTotal = 0;

  /* --gpuTime=true: per-phase wall-time accumulators (printed at the
     end of the run) to locate where a GPU step spends its time */
  config const gpuTime = false;
  var tDown, tHostBC, tUp: real;      // gpuStageBCs phases
  var tRHS, tKern, tDt: real;         // device-kernel phases
  use Time;

  /* ---- flattened kernel iteration -----------------------------------
     Chapel 2.8 parallelises a GPU forall only over a 1D iteration
     space: a forall over a 3D domain runs ~50-100x below memory
     bandwidth (measured 12 vs 1170 GB/s on an A100).  Every device
     kernel therefore iterates 0..#D.size and recovers (i,j,k) from
     the flat index with these helpers (a record of plain ints — a
     kernel must never query a domain object, see packKernel). */
  record FlatDom {
    var lo1, lo2, lo3: int;
    var n23, n3v: int;
    var size: int;
  }

  proc mkFlat(D): FlatDom {
    var fd: FlatDom;
    fd.lo1 = D.dim(0).low;
    fd.lo2 = D.dim(1).low;
    fd.lo3 = D.dim(2).low;
    fd.n3v = D.dim(2).size;
    fd.n23 = D.dim(1).size*fd.n3v;
    fd.size = D.size;
    return fd;
  }

  inline proc unflat(fd: FlatDom, q: int): (int, int, int) {
    const i = fd.lo1 + q/fd.n23;
    const r = q%fd.n23;
    return (i, fd.lo2 + r/fd.n3v, fd.lo3 + r%fd.n3v);
  }

  /* ---- two-hop slab copies ------------------------------------------
     device <-> local host buffer <-> (possibly remote) StencilDist
     array; the intermediate buffer keeps every leg of the transfer a
     plain local/DMA copy.  Only use these for slabs that are
     CONTIGUOUS on the device (x1 planes spanning the full padded
     x2/x3 extent, or whole arrays) — strided slabs go through the
     pack kernels below. */
  proc slabToHost(ref hA, const ref dA, dom) {
    if dom.size == 0 then return;
    var tmp: [dom] hA.eltType;
    tmp = dA[dom];
    hA[dom] = tmp;
  }

  proc slabToDev(ref dA, const ref hA, dom) {
    if dom.size == 0 then return;
    var tmp: [dom] hA.eltType;
    tmp = hA[dom];
    dA[dom] = tmp;
  }

  /* ---- packed slab copies (strided x2/x3 ghost planes) ------------- */

  inline proc linIdx(dom, i: int, j: int, k: int): int {
    const d2 = dom.dim(1), d3 = dom.dim(2);
    return ((i - dom.dim(0).low)*d2.size + (j - d2.low))*d3.size
           + (k - d3.low);
  }

  /* flat iteration (see FlatDom): the packed buffer index IS the flat
     loop index */
  proc packKernel(const ref dA, ref buf, dom) {
    const fd = mkFlat(dom);
    forall q in 0..#fd.size {
      const (i, j, k) = unflat(fd, q);
      buf[q] = dA[i, j, k];
    }
  }

  proc unpackKernel(ref dA, const ref buf, dom) {
    const fd = mkFlat(dom);
    forall q in 0..#fd.size {
      const (i, j, k) = unflat(fd, q);
      dA[i, j, k] = buf[q];
    }
  }

  /* V-slab staging through the block's pack buffers.  The kernels are
     launched with the block's fields referenced directly inside the
     `on b.dev` block (routing device arrays through host-side ref
     formals hands the kernel a host-frame pointer — illegal access). */
  proc pkVToHost(b: unmanaged GpuBlock, dom, param d3: bool) {
    if dom.size == 0 then return;
    on b.dev {
      if d3 then packKernel(b.V, b.pk3, dom);
             else packKernel(b.V, b.pk2, dom);
    }
    var hb: [0..#dom.size] StateVec;
    if d3 then hb = b.pk3[0..#dom.size];
          else hb = b.pk2[0..#dom.size];
    forall (i, j, k) in dom do
      V[i, j, k] = hb[linIdx(dom, i, j, k)];
  }

  proc pkVToDev(b: unmanaged GpuBlock, dom, param d3: bool) {
    if dom.size == 0 then return;
    var hb: [0..#dom.size] StateVec;
    forall (i, j, k) in dom do
      hb[linIdx(dom, i, j, k)] = V[i, j, k];
    if d3 then b.pk3[0..#dom.size] = hb;
          else b.pk2[0..#dom.size] = hb;
    on b.dev {
      if d3 then unpackKernel(b.V, b.pk3, dom);
             else unpackKernel(b.V, b.pk2, dom);
    }
  }

  /* ---- setup --------------------------------------------------------
     called once, after the initial condition / restart state is
     complete on the host (V, U, solveMask, forcing tables, PHI). */
  proc gpuInit() {
    if !gpuEnabled then return;

    coforall loc in Locales do on loc {
      const nDev = here.gpus.size;
      if nDev == 0 then
        halt("GPU build, but no GPUs are visible on ", here.name,
             " — run on a GPU node (and with one co-locale per GPU ",
             "for multi-locale runs)");

      const myX1 = DInt.localSubdomain().dim(0);
      var rs: list(range);
      for r in chunks(myX1, min(nDev, myX1.size)) do rs.pushBack(r);

      ref myList = locBlocks[here.id];
      for i in 0..#rs.size {
        const g = here.gpus[i];
        const r = rs[i];
        if r.size < ng1 then
          halt("GPU block too thin in x1 (", r.size, " cells < ng1=", ng1,
               "): reduce the GPU/locale count or increase nx1");
        var b: unmanaged GpuBlock?;
        on g {
          const DA  = {r.low-ng1..r.high+ng1, 1-ng2..nx2+ng2,
                       1-ng3..nx3+ng3};
          const DI  = {r.low..r.high, 1..nx2, 1..nx3};
          const DF1 = {r.low..r.high+1, 1..nx2, 1..nx3};
          const DF2 = {r.low..r.high, 1..nx2+1, 1..nx3};
          const DF3 = {r.low..r.high, 1..nx2, 1..nx3+1};
          const DSG = if sgFourPiG > 0.0 then DA
                      else {1..0, 1..0, 1..0};
          const DPk2 = if act2
                       then {0..#(r.size*ng2*(nx3 + 2*ng3))}
                       else {0..#0};
          const DPk3 = if act3
                       then {0..#(r.size*nx2*ng3)}
                       else {0..#0};
          b = new unmanaged GpuBlock(dev=g, lo=r.low, hi=r.high,
                                     DA=DA, DI=DI, DF1=DF1, DF2=DF2,
                                     DF3=DF3, DSG=DSG,
                                     DPk2=DPk2, DPk3=DPk3);
        }
        myList.pushBack(b!);
      }
    }
    nGpuTotal = + reduce ([l in LocD] locBlocks[l].size);

    // full initial upload (state, mask, forcing tables, potential);
    // whole-array assignments are contiguous: one DMA each
    coforall loc in Locales do on loc {
      for b in locBlocks[here.id] {
        var tmp: [b.DA] StateVec;
        tmp = V[b.DA];  b.V = tmp;
        tmp = U[b.DA];  b.U = tmp;
        var tmpb: [b.DA] bool;
        tmpb = solveMask[b.DA];  b.mask = tmpb;
        if sgFourPiG > 0.0 {
          var tmpp: [b.DSG] real;
          tmpp = SelfGravity.PHI[b.DSG];  b.PHI = tmpp;
        }
        b.fkv = Forcing.kv;  b.fe1 = Forcing.e1v;  b.fe2 = Forcing.e2v;
        b.fph = Forcing.phs; b.fa1 = Forcing.am1;  b.fa2 = Forcing.am2;
      }
    }

    // self-test: make sure foralls on the blocks really launch kernels
    // (an ineligible loop silently falls back to the host and would be
    // orders of magnitude slower)
    startGpuDiagnostics();
    coforall loc in Locales do on loc do
      coforall b in locBlocks[here.id] do on b.dev {
        const fd = mkFlat(b.DI);
        forall q in 0..#fd.size {
          const (i, j, k) = unflat(fd, q);
          b.RHS[i, j, k] = b.U[i, j, k];
        }
      }
    stopGpuDiagnostics();
    const d = getGpuDiagnostics();
    var launches = 0;
    for x in d do launches += x.kernel_launch: int;
    if launches == 0 then
      halt("GPU self-test failed: no kernels were launched");
    writeln("  gpu        : ", nGpuTotal, " device(s), ",
            "self-test ok (", launches, " kernel launches)");
  }

  /* ---- device-side x2/x3 boundary conditions ------------------------
     x2/x3 are never split across blocks, so every standard BC type is
     block-local and can run as a device kernel — the dominant cost of
     the host-staged shell exchange (userdef and shear-periodic sides
     still fall back to the host path).  The formulas and the
     x1->x2->x3 ordering mirror Boundary.applySide exactly. */

  inline proc stdBCdir(dir: int): bool {   // dir 1 = x2, 2 = x3
    const cmin = bcCode(2*dir), cmax = bcCode(2*dir + 1);
    inline proc ok(c: int): bool do
      return c == BC_PERIODIC || c == BC_ZEROGRAD || c == BC_REFLECT
          || c == BC_AXIS || c == BC_INFLOW || c == BC_OUT_DIODE
          || c == BC_IN_DIODE;
    return ok(cmin) && ok(cmax);
  }

  /* fill one x2/x3 ghost slab of a block on its device.
     side: 2:x2min 3:x2max 4:x3min 5:x3max (as in Boundary). */
  proc devSideBC(b: unmanaged GpuBlock, side: int) {
    const code = bcCode(side);
    const nv = IVX1 + side/2;             // normal velocity slot
    const ivAxis = if geom == Geom.polar then IVX2 else IVX3;
    const iso = eosCode == EOS_ISO;
    const l = b.lo, h = b.hi;
    const f2 = 1-ng2..nx2+ng2;
    // slab spans the full padded extent of the earlier directions
    // (corners become valid in the x1 -> x2 -> x3 pass order)
    const dom = if side == 2 then {l-ng1..h+ng1, 1-ng2..0, 1-ng3..nx3+ng3}
           else if side == 3 then {l-ng1..h+ng1, nx2+1..nx2+ng2,
                                   1-ng3..nx3+ng3}
           else if side == 4 then {l-ng1..h+ng1, f2, 1-ng3..0}
           else                   {l-ng1..h+ng1, f2, nx3+1..nx3+ng3};
    on b.dev {
      const fd = mkFlat(dom);
      forall q in 0..#fd.size {
        const (i, j, k) = unflat(fd, q);
        var w: StateVec;
        if code == BC_PERIODIC {
          w = if side == 2 then b.V[i, j+nx2, k]
         else if side == 3 then b.V[i, j-nx2, k]
         else if side == 4 then b.V[i, j, k+nx3]
         else                   b.V[i, j, k-nx3];
        } else if code == BC_REFLECT || code == BC_AXIS {
          w = if side == 2 then b.V[i, 1-j, k]
         else if side == 3 then b.V[i, 2*nx2+1-j, k]
         else if side == 4 then b.V[i, j, 1-k]
         else                   b.V[i, j, 2*nx3+1-k];
          w(nv) = -w(nv);
          if code == BC_AXIS then w(ivAxis) = -w(ivAxis);
        } else if code == BC_INFLOW {
          w = mkPrim(inRho, inVx1, inVx2, inVx3, inPrs);
        } else {                            // zero-gradient family
          w = if side == 2 then b.V[i, 1, k]
         else if side == 3 then b.V[i, nx2, k]
         else if side == 4 then b.V[i, j, 1]
         else                   b.V[i, j, nx3];
          if code == BC_OUT_DIODE {
            if side % 2 == 0 then w(nv) = min(w(nv), 0.0);
                             else w(nv) = max(w(nv), 0.0);
          } else if code == BC_IN_DIODE {
            if side % 2 == 0 then w(nv) = max(w(nv), 0.0);
                             else w(nv) = min(w(nv), 0.0);
          }
        }
        if iso then w(IPRS) = w(IRHO)*cs2At(i, j, k);
        b.V[i, j, k] = w;
      }
    }
  }

  /* ---- per-stage ghost/halo refresh ---------------------------------
     1. pull each block's ng-deep interior edge shell to the host (the
        source cells for boundary conditions and halo exchange),
     2. run the problem hook + host boundary conditions (they also
        refresh the StencilDist fluff caches),
     3. push the padded ghost shell of each block back to its device. */
  proc gpuStageBCs(t: real) {
    if !gpuEnabled then return;
    var sw: stopwatch;
    if gpuTime then sw.start();

    /* standard x2/x3 BCs run on the devices (block-local; the big
       staging win) — userdef/shear sides fall back to the host path */
    const dev2 = act2 && stdBCdir(1);
    const dev3 = act3 && stdBCdir(2);

    coforall loc in Locales do on loc {
      coforall b in locBlocks[here.id] {
        const l = b.lo, h = b.hi;
        const f3 = 1-ng3..nx3+ng3;
        // x1 planes over the full padded extent are contiguous on the
        // device (stale device ghosts written to the host are
        // recomputed by applyBCs below); x2/x3 slabs are strided and
        // go through the pack kernels
        slabToHost(V, b.V, {l..l+ng1-1, 1-ng2..nx2+ng2, f3});
        slabToHost(V, b.V, {h-ng1+1..h, 1-ng2..nx2+ng2, f3});
        if act2 && !dev2 {
          pkVToHost(b, {l..h, 1..ng2, f3}, false);
          pkVToHost(b, {l..h, nx2-ng2+1..nx2, f3}, false);
        }
        if act3 && !dev3 {
          pkVToHost(b, {l..h, 1..nx2, 1..ng3}, true);
          pkVToHost(b, {l..h, 1..nx2, nx3-ng3+1..nx3}, true);
        }
      }
    }
    if gpuTime { sw.stop(); tDown += sw.elapsed(); sw.clear(); sw.start(); }

    Problems.problemInternalBC(t);
    Boundary.applyBCs(t, refreshU=false, doX2=!dev2, doX3=!dev3);
    if gpuTime { sw.stop(); tHostBC += sw.elapsed(); sw.clear(); sw.start(); }

    coforall loc in Locales do on loc {
      coforall b in locBlocks[here.id] {
        const l = b.lo, h = b.hi;
        const f3 = 1-ng3..nx3+ng3;
        // x1 slabs span the full padded x2/x3 extent (corners included;
        // with device-side x2/x3 BCs their corner values are stale
        // here and overwritten by devSideBC below, preserving the
        // x1 -> x2 -> x3 ordering of the host pass)
        slabToDev(b.V, V, {l-ng1..l-1, 1-ng2..nx2+ng2, f3});
        slabToDev(b.V, V, {h+1..h+ng1, 1-ng2..nx2+ng2, f3});
        if act2 && !dev2 {
          pkVToDev(b, {l..h, 1-ng2..0, f3}, false);
          pkVToDev(b, {l..h, nx2+1..nx2+ng2, f3}, false);
        }
        if act3 && !dev3 {
          pkVToDev(b, {l..h, 1..nx2, 1-ng3..0}, true);
          pkVToDev(b, {l..h, 1..nx2, nx3+1..nx3+ng3}, true);
        }
        if dev2 { devSideBC(b, 2); devSideBC(b, 3); }
        if dev3 { devSideBC(b, 4); devSideBC(b, 5); }
      }
    }
    if gpuTime { sw.stop(); tUp += sw.elapsed(); }
  }

  proc gpuPrintTimers() {
    if gpuTime then
      writeln("gpuTime: shell-down ", tDown, " s, host BCs ", tHostBC,
              " s, shell-up ", tUp, " s, RHS kernels ", tRHS,
              " s, update kernels ", tKern, " s, dt ", tDt, " s");
  }

  /* ---- host <-> device state movement for the host-side features ---- */

  /* interior primitives to the host (output, particles, self-gravity);
     the whole padded block is one contiguous DMA — cheaper than a
     strided interior copy.  The x2/x3 ghost slabs come along too:
     with device-side BCs they are only valid on the device, and the
     host particle advection interpolates into them. */
  proc gpuDownV() {
    if !gpuEnabled then return;
    coforall loc in Locales do on loc do
      coforall b in locBlocks[here.id] {
        var tmp: [b.DA] StateVec;
        tmp = b.V;
        V[b.DI] = tmp[b.DI];
        const l = b.lo, h = b.hi;
        const f3 = 1-ng3..nx3+ng3;
        if act2 {
          V[{l..h, 1-ng2..0, f3}] = tmp[{l..h, 1-ng2..0, f3}];
          V[{l..h, nx2+1..nx2+ng2, f3}] =
            tmp[{l..h, nx2+1..nx2+ng2, f3}];
        }
        if act3 {
          V[{l..h, 1..nx2, 1-ng3..0}] = tmp[{l..h, 1..nx2, 1-ng3..0}];
          V[{l..h, 1..nx2, nx3+1..nx3+ng3}] =
            tmp[{l..h, 1..nx2, nx3+1..nx3+ng3}];
        }
      }
  }

  /* interior conservatives to the host (restart, FARGO) */
  proc gpuDownU() {
    if !gpuEnabled then return;
    coforall loc in Locales do on loc do
      coforall b in locBlocks[here.id] {
        var tmp: [b.DA] StateVec;
        tmp = b.U;
        U[b.DI] = tmp[b.DI];
      }
  }

  /* full re-upload after a host-side modification (FARGO, restart) */
  proc gpuUpAll() {
    if !gpuEnabled then return;
    coforall loc in Locales do on loc do
      coforall b in locBlocks[here.id] {
        var tmp: [b.DA] StateVec;
        tmp = V[b.DA];  b.V = tmp;
        tmp = U[b.DA];  b.U = tmp;
      }
  }

  /* refreshed potential to the devices (after the host CG solve) */
  proc gpuUpPhi() {
    if !gpuEnabled then return;
    coforall loc in Locales do on loc do
      coforall b in locBlocks[here.id] {
        var tmp: [b.DSG] real;
        tmp = SelfGravity.PHI[b.DSG];
        b.PHI = tmp;
      }
  }

  /* refreshed OU amplitudes to the devices (once per step) */
  proc gpuUpForcing() {
    if !gpuEnabled then return;
    coforall loc in Locales do on loc do
      coforall b in locBlocks[here.id] {
        b.fa1 = Forcing.am1;
        b.fa2 = Forcing.am2;
      }
  }
}
