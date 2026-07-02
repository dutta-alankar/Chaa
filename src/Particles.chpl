/* Particles.chpl — distributed Lagrangian tracer particles.
 *
 * Massless tracers advected with the gas: velocity is tri-linearly
 * interpolated from cell centres (any grid law, via the inverse of the
 * closed-form coordinate maps) and positions advance with RK2.
 *
 * Fully distributed, owner-computes storage: every locale holds the
 * particles that live inside its block of the (StencilDist) grid in a
 * local bag, advances them with purely local field reads (interpolation
 * stencils at block edges are served by the fluff caches), and hands
 * particles that cross a block boundary to the receiving locale after
 * every step.  Positions are gathered only for output.
 *
 * Initial positions come from the problem hook `problemParticleInit`
 * (see Problems.chpl); the default is a uniform random scatter over the
 * whole domain (deterministic for a given --partSeed).
 * Periodic and shear-periodic sides wrap; other sides clamp.
 * Enabled with --nParticles > 0; positions are written next to every
 * field dump as <problem>.particles.NNNN.txt.
 */
module Particles {
  use Params, Grid, State, Problems;
  use Math, Random, IO, List, BlockDist;

  record TracerParticle {
    var id: int;
    var pos: 3*real;
  }

  /* one particle bag per locale, each stored on its locale */
  const LocDom = {0..#numLocales} dmapped new blockDist({0..#numLocales});
  var bag:    [LocDom] list(TracerParticle);
  var outbox: [LocDom] list(TracerParticle);   // emigrants of the step

  /* fractional face index of coordinate x along dir */
  inline proc faceIndexOf(dir: int, x: real): real {
    if dir == 0 then
      return lawIndex(gridCode(0), x1min, x1max, nx1, stretchX1,
                      stretchUniX1, x);
    if dir == 1 then
      return lawIndex(gridCode(1), x2min, x2max, nx2, stretchX2,
                      stretchUniX2, x);
    return lawIndex(gridCode(2), x3min, x3max, nx3, stretchX3,
                    stretchUniX3, x);
  }

  /* interior cell containing coordinate x along dir */
  inline proc cellOf(dir: int, x: real, n: int, active: bool): int {
    if !active then return 1;
    var i = floor(faceIndexOf(dir, x)): int;
    if i < 1 then i = 1;
    if i > n then i = n;
    return i;
  }

  /* locale owning the grid block that contains pos */
  inline proc ownerOf(pos: 3*real): int {
    const idx = (cellOf(0, pos(0), nx1, act1),
                 cellOf(1, pos(1), nx2, act2),
                 cellOf(2, pos(2), nx3, act3));
    return DInt.distribution.idxToLocale(idx).id;
  }

  proc initParticles() {
    if nParticles <= 0 then return;
    var pos0: [0..#nParticles] 3*real;
    if !problemParticleInit(pos0) {
      /* default: uniform random scatter over the whole domain */
      var rng = new randomStream(real, seed = partSeed);
      for p in 0..#nParticles {
        pos0[p](0) = x1min + rng.next()*(x1max - x1min);
        pos0[p](1) = if act2 then x2min + rng.next()*(x2max - x2min)
                             else x2c(1);
        pos0[p](2) = if act3 then x3min + rng.next()*(x3max - x3min)
                             else x3c(1);
      }
    }
    for p in 0..#nParticles do wrap(pos0[p], 0.0);
    /* scatter to the owning locales */
    coforall loc in Locales do on loc {
      for p in 0..#nParticles do
        if ownerOf(pos0[p]) == here.id then
          bag[here.id].pushBack(new TracerParticle(p, pos0[p]));
    }
  }

  inline proc centerOf(dir: int, q: int): real {
    if dir == 0 then return x1c(q);
    if dir == 1 then return x2c(q);
    return x3c(q);
  }

  /* base cell index and linear weight for interpolation along dir */
  inline proc locate(dir: int, x: real, n: int, active: bool): (int, real) {
    if !active then return (1, 0.0);
    var i = cellOf(dir, x, n, active);
    // interpolate between the bracketing cell centres
    if x < centerOf(dir, i) && i > 0 then i -= 1;
    const c0 = centerOf(dir, i), c1 = centerOf(dir, i+1);
    var t = (x - c0)/(c1 - c0);
    if t < 0.0 then t = 0.0;
    if t > 1.0 then t = 1.0;
    return (i, t);
  }

  proc velAt(pos: 3*real): 3*real {
    const (i, tx) = locate(0, pos(0), nx1, act1);
    const (j, ty) = locate(1, pos(1), nx2, act2);
    const (k, tz) = locate(2, pos(2), nx3, act3);
    var v: 3*real;
    for param c in 0..2 {
      const c00 = (1.0-tx)*V[i,j,k](IVX1+c)   + tx*V[i+1,j,k](IVX1+c);
      const c10 = if act2
        then (1.0-tx)*V[i,j+1,k](IVX1+c) + tx*V[i+1,j+1,k](IVX1+c)
        else c00;
      const c01 = if act3
        then (1.0-tx)*V[i,j,k+1](IVX1+c) + tx*V[i+1,j,k+1](IVX1+c)
        else c00;
      const c11 = if act2 && act3
        then (1.0-tx)*V[i,j+1,k+1](IVX1+c) + tx*V[i+1,j+1,k+1](IVX1+c)
        else c10;
      const cy0 = (1.0-ty)*c00 + ty*c10;
      const cy1 = (1.0-ty)*c01 + ty*c11;
      v(c) = (1.0-tz)*cy0 + tz*cy1;
    }
    return v;
  }

  inline proc wrap(ref pos: 3*real, t: real) {
    /* shear-periodic x1 sides: the image of a point leaving through
       x1max sits at (x-Lx, y - q*Omega*Lx*t), matching Boundary.chpl */
    if act1 && bcCode(0) == BC_SHEAR {
      const Lx = x1max - x1min, off = shearQ*omegaRot*Lx*t;
      if pos(0) >= x1max { pos(0) -= Lx; pos(1) -= off; }
      else if pos(0) < x1min { pos(0) += Lx; pos(1) += off; }
    }
    const lims = ((x1min, x1max, act1, bcCode(0)),
                  (x2min, x2max, act2, bcCode(2)),
                  (x3min, x3max, act3, bcCode(4)));
    for param d in 0..2 {
      const (lo, hi, act, bc) = lims(d);
      if !act then continue;
      if bc == BC_PERIODIC || (d == 1 && bcCode(0) == BC_SHEAR) {
        const L = hi - lo;
        pos(d) = lo + mod(pos(d) - lo, L);
      } else if !(d == 0 && bc == BC_SHEAR) {
        if pos(d) < lo then pos(d) = lo;
        if pos(d) > hi then pos(d) = hi;
      }
    }
  }

  /* RK2 (midpoint) advection in the frozen velocity field of the step.
     Each locale advances its own particles (local reads, block-edge
     stencils served by the fluff caches), then emigrants migrate. */
  proc advanceParticles(dt: real, t: real) {
    if nParticles <= 0 then return;
    coforall loc in Locales do on loc {
      var arr = bag[here.id].toArray();
      forall pt in arr {
        const v1 = velAt(pt.pos);
        var mid = pt.pos;
        for param c in 0..2 do mid(c) += 0.5*dt*v1(c);
        wrap(mid, t + 0.5*dt);
        const v2 = velAt(mid);
        for param c in 0..2 do pt.pos(c) += dt*v2(c);
        wrap(pt.pos, t + dt);
      }
      var keep, move: list(TracerParticle);
      for pt in arr do
        if ownerOf(pt.pos) == here.id then keep.pushBack(pt);
                                      else move.pushBack(pt);
      bag[here.id] = keep;
      outbox[here.id] = move;
    }
    if numLocales > 1 then
      coforall loc in Locales do on loc {
        for src in 0..#numLocales {
          if src == here.id then continue;
          for pt in outbox[src] do
            if ownerOf(pt.pos) == here.id then bag[here.id].pushBack(pt);
        }
      }
  }

  proc writeParticles(path: string) throws {
    if nParticles <= 0 then return;
    /* gather positions by particle id (order-independent of ownership) */
    var all: [0..#nParticles] 3*real;
    var nGot = 0;
    for l in 0..#numLocales {
      const arr = bag[l].toArray();
      for pt in arr do all[pt.id] = pt.pos;
      nGot += arr.size;
    }
    if nGot != nParticles then
      halt("particle bookkeeping lost particles: ", nGot, " of ",
           nParticles);
    var f = open(path, ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("# id  x  y  z   (physically mapped positions)");
    for p in 0..#nParticles do
      w.writef("%i %.12er %.12er %.12er\n",
               p, all[p](0), all[p](1), all[p](2));
    w.close();
  }
}
