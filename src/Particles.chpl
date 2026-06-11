/* Particles.chpl — Lagrangian tracer particles.
 *
 * Massless tracers advected with the gas: velocity is tri-linearly
 * interpolated from cell centres (any grid law, via the inverse of the
 * closed-form coordinate maps) and positions advance with RK2.
 * Periodic sides wrap; other sides clamp to the domain.
 * Enabled with --nParticles > 0; positions are written next to every
 * field dump as <problem>.particles.NNNN.txt.
 */
module Particles {
  use Params, Grid, State;
  use Math, Random, IO;

  const PD = {0..#max(nParticles, 1)};
  var ppos: [PD] 3*real;

  proc initParticles() {
    if nParticles <= 0 then return;
    var rng = new randomStream(real, seed = partSeed);
    for p in PD {
      ppos[p](0) = x1min + rng.next()*(x1max - x1min);
      ppos[p](1) = if act2 then x2min + rng.next()*(x2max - x2min)
                           else x2c(1);
      ppos[p](2) = if act3 then x3min + rng.next()*(x3max - x3min)
                           else x3c(1);
    }
  }

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

  inline proc centerOf(dir: int, q: int): real {
    if dir == 0 then return x1c(q);
    if dir == 1 then return x2c(q);
    return x3c(q);
  }

  /* base cell index and linear weight for interpolation along dir */
  inline proc locate(dir: int, x: real, n: int, active: bool): (int, real) {
    if !active then return (1, 0.0);
    var i = floor(faceIndexOf(dir, x)): int;
    if i < 1 then i = 1;
    if i > n then i = n;
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

  inline proc wrap(ref pos: 3*real) {
    const lims = ((x1min, x1max, act1, bcCode(0)),
                  (x2min, x2max, act2, bcCode(2)),
                  (x3min, x3max, act3, bcCode(4)));
    for param d in 0..2 {
      const (lo, hi, act, bc) = lims(d);
      if !act then continue;
      if bc == BC_PERIODIC {
        const L = hi - lo;
        pos(d) = lo + mod(pos(d) - lo, L);
      } else {
        if pos(d) < lo then pos(d) = lo;
        if pos(d) > hi then pos(d) = hi;
      }
    }
  }

  /* RK2 (midpoint) advection in the frozen velocity field of the step */
  proc advanceParticles(dt: real) {
    if nParticles <= 0 then return;
    forall p in PD {
      const v1 = velAt(ppos[p]);
      var mid = ppos[p];
      for param c in 0..2 do mid(c) += 0.5*dt*v1(c);
      wrap(mid);
      const v2 = velAt(mid);
      var fin = ppos[p];
      for param c in 0..2 do fin(c) += dt*v2(c);
      wrap(fin);
      ppos[p] = fin;
    }
  }

  proc writeParticles(path: string) throws {
    if nParticles <= 0 then return;
    var f = open(path, ioMode.cw);
    var w = f.writer(locking=false);
    w.writeln("# id  x  y  z   (physically mapped positions)");
    for p in 0..#nParticles do
      w.writef("%i %.12er %.12er %.12er\n",
               p, ppos[p](0), ppos[p](1), ppos[p](2));
    w.close();
  }
}
