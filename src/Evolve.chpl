/* Evolve.chpl — strong-stability-preserving Runge-Kutta time stepping.
 *
 * The per-cell update loops are written once as generic kernel procs:
 * the CPU path passes the distributed arrays, the GPU path a block's
 * device arrays (the foralls then compile to device kernels).
 */
module Evolve {
  use Params, Grid, State, Eos, Hydro, Boundary, Problems;
  use Math;
  import CompileParams.gpuEnabled;
  import Gpu;
  import Time;

  /* enforce floors (and the isothermal pressure constraint) and refresh
     primitives after a conservative update; the body is shared between
     the CPU (domain forall) and GPU (flattened forall) drivers */
  inline proc floorsBody(ref Vv, ref Uu, i: int, j: int, k: int) {
    var w = cons2prim(Uu[i,j,k]);
    if eosCode == EOS_ISO then
      w(IPRS) = w(IRHO)*cs2At(i, j, k);
    Vv[i,j,k] = w;
    Uu[i,j,k] = prim2cons(w);  // keeps U consistent if a floor triggered
  }

  /* onDev=true: device-block call, flattened iteration (see
     Gpu.FlatDom); onDev=false: host/distributed arrays, domain forall
     (a GPU build still calls this on the host for restart/FARGO) */
  proc floorsKernel(ref Vv, ref Uu, D, param onDev: bool) {
    if onDev {
      const fd = Gpu.mkFlat(D);
      forall q in 0..#fd.size {
        const (i, j, k) = Gpu.unflat(fd, q);
        floorsBody(Vv, Uu, i, j, k);
      }
    } else {
      forall (i, j, k) in D do
        floorsBody(Vv, Uu, i, j, k);
    }
  }

  /* host-array version (also used before the GPU blocks exist) */
  proc applyFloorsAndPrims() {
    floorsKernel(V, U, DInt, false);
  }

  /*  fromU0=false:  U <- cA*U0 + cB*(U  + dt*L(U))
      fromU0=true :  U <- cA*U0 + cB*(U0 + dt*L(U))   (vl2 corrector)  */
  inline proc stageUpdateBody(ref Uu, const ref Uz, const ref Rr,
                              const ref Mm, idx: 3*int, cA: real,
                              cB: real, dt: real, param fromU0: bool) {
    if Mm[idx] {
      for param v in 0..NTOT-1 {
        const base = if fromU0 then Uz[idx](v) else Uu[idx](v);
        Uu[idx](v) = cA*Uz[idx](v) + cB*(base + dt*Rr[idx](v));
      }
    }
  }

  proc stageUpdateKernel(ref Uu, const ref Uz, const ref Rr,
                         const ref Mm, D, cA: real, cB: real, dt: real,
                         param fromU0: bool) {
    if gpuEnabled {      // only called on device blocks in a GPU build
      const fd = Gpu.mkFlat(D);
      forall q in 0..#fd.size do
        stageUpdateBody(Uu, Uz, Rr, Mm, Gpu.unflat(fd, q),
                        cA, cB, dt, fromU0);
    } else {
      forall idx in D do
        stageUpdateBody(Uu, Uz, Rr, Mm, idx, cA, cB, dt, fromU0);
    }
  }

  /*  U <- cA*U0 + cB*(U + dt*L(U))  evaluated at stage time tStage  */
  proc stage(cA: real, cB: real, dt: real, tStage: real) {
    computeRHS(tStage);
    if gpuEnabled {
      var sw: Time.stopwatch;
      if Gpu.gpuTime then sw.start();
      coforall loc in Locales do on loc do
        coforall b in Gpu.locBlocks[here.id] do on b.dev {
          stageUpdateKernel(b.U, b.U0, b.RHS, b.mask, b.DI,
                            cA, cB, dt, false);
          floorsKernel(b.V, b.U, b.DI, true);
        }
      if Gpu.gpuTime { sw.stop(); Gpu.tKern += sw.elapsed(); }
      Gpu.gpuStageBCs(tStage);
    } else {
      stageUpdateKernel(U, U0, RHS, solveMask, DInt, cA, cB, dt, false);
      applyFloorsAndPrims();
      problemInternalBC(tStage);
      applyBCs(tStage);
    }
  }

  proc advance(dt: real, t: real) {
    if gpuEnabled {
      coforall loc in Locales do on loc do
        coforall b in Gpu.locBlocks[here.id] do on b.dev do
          b.U0 = b.U;
    } else {
      U0 = U;
    }
    select tiCode {
      when TI_EULER {
        stage(0.0, 1.0, dt, t);
      }
      when TI_RK2 {
        stage(0.0, 1.0, dt, t);
        stage(0.5, 0.5, dt, t + dt);
      }
      when TI_RK3 {
        stage(0.0, 1.0, dt, t);
        stage(0.75, 0.25, dt, t + dt);
        stage(1.0/3.0, 2.0/3.0, dt, t + 0.5*dt);
      }
      when TI_VL2 {
        // predictor-corrector (van Leer / midpoint):
        //   U* = U0 + dt/2 L(U0);   U = U0 + dt L(U*)
        stage(0.0, 1.0, 0.5*dt, t);
        computeRHS(t + 0.5*dt);
        if gpuEnabled {
          coforall loc in Locales do on loc do
            coforall b in Gpu.locBlocks[here.id] do on b.dev {
              stageUpdateKernel(b.U, b.U0, b.RHS, b.mask, b.DI,
                                0.0, 1.0, dt, true);
              floorsKernel(b.V, b.U, b.DI, true);
            }
          Gpu.gpuStageBCs(t + dt);
        } else {
          stageUpdateKernel(U, U0, RHS, solveMask, DInt,
                            0.0, 1.0, dt, true);
          applyFloorsAndPrims();
          problemInternalBC(t + dt);
          applyBCs(t + dt);
        }
      }
      otherwise do halt("unknown integrator");
    }
    if coolLambda0 > 0.0 {
      if gpuEnabled {
        coforall loc in Locales do on loc do
          coforall b in Gpu.locBlocks[here.id] do on b.dev do
            coolKernel(b.V, b.U, b.mask, b.DI, dt);
        Gpu.gpuStageBCs(t + dt);
      } else {
        coolKernel(V, U, solveMask, DInt, dt);
        applyBCs(t + dt);
      }
    }
  }

  /* optically thin cooling, edot = -rho^2 Lambda0 T^alpha, integrated
     exactly over the step (Townsend 2009 for a power law); operator
     split, hence unconditionally stable */
  inline proc coolBody(ref Vv, ref Uu, const ref Mm,
                       i: int, j: int, k: int, dt: real) {
    if !Mm[i,j,k] then return;
    var w = Vv[i,j,k];
    const T0 = w(IPRS)/w(IRHO);
    const kfac = (gam - 1.0)*w(IRHO)*coolLambda0;
    var T1: real;
    if abs(coolAlpha - 1.0) < 1.0e-12 {
      T1 = T0*exp(-kfac*dt);
    } else {
      const ex = 1.0 - coolAlpha;
      const arg = T0**ex - ex*kfac*dt;
      T1 = if arg > 0.0 then arg**(1.0/ex) else coolTfloor;
    }
    w(IPRS) = w(IRHO)*max(T1, coolTfloor);
    Vv[i,j,k] = w;
    Uu[i,j,k] = prim2cons(w);
  }

  proc coolKernel(ref Vv, ref Uu, const ref Mm, D, dt: real) {
    if gpuEnabled {      // only called on device blocks in a GPU build
      const fd = Gpu.mkFlat(D);
      forall q in 0..#fd.size {
        const (i, j, k) = Gpu.unflat(fd, q);
        coolBody(Vv, Uu, Mm, i, j, k, dt);
      }
    } else {
      forall (i, j, k) in D do
        coolBody(Vv, Uu, Mm, i, j, k, dt);
    }
  }
}
