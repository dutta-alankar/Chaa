/* Evolve.chpl — strong-stability-preserving Runge-Kutta time stepping. */
module Evolve {
  use Params, Grid, State, Eos, Hydro, Boundary, Problems;
  use Math;

  /* enforce floors (and the isothermal pressure constraint) and refresh
     primitives after a conservative update */
  proc applyFloorsAndPrims() {
    forall (i, j, k) in DInt {
      var w = cons2prim(U[i,j,k]);
      if eosCode == EOS_ISO then
        w(IPRS) = w(IRHO)*cs2At(i, j, k);
      V[i,j,k] = w;
      U[i,j,k] = prim2cons(w);  // keeps U consistent if a floor triggered
    }
  }

  /*  U <- cA*U0 + cB*(U + dt*L(U))  evaluated at stage time tStage  */
  proc stage(cA: real, cB: real, dt: real, tStage: real) {
    computeRHS(tStage);
    forall idx in DInt {
      if solveMask[idx] {
        for param v in 0..NTOT-1 do
          U[idx](v) = cA*U0[idx](v) + cB*(U[idx](v) + dt*RHS[idx](v));
      }
    }
    applyFloorsAndPrims();
    problemInternalBC(tStage);
    applyBCs(tStage);
  }

  proc advance(dt: real, t: real) {
    U0 = U;
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
        forall idx in DInt {
          if solveMask[idx] {
            for param v in 0..NTOT-1 do
              U[idx](v) = U0[idx](v) + dt*RHS[idx](v);
          }
        }
        applyFloorsAndPrims();
        problemInternalBC(t + dt);
        applyBCs(t + dt);
      }
      otherwise do halt("unknown integrator");
    }
    if coolLambda0 > 0.0 {
      applyCooling(dt);
      applyBCs(t + dt);
    }
  }

  /* optically thin cooling, edot = -rho^2 Lambda0 T^alpha, integrated
     exactly over the step (Townsend 2009 for a power law); operator
     split, hence unconditionally stable */
  proc applyCooling(dt: real) {
    forall (i, j, k) in DInt {
      if !solveMask[i,j,k] then continue;
      var w = V[i,j,k];
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
      V[i,j,k] = w;
      U[i,j,k] = prim2cons(w);
    }
  }
}
