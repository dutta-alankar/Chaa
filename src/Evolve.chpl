/* Evolve.chpl — strong-stability-preserving Runge-Kutta time stepping. */
module Evolve {
  use Params, State, Eos, Hydro, Boundary, Problems;

  /* enforce floors and refresh primitives after a conservative update */
  proc applyFloorsAndPrims() {
    forall idx in DInt {
      const w = cons2prim(U[idx]);
      V[idx] = w;
      U[idx] = prim2cons(w);   // keeps U consistent if a floor triggered
    }
  }

  /*  U <- cA*U0 + cB*(U + dt*L(U))  evaluated at stage time tStage  */
  proc stage(cA: real, cB: real, dt: real, tStage: real) {
    computeRHS(tStage);
    forall idx in DInt {
      if solveMask[idx] {
        for param v in 0..NVAR-1 do
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
      otherwise do halt("unknown integrator");
    }
  }
}
