/* Hydro.chpl — finite-volume right-hand side:
 *   dimension-by-dimension sweeps (reconstruct -> Riemann -> divergence),
 *   well-balanced geometric source terms, constant gravity, and explicit
 *   viscous fluxes (full stress tensor in Cartesian; the dominant
 *   tau_{R,phi} term in cylindrical/polar for Taylor-Couette-like flows).
 */
module Hydro {
  use Params, Grid, State, Eos, Recon, Riemann;
  use Math;

  inline proc dirOffset(dir: int): 3*int {
    if dir == 0 then return (1, 0, 0);
    if dir == 1 then return (0, 1, 0);
    return (0, 0, 1);
  }

  proc computeRHS(t: real) {
    forall idx in DInt {
      var z: StateVec;
      RHS[idx] = z;
    }
    if act1 then sweepDir(0);
    if act2 then sweepDir(1);
    if act3 then sweepDir(2);
    addSources();
  }

  proc sweepDir(dir: int) {
    const e = dirOffset(dir);

    const DF = if dir == 0 then DAll[1..nx1+1, 1..nx2, 1..nx3]
          else if dir == 1 then DAll[1..nx1, 1..nx2+1, 1..nx3]
          else                  DAll[1..nx1, 1..nx2, 1..nx3+1];

    const dxd = if dir == 0 then dx1 else if dir == 1 then dx2 else dx3;

    /* face fluxes (face idx = left face of cell idx along dir) */
    forall idx in DF {
      var wL, wR: StateVec;
      if reconCode == RECON_PPM then
        faceStatesPPM(V[idx-e-e-e], V[idx-e-e], V[idx-e],
                      V[idx], V[idx+e], V[idx+e+e], wL, wR);
      else
        faceStates(V[idx-e-e], V[idx-e], V[idx], V[idx+e], dxd, wL, wR);
      var f = riemannFlux(wL, wR, dir);
      if mu > 0.0 then addViscFlux(idx, dir, f);
      if kappa > 0.0 then addConductionFlux(idx, dir, f);
      FLX[idx] = f;
    }
    FLX.updateFluff();

    /* flux divergence with geometry-aware face areas and volumes */
    forall (i, j, k) in DInt {
      if !solveMask[i, j, k] then continue;
      if dir == 0 {
        const w = invV1(i);
        const aL = fA1(i), aR = fA1(i+1);
        for param v in 0..NVAR-1 do
          RHS[i,j,k](v) += (aL*FLX[i,j,k](v) - aR*FLX[i+1,j,k](v))*w;
      } else if dir == 1 {
        const w = invV2(j)*g2(i);
        const aL = fA2(j), aR = fA2(j+1);
        for param v in 0..NVAR-1 do
          RHS[i,j,k](v) += (aL*FLX[i,j,k](v) - aR*FLX[i,j+1,k](v))*w;
      } else {
        const w = invV3()*g3(i, j);
        for param v in 0..NVAR-1 do
          RHS[i,j,k](v) += (FLX[i,j,k](v) - FLX[i,j,k+1](v))*w;
      }
    }
  }

  /* geometric (curvilinear), gravitational and curvilinear-viscous
     source terms, evaluated at cell centres */
  proc addSources() {
    forall (i, j, k) in DInt {
      if !solveMask[i, j, k] then continue;
      const w = V[i, j, k];
      var s: StateVec;

      select geom {
        when Geom.cylindrical {              // x1=R, x3=phi
          const iR = 1.0/rGeo(i);
          s(IMX1) += (w(IPRS) + w(IRHO)*w(IVX3)**2)*iR;
          s(IMX3) += -w(IRHO)*w(IVX1)*w(IVX3)*iR;
        }
        when Geom.polar {                    // x1=R, x2=phi
          const iR = 1.0/rGeo(i);
          s(IMX1) += (w(IPRS) + w(IRHO)*w(IVX2)**2)*iR;
          s(IMX2) += -w(IRHO)*w(IVX1)*w(IVX2)*iR;
        }
        when Geom.spherical {                // x1=r, x2=theta, x3=phi
          const ir = 1.0/rGeo(i);
          // the cot(theta) pressure term compensates the sin(theta)
          // area factor of the theta flux operator; without a theta
          // dimension the run is equatorially symmetric (cot == 0)
          const ct = if act2 then cotGeo(j) else 0.0;
          s(IMX1) += (2.0*w(IPRS) + w(IRHO)*(w(IVX2)**2 + w(IVX3)**2))*ir;
          s(IMX2) += (ct*(w(IPRS) + w(IRHO)*w(IVX3)**2)
                      - w(IRHO)*w(IVX1)*w(IVX2))*ir;
          s(IMX3) += -w(IRHO)*w(IVX3)*(w(IVX1) + w(IVX2)*ct)*ir;
        }
        otherwise { }
      }

      if grav1 != 0.0 || grav2 != 0.0 || grav3 != 0.0 {
        s(IMX1) += w(IRHO)*grav1;
        s(IMX2) += w(IRHO)*grav2;
        s(IMX3) += w(IRHO)*grav3;
        s(IENG) += w(IRHO)*(w(IVX1)*grav1 + w(IVX2)*grav2 + w(IVX3)*grav3);
      }

      /* central point-mass gravity g = -GM r_vec / (r^2+eps^2)^(3/2),
         with the mass at the coordinate origin (Cartesian: at the
         configured centre cen1..cen3) */
      if gravCentral > 0.0 {
        var g: 3*real = (0.0, 0.0, 0.0);
        select geom {
          when Geom.spherical {
            const r = x1c(i);
            g(0) = -gravCentral*r/(r*r + gravEps*gravEps)**1.5;
          }
          when Geom.polar {                  // x1=R, x3=z
            const R = x1c(i);
            const z = if act3 then x3c(k) else 0.0;
            const ir3 = 1.0/(R*R + z*z + gravEps*gravEps)**1.5;
            g(0) = -gravCentral*R*ir3;
            g(2) = -gravCentral*z*ir3;
          }
          when Geom.cylindrical {            // x1=R, x2=z
            const R = x1c(i);
            const z = if act2 then x2c(j) else 0.0;
            const ir3 = 1.0/(R*R + z*z + gravEps*gravEps)**1.5;
            g(0) = -gravCentral*R*ir3;
            g(1) = -gravCentral*z*ir3;
          }
          when Geom.cartesian {
            const p = physPos(i, j, k);
            const d = (p(0)-cen1, p(1)-cen2, p(2)-cen3);
            const ir3 = 1.0/(d(0)**2 + d(1)**2 + d(2)**2
                             + gravEps*gravEps)**1.5;
            for param c in 0..2 do g(c) = -gravCentral*d(c)*ir3;
          }
        }
        for param c in 0..2 {
          s(IMX1 + c) += w(IRHO)*g(c);
          s(IENG)     += w(IRHO)*w(IVX1 + c)*g(c);
        }
      }

      /* residual viscous source for the azimuthal momentum in
         cylindrical/polar:  (1/R^2) d(R^2 tau)/dR
                           = (1/R) d(R tau)/dR + tau/R              */
      if mu > 0.0 && (geom == Geom.cylindrical || geom == Geom.polar) {
        const ivp = if geom == Geom.cylindrical then IVX3 else IVX2;
        const Rc = x1c(i);
        const dvp = (V[i+1,j,k](ivp) - V[i-1,j,k](ivp))/(2.0*dx1);
        const tauC = mu*(dvp - w(ivp)/Rc);
        s(ivp + 0) += tauC/Rc;          // ivp slot == momentum slot
      }

      for param v in 0..NVAR-1 do
        RHS[i,j,k](v) += s(v);
    }
  }

  /* viscous flux at a face; subtracted from the advective flux.
     Cartesian: full Navier-Stokes stress tensor.
     Cylindrical/polar: tau_{R,phi} only (radial sweeps), which is the
     term that drives Taylor-Couette flow.                          */
  inline proc addViscFlux(idx: 3*int, dir: int, ref f: StateVec) {
    if geom == Geom.cartesian {
      const e = dirOffset(dir);
      const cm = idx - e;
      const dxs = (dx1, dx2, dx3);

      // grad(c)(d) = d v_c / d x_d at the face
      var grad: 3*(3*real);
      for param c in 0..2 {
        grad(c)(dir) = (V[idx](IVX1+c) - V[cm](IVX1+c))/dxs(dir);
        for param d in 0..2 {
          if d != dir {
            const active = (d == 0 && act1) || (d == 1 && act2)
                        || (d == 2 && act3);
            if active {
              const et = dirOffset(d);
              grad(c)(d) = (V[cm+et](IVX1+c) - V[cm-et](IVX1+c)
                          + V[idx+et](IVX1+c) - V[idx-et](IVX1+c))
                         / (4.0*dxs(d));
            } else {
              grad(c)(d) = 0.0;
            }
          }
        }
      }

      const divv = grad(0)(0) + grad(1)(1) + grad(2)(2);
      for param c in 0..2 {
        var tau = mu*(grad(c)(dir) + grad(dir)(c));
        if c == dir then tau -= (2.0/3.0)*mu*divv;
        const vf = 0.5*(V[cm](IVX1+c) + V[idx](IVX1+c));
        f(IMX1 + c) -= tau;
        f(IENG)     -= vf*tau;
      }
    } else if (geom == Geom.cylindrical || geom == Geom.polar) && dir == 0 {
      const ivp = if geom == Geom.cylindrical then IVX3 else IVX2;
      const (i, j, k) = idx;
      const Rf = x1f(i);
      if abs(Rf) > 1.0e-14 {
        const vm = V[i-1,j,k](ivp), vc = V[i,j,k](ivp);
        const vf = 0.5*(vm + vc);
        const tau = mu*((vc - vm)/dx1 - vf/Rf);
        f(IMX1 + (ivp - IVX1)) -= tau;
        f(IENG) -= vf*tau;
      }
    }
  }

  /* explicit thermal conduction: energy flux -kappa dT/dl across the
     face, with T = p/rho and dl the physical spacing along the sweep */
  inline proc addConductionFlux(idx: 3*int, dir: int, ref f: StateVec) {
    const (i, j, k) = idx;
    const wm = V[idx - dirOffset(dir)], wc = V[idx];
    const Tm = wm(IPRS)/wm(IRHO), Tc = wc(IPRS)/wc(IRHO);
    const dl = if dir == 0 then dl1()
          else if dir == 1 then dl2(i)
          else                  dl3(i, j);
    f(IENG) -= kappa*(Tc - Tm)/dl;
  }

  /* CFL time step: advective everywhere + explicit-diffusion limit */
  proc computeDt(): real {
    const dtmin = min reduce ([idx in DInt] cellDt(idx));
    return min(dtmin, dtMax);
  }

  inline proc cellDt(idx: 3*int): real {
    if !solveMask[idx] then return 1.0e30;
    const (i, j, k) = idx;
    const w = V[idx];
    const cs = soundSpeed(w);
    var dta = 1.0e30, dlmin = 1.0e30;
    if act1 {
      dta = min(dta, dl1()/(abs(w(IVX1)) + cs));
      dlmin = min(dlmin, dl1());
    }
    if act2 {
      dta = min(dta, dl2(i)/(abs(w(IVX2)) + cs));
      dlmin = min(dlmin, dl2(i));
    }
    if act3 {
      dta = min(dta, dl3(i, j)/(abs(w(IVX3)) + cs));
      dlmin = min(dlmin, dl3(i, j));
    }
    var dt = cfl*dta;
    if mu > 0.0 then
      dt = min(dt, cflVisc*0.5*dlmin**2*w(IRHO)/(mu*ndim));
    if kappa > 0.0 {
      // thermal diffusivity chi = kappa*(gam-1)/rho
      const chi = kappa*(gam - 1.0)/w(IRHO);
      dt = min(dt, cflVisc*0.5*dlmin**2/(chi*ndim));
    }
    return dt;
  }
}
