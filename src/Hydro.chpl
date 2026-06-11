/* Hydro.chpl — finite-volume right-hand side:
 *   dimension-by-dimension sweeps (reconstruct -> Riemann -> divergence),
 *   well-balanced geometric source terms, gravity (constant + central
 *   point mass), turbulence forcing, and explicit diffusion: viscosity
 *   (full stress tensor in Cartesian; tau_{R,phi} in cylindrical/polar),
 *   thermal conduction and passive-scalar diffusion.
 *
 * All spacings are evaluated locally, so the sweeps work unchanged on
 * uniform, logarithmic and geometrically stretched grids.
 */
module Hydro {
  use Params, Grid, State, Eos, Recon, Riemann, Forcing;
  use Math;

  inline proc dirOffset(dir: int): 3*int {
    if dir == 0 then return (1, 0, 0);
    if dir == 1 then return (0, 1, 0);
    return (0, 0, 1);
  }

  /* physical centre-to-centre distance across face idx along dir */
  inline proc physDc(dir: int, i: int, j: int, k: int): real {
    if dir == 0 then return x1c(i) - x1c(i-1);
    if dir == 1 {
      const d = x2c(j) - x2c(j-1);
      if geom == Geom.polar || geom == Geom.spherical then
        return x1c(i)*d;
      return d;
    }
    const d = x3c(k) - x3c(k-1);
    if geom == Geom.spherical then return x1c(i)*sin(x2c(j))*d;
    if geom == Geom.polar then return d;       // x3 = z
    return d;
  }

  proc computeRHS(t: real) {
    forall idx in DInt {
      var z: StateVec;
      RHS[idx] = z;
    }
    if act1 then sweepDir(0);
    if act2 then sweepDir(1);
    if act3 then sweepDir(2);
    addSources(t);
  }

  proc sweepDir(dir: int) {
    const e = dirOffset(dir);

    const DF = if dir == 0 then DAll[1..nx1+1, 1..nx2, 1..nx3]
          else if dir == 1 then DAll[1..nx1, 1..nx2+1, 1..nx3]
          else                  DAll[1..nx1, 1..nx2, 1..nx3+1];

    /* face fluxes (face idx = left face of cell idx along dir) */
    forall idx in DF {
      var wL, wR: StateVec;
      if reconCode == RECON_PPM || reconCode == RECON_WENOZ then
        faceStates6(V[idx-e-e-e], V[idx-e-e], V[idx-e],
                    V[idx], V[idx+e], V[idx+e+e], wL, wR);
      else
        faceStates(V[idx-e-e], V[idx-e], V[idx], V[idx+e],
                   dcAt(dir, idx(dir)), wL, wR);
      var f = riemannFlux(wL, wR, dir);
      if mu > 0.0 then addViscFlux(idx, dir, f);
      if kappa > 0.0 then addConductionFlux(idx, dir, f);
      if scDiff > 0.0 then addScalarDiffFlux(idx, dir, f);
      FLX[idx] = f;
    }
    FLX.updateFluff();

    /* flux divergence with geometry-aware face areas and volumes */
    forall (i, j, k) in DInt {
      if !solveMask[i, j, k] then continue;
      if dir == 0 {
        const w = invV1(i);
        const aL = fA1(i), aR = fA1(i+1);
        for param v in 0..NTOT-1 do
          RHS[i,j,k](v) += (aL*FLX[i,j,k](v) - aR*FLX[i+1,j,k](v))*w;
      } else if dir == 1 {
        const w = invV2(j)*g2(i);
        const aL = fA2(j), aR = fA2(j+1);
        for param v in 0..NTOT-1 do
          RHS[i,j,k](v) += (aL*FLX[i,j,k](v) - aR*FLX[i,j+1,k](v))*w;
      } else {
        const w = invV3(k)*g3(i, j);
        for param v in 0..NTOT-1 do
          RHS[i,j,k](v) += (FLX[i,j,k](v) - FLX[i,j,k+1](v))*w;
      }
    }
  }

  /* geometric (curvilinear), gravitational, forcing and
     curvilinear-viscous source terms, evaluated at cell centres */
  proc addSources(t: real) {
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

      /* Ornstein-Uhlenbeck spectral forcing (turbulence driving) */
      if forceAmp > 0.0 {
        const a = forceAccel(i, j, k);
        for param c in 0..2 {
          s(IMX1 + c) += w(IRHO)*a(c);
          s(IENG)     += w(IRHO)*w(IVX1 + c)*a(c);
        }
      }

      /* residual viscous source for the azimuthal momentum in
         cylindrical/polar:  (1/R^2) d(R^2 tau)/dR
                           = (1/R) d(R tau)/dR + tau/R              */
      if mu > 0.0 && (geom == Geom.cylindrical || geom == Geom.polar) {
        const ivp = if geom == Geom.cylindrical then IVX3 else IVX2;
        const Rc = x1c(i);
        const dvp = (V[i+1,j,k](ivp) - V[i-1,j,k](ivp))
                  / (x1c(i+1) - x1c(i-1));
        const tauC = mu*(dvp - w(ivp)/Rc);
        s(ivp + 0) += tauC/Rc;          // ivp slot == momentum slot
      }

      for param v in 0..NTOT-1 do
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

      // grad(c)(d) = d v_c / d x_d at the face
      var grad: 3*(3*real);
      for param c in 0..2 {
        grad(c)(dir) = (V[idx](IVX1+c) - V[cm](IVX1+c))
                     / dcAt(dir, idx(dir));
        for param d in 0..2 {
          if d != dir {
            const active = (d == 0 && act1) || (d == 1 && act2)
                        || (d == 2 && act3);
            if active {
              const et = dirOffset(d);
              const q = idx(d);
              const span = 2.0*(centerCoord(d, q+1) - centerCoord(d, q-1));
              grad(c)(d) = (V[cm+et](IVX1+c) - V[cm-et](IVX1+c)
                          + V[idx+et](IVX1+c) - V[idx-et](IVX1+c))/span;
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
        const tau = mu*((vc - vm)/(x1c(i) - x1c(i-1)) - vf/Rf);
        f(IMX1 + (ivp - IVX1)) -= tau;
        f(IENG) -= vf*tau;
      }
    }
  }

  inline proc centerCoord(dir: int, q: int): real {
    if dir == 0 then return x1c(q);
    if dir == 1 then return x2c(q);
    return x3c(q);
  }

  /* explicit thermal conduction: energy flux -kappa dT/dl across the
     face, with T = p/rho and dl the physical centre-to-centre spacing */
  inline proc addConductionFlux(idx: 3*int, dir: int, ref f: StateVec) {
    const (i, j, k) = idx;
    const wm = V[idx - dirOffset(dir)], wc = V[idx];
    const Tm = wm(IPRS)/wm(IRHO), Tc = wc(IPRS)/wc(IRHO);
    f(IENG) -= kappa*(Tc - Tm)/physDc(dir, i, j, k);
  }

  /* passive-scalar diffusion: F = -rho_face * D * ds/dl */
  inline proc addScalarDiffFlux(idx: 3*int, dir: int, ref f: StateVec) {
    const (i, j, k) = idx;
    const wm = V[idx - dirOffset(dir)], wc = V[idx];
    const dl = physDc(dir, i, j, k);
    const rhoF = 0.5*(wm(IRHO) + wc(IRHO));
    for param sl in ISC..NTOT-1 do
      f(sl) -= rhoF*scDiff*(wc(sl) - wm(sl))/dl;
  }

  /* CFL time step: advective everywhere + explicit-diffusion limits */
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
      dta = min(dta, dl1(i)/(abs(w(IVX1)) + cs));
      dlmin = min(dlmin, dl1(i));
    }
    if act2 {
      dta = min(dta, dl2(i, j)/(abs(w(IVX2)) + cs));
      dlmin = min(dlmin, dl2(i, j));
    }
    if act3 {
      dta = min(dta, dl3(i, j, k)/(abs(w(IVX3)) + cs));
      dlmin = min(dlmin, dl3(i, j, k));
    }
    var dt = cfl*dta;
    if mu > 0.0 then
      dt = min(dt, cflVisc*0.5*dlmin**2*w(IRHO)/(mu*ndim));
    if kappa > 0.0 {
      // thermal diffusivity chi = kappa*(gam-1)/rho
      const chi = kappa*(gam - 1.0)/w(IRHO);
      dt = min(dt, cflVisc*0.5*dlmin**2/(chi*ndim));
    }
    if scDiff > 0.0 then
      dt = min(dt, cflVisc*0.5*dlmin**2/(scDiff*ndim));
    return dt;
  }
}
