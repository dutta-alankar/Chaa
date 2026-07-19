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
  use Params, Grid, State, Eos, Recon, Riemann, Forcing, SelfGravity;
  use Problems;
  use Math;
  import Fargo.wBg;
  import CompileParams.gpuEnabled;
  import Gpu;
  import Time;

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
    if gpuEnabled {
      /* one padded block per GPU: the sweeps run as device kernels and
         every face flux a block needs is computed locally (faces are
         recomputed on both sides of a block cut — no flux exchange) */
      var sw: Time.stopwatch;
      if Gpu.gpuTime then sw.start();
      coforall loc in Locales do on loc do
        coforall b in Gpu.locBlocks[here.id] do on b.dev {
          const fdI = Gpu.mkFlat(b.DI);
          forall q in 0..#fdI.size {
            const (i, j, k) = Gpu.unflat(fdI, q);
            var z: StateVec;
            b.RHS[i, j, k] = z;
          }
          /* the flux evaluation is split into a reconstruction and a
             solve kernel staged through WL/WR, and each is
             compile-time specialized to the active scheme: a fused
             all-scheme kernel needs 255 registers + ~16 KB of local
             spills per thread (~50x slowdown) */
          if act1 {
            reconKernel(b.V, b.WL, b.WR, b.DF1, 0);
            solveKernel(b.V, b.WL, b.WR, b.FLX, b.DF1, 0);
            divKernel(b.FLX, b.RHS, b.mask, b.DI, 0);
          }
          if act2 {
            reconKernel(b.V, b.WL, b.WR, b.DF2, 1);
            solveKernel(b.V, b.WL, b.WR, b.FLX, b.DF2, 1);
            divKernel(b.FLX, b.RHS, b.mask, b.DI, 1);
          }
          if act3 {
            reconKernel(b.V, b.WL, b.WR, b.DF3, 2);
            solveKernel(b.V, b.WL, b.WR, b.FLX, b.DF3, 2);
            divKernel(b.FLX, b.RHS, b.mask, b.DI, 2);
          }
          sourcesKernel(true, b.V, b.RHS, b.mask, b.PHI,
                        b.fkv, b.fe1, b.fe2, b.fph, b.fa1, b.fa2,
                        Forcing.nModes, b.DI, t);
        }
      if Gpu.gpuTime { sw.stop(); Gpu.tRHS += sw.elapsed(); }
    } else {
      forall idx in DInt {
        var z: StateVec;
        RHS[idx] = z;
      }
      if act1 then sweepDir(0);
      if act2 then sweepDir(1);
      if act3 then sweepDir(2);
      sourcesKernel(false, V, RHS, solveMask, SelfGravity.PHI,
                    Forcing.kv, Forcing.e1v, Forcing.e2v, Forcing.phs,
                    Forcing.am1, Forcing.am2, Forcing.nModes, DInt, t);
    }
  }

  proc sweepDir(dir: int) {
    const DF = if dir == 0 then DAll[1..nx1+1, 1..nx2, 1..nx3]
          else if dir == 1 then DAll[1..nx1, 1..nx2+1, 1..nx3]
          else                  DAll[1..nx1, 1..nx2, 1..nx3+1];

    fluxKernel(V, FLX, DF, dir);
    FLX.updateFluff();
    divKernel(FLX, RHS, solveMask, DInt, dir);
  }

  /* face fluxes (face idx = left face of cell idx along dir); written
     once, generically: the CPU path passes the distributed arrays, the
     GPU path a block's device arrays (then the forall is a kernel).
     The heavy scheme bodies (ppm/wenoz/limo3, the exact solver) are
     non-inline device functions so the kernel's parameter list stays
     within the hardware's 32 KB limit (see Riemann.chpl). */
  proc fluxKernel(const ref Vv, ref Ff, DF, dir: int) {
    const e = dirOffset(dir);
    const acts = (act1, act2, act3);
    // viscous-flux mode: 0 off, 1 Cartesian stress, 2 cyl/polar tau_Rphi
    const vmode = if mu <= 0.0 then 0
                  else if geom == Geom.cartesian then 1
                  else if (geom == Geom.cylindrical || geom == Geom.polar)
                       && dir == 0 then 2
                  else 0;
    const ivp = if geom == Geom.cylindrical then IVX3 else IVX2;
    const needDl = kappa > 0.0 || scDiff > 0.0;
    forall idx in DF {
      var wL, wR: StateVec;
      /* per-face metric factors, evaluated once at kernel level: the
         closed-form grid law expands to many captured constants, so
         the diffusion helpers below receive them as plain scalars
         (keeps the GPU kernel-parameter list within hardware limits) */
      const dc = dcAt(dir, idx(dir));
      const dl = if needDl then physDc(dir, idx(0), idx(1), idx(2))
                 else 0.0;
      if reconCode == RECON_PPM || reconCode == RECON_WENOZ then
        faceStates6(Vv[idx-e-e-e], Vv[idx-e-e], Vv[idx-e],
                    Vv[idx], Vv[idx+e], Vv[idx+e+e], wL, wR);
      else
        faceStates(Vv[idx-e-e], Vv[idx-e], Vv[idx], Vv[idx+e],
                   dc, wL, wR);
      var f: StateVec;
      if useFargo && dir == 1 {
        // FARGO: solve in the frame comoving with the background w(R),
        // then transform the flux back to total momentum/energy
        const wf = wBg(idx(0));
        wL(IVX2) -= wf;
        wR(IVX2) -= wf;
        f = riemannFlux(wL, wR, dir);
        f(IENG) += wf*f(IMX2) + 0.5*wf*wf*f(IRHO);
        f(IMX2) += wf*f(IRHO);
      } else {
        f = riemannFlux(wL, wR, dir);
      }
      if vmode == 1 {
        var span: 3*real;
        for param d in 0..2 do
          if d != dir && acts(d) then
            span(d) = 2.0*(centerCoord(d, idx(d)+1)
                           - centerCoord(d, idx(d)-1));
        addViscFluxCart(Vv, idx, dir, dc, span, acts, mu, f);
      } else if vmode == 2 {
        addViscFluxRPhi(Vv, idx, ivp, x1f(idx(0)), dc, mu, f);
      }
      if kappa > 0.0 then addConductionFlux(Vv, idx, dir, dl, kappa, f);
      if scDiff > 0.0 then addScalarDiffFlux(Vv, idx, dir, dl, scDiff, f);
      Ff[idx] = f;
    }
  }

  /* ---- two-pass GPU flux pipeline -----------------------------------
     reconstruction and Riemann solve as separate, compile-time-
     specialized kernels staged through the WL/WR face-state arrays;
     used only by the GPU path (the CPU path keeps the fused
     fluxKernel above, whose behaviour is unchanged). */
  proc reconKernel(const ref Vv, ref WLa, ref WRa, DF, dir: int) {
    select reconCode {
      when RECON_CONST  do reconP(RECON_CONST,  Vv, WLa, WRa, DF, dir);
      when RECON_LIMO3  do reconP(RECON_LIMO3,  Vv, WLa, WRa, DF, dir);
      when RECON_PPM    do reconP(RECON_PPM,    Vv, WLa, WRa, DF, dir);
      when RECON_WENOZ  do reconP(RECON_WENOZ,  Vv, WLa, WRa, DF, dir);
      otherwise         do reconP(RECON_LINEAR, Vv, WLa, WRa, DF, dir);
    }
  }

  proc reconP(param rc: int, const ref Vv, ref WLa, ref WRa, DF,
              dir: int) {
    const e = dirOffset(dir);
    const fd = Gpu.mkFlat(DF);
    forall q in 0..#fd.size {
      const idx = Gpu.unflat(fd, q);
      var wL, wR: StateVec;
      if rc == RECON_PPM || rc == RECON_WENOZ then
        faceStates6P(rc, Vv[idx-e-e-e], Vv[idx-e-e], Vv[idx-e],
                     Vv[idx], Vv[idx+e], Vv[idx+e+e], wL, wR);
      else
        faceStatesP(rc, Vv[idx-e-e], Vv[idx-e], Vv[idx], Vv[idx+e],
                    dcAt(dir, idx(dir)), wL, wR);
      WLa[idx] = wL;
      WRa[idx] = wR;
    }
  }

  proc solveKernel(const ref Vv, const ref WLa, const ref WRa, ref Ff,
                   DF, dir: int) {
    select rsCode {
      when RS_LLF   do solveP(RS_LLF,   Vv, WLa, WRa, Ff, DF, dir);
      when RS_HLL   do solveP(RS_HLL,   Vv, WLa, WRa, Ff, DF, dir);
      when RS_EXACT do solveP(RS_EXACT, Vv, WLa, WRa, Ff, DF, dir);
      otherwise     do solveP(RS_HLLC,  Vv, WLa, WRa, Ff, DF, dir);
    }
  }

  proc solveP(param rs: int, const ref Vv, const ref WLa,
              const ref WRa, ref Ff, DF, dir: int) {
    const acts = (act1, act2, act3);
    const vmode = if mu <= 0.0 then 0
                  else if geom == Geom.cartesian then 1
                  else if (geom == Geom.cylindrical || geom == Geom.polar)
                       && dir == 0 then 2
                  else 0;
    const ivp = if geom == Geom.cylindrical then IVX3 else IVX2;
    const needDl = kappa > 0.0 || scDiff > 0.0;
    const fd = Gpu.mkFlat(DF);
    forall q in 0..#fd.size {
      const idx = Gpu.unflat(fd, q);
      var wL = WLa[idx], wR = WRa[idx];
      const dc = dcAt(dir, idx(dir));
      const dl = if needDl then physDc(dir, idx(0), idx(1), idx(2))
                 else 0.0;
      var f: StateVec;
      if useFargo && dir == 1 {
        const wf = wBg(idx(0));
        wL(IVX2) -= wf;
        wR(IVX2) -= wf;
        f = riemannFluxP(rs, wL, wR, dir);
        f(IENG) += wf*f(IMX2) + 0.5*wf*wf*f(IRHO);
        f(IMX2) += wf*f(IRHO);
      } else {
        f = riemannFluxP(rs, wL, wR, dir);
      }
      if vmode == 1 {
        var span: 3*real;
        for param d in 0..2 do
          if d != dir && acts(d) then
            span(d) = 2.0*(centerCoord(d, idx(d)+1)
                           - centerCoord(d, idx(d)-1));
        addViscFluxCart(Vv, idx, dir, dc, span, acts, mu, f);
      } else if vmode == 2 {
        addViscFluxRPhi(Vv, idx, ivp, x1f(idx(0)), dc, mu, f);
      }
      if kappa > 0.0 then addConductionFlux(Vv, idx, dir, dl, kappa, f);
      if scDiff > 0.0 then addScalarDiffFlux(Vv, idx, dir, dl, scDiff, f);
      Ff[idx] = f;
    }
  }

  /* flux divergence with geometry-aware face areas and volumes;
     the body is shared between the CPU (domain forall) and GPU
     (flattened forall — see Gpu.FlatDom) drivers */
  inline proc divBody(const ref Ff, ref Rr, const ref Mm,
                      i: int, j: int, k: int, dir: int) {
    if !Mm[i, j, k] then return;
    if dir == 0 {
      const w = invV1(i);
      const aL = fA1(i), aR = fA1(i+1);
      for param v in 0..NTOT-1 do
        Rr[i,j,k](v) += (aL*Ff[i,j,k](v) - aR*Ff[i+1,j,k](v))*w;
    } else if dir == 1 {
      const w = invV2(j)*g2(i);
      const aL = fA2(j), aR = fA2(j+1);
      for param v in 0..NTOT-1 do
        Rr[i,j,k](v) += (aL*Ff[i,j,k](v) - aR*Ff[i,j+1,k](v))*w;
    } else {
      const w = invV3(k)*g3(i, j);
      for param v in 0..NTOT-1 do
        Rr[i,j,k](v) += (Ff[i,j,k](v) - Ff[i,j,k+1](v))*w;
    }
  }

  proc divKernel(const ref Ff, ref Rr, const ref Mm, D, dir: int) {
    if gpuEnabled {
      const fd = Gpu.mkFlat(D);
      forall q in 0..#fd.size {
        const (i, j, k) = Gpu.unflat(fd, q);
        divBody(Ff, Rr, Mm, i, j, k, dir);
      }
    } else {
      forall (i, j, k) in D do
        divBody(Ff, Rr, Mm, i, j, k, dir);
    }
  }

  /* geometric (curvilinear), gravitational, forcing and
     curvilinear-viscous source terms, evaluated at cell centres.
     onGpu=true compiles out the host-only problemBodyForce hook (a GPU
     run with the hook enabled is rejected at startup). */
  proc sourcesKernel(param onGpu: bool, const ref Vv, ref Rr,
                     const ref Mm, const ref PH,
                     const ref fkv, const ref fe1, const ref fe2,
                     const ref fph, const ref fa1, const ref fa2,
                     nm: int, D, t: real) {
    if onGpu {
      const fd = Gpu.mkFlat(D);
      forall q in 0..#fd.size {
        const (i, j, k) = Gpu.unflat(fd, q);
        srcBody(onGpu, Vv, Rr, Mm, PH, fkv, fe1, fe2, fph, fa1, fa2,
                nm, i, j, k, t);
      }
    } else {
      forall (i, j, k) in D do
        srcBody(onGpu, Vv, Rr, Mm, PH, fkv, fe1, fe2, fph, fa1, fa2,
                nm, i, j, k, t);
    }
  }

  inline proc srcBody(param onGpu: bool, const ref Vv, ref Rr,
                      const ref Mm, const ref PH,
                      const ref fkv, const ref fe1, const ref fe2,
                      const ref fph, const ref fa1, const ref fa2,
                      nm: int, i: int, j: int, k: int, t: real) {
    {
      if !Mm[i, j, k] then return;
      const w = Vv[i, j, k];
      var s: StateVec;

      // geometry if-chains, not `select`: Chapel 2.8 miscompiles
      // `select` on an enum const inside inlined procs in GPU kernels
      if geom == Geom.cylindrical {          // x1=R, x3=phi
        const iR = 1.0/rGeo(i);
        s(IMX1) += (w(IPRS) + w(IRHO)*w(IVX3)**2)*iR;
        s(IMX3) += -w(IRHO)*w(IVX1)*w(IVX3)*iR;
      } else if geom == Geom.polar {         // x1=R, x2=phi
        const iR = 1.0/rGeo(i);
        s(IMX1) += (w(IPRS) + w(IRHO)*w(IVX2)**2)*iR;
        s(IMX2) += -w(IRHO)*w(IVX1)*w(IVX2)*iR;
      } else if geom == Geom.spherical {     // x1=r, x2=theta, x3=phi
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
        if geom == Geom.spherical {
          const r = x1c(i);
          g(0) = -gravCentral*r/(r*r + gravEps*gravEps)**1.5;
        } else if geom == Geom.polar {       // x1=R, x3=z
          const R = x1c(i);
          const z = if act3 then x3c(k) else 0.0;
          const ir3 = 1.0/(R*R + z*z + gravEps*gravEps)**1.5;
          g(0) = -gravCentral*R*ir3;
          g(2) = -gravCentral*z*ir3;
        } else if geom == Geom.cylindrical { // x1=R, x2=z
          const R = x1c(i);
          const z = if act2 then x2c(j) else 0.0;
          const ir3 = 1.0/(R*R + z*z + gravEps*gravEps)**1.5;
          g(0) = -gravCentral*R*ir3;
          g(1) = -gravCentral*z*ir3;
        } else {                             // cartesian
          const p = physPos(i, j, k);
          const d = (p(0)-cen1, p(1)-cen2, p(2)-cen3);
          const ir3 = 1.0/(d(0)**2 + d(1)**2 + d(2)**2
                           + gravEps*gravEps)**1.5;
          for param c in 0..2 do g(c) = -gravCentral*d(c)*ir3;
        }
        for param c in 0..2 {
          s(IMX1 + c) += w(IRHO)*g(c);
          s(IENG)     += w(IRHO)*w(IVX1 + c)*g(c);
        }
      }

      /* shearing box: Coriolis + tidal sources (Cartesian frame
         rotating at omegaRot, shear parameter q = shearQ) */
      if omegaRot > 0.0 && geom == Geom.cartesian {
        const x = x1c(i);
        s(IMX1) += w(IRHO)*(2.0*omegaRot*w(IVX2)
                            + 2.0*shearQ*omegaRot*omegaRot*x);
        s(IMX2) += -2.0*omegaRot*w(IRHO)*w(IVX1);
        s(IENG) += w(IRHO)*w(IVX1)*2.0*shearQ*omegaRot*omegaRot*x;
      }

      /* self-gravity g = -grad(Phi) (Poisson solve, once per step) */
      if sgFourPiG > 0.0 {
        const g = sgAccel(PH, i, j, k);
        for param c in 0..2 {
          s(IMX1 + c) += w(IRHO)*g(c);
          s(IENG)     += w(IRHO)*w(IVX1 + c)*g(c);
        }
      }

      /* problem-defined body force (custom potentials/forces) */
      if !onGpu && problemHasBodyForce {
        const a = problemBodyForce(i, j, k, t);
        for param c in 0..2 {
          s(IMX1 + c) += w(IRHO)*a(c);
          s(IENG)     += w(IRHO)*w(IVX1 + c)*a(c);
        }
      }

      /* Ornstein-Uhlenbeck spectral forcing (turbulence driving) */
      if forceAmp > 0.0 {
        const a = forceAccel(fkv, fe1, fe2, fph, fa1, fa2, nm, i, j, k);
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
        const dvp = (Vv[i+1,j,k](ivp) - Vv[i-1,j,k](ivp))
                  / (x1c(i+1) - x1c(i-1));
        const tauC = mu*(dvp - w(ivp)/Rc);
        s(ivp + 0) += tauC/Rc;          // ivp slot == momentum slot
      }

      for param v in 0..NTOT-1 do
        Rr[i,j,k](v) += s(v);
    }
  }

  /* viscous flux at a face; subtracted from the advective flux.
     Cartesian: full Navier-Stokes stress tensor.
     Cylindrical/polar: tau_{R,phi} only (radial sweeps), which is the
     term that drives Taylor-Couette flow.                          */
  /* Not inline, and every metric/runtime constant arrives as an
     argument (see fluxKernel): self-contained device functions add
     nothing to the GPU kernel-parameter list. */
  proc addViscFluxCart(const ref Vv, idx: 3*int, dir: int, dc: real,
                       span: 3*real, acts: 3*bool, muv: real,
                       ref f: StateVec) {
    const e = dirOffset(dir);
    const cm = idx - e;

    // grad(c)(d) = d v_c / d x_d at the face
    var grad: 3*(3*real);
    for param c in 0..2 {
      grad(c)(dir) = (Vv[idx](IVX1+c) - Vv[cm](IVX1+c))/dc;
      for param d in 0..2 {
        if d != dir {
          if acts(d) {
            const et = dirOffset(d);
            grad(c)(d) = (Vv[cm+et](IVX1+c) - Vv[cm-et](IVX1+c)
                        + Vv[idx+et](IVX1+c) - Vv[idx-et](IVX1+c))
                        / span(d);
          } else {
            grad(c)(d) = 0.0;
          }
        }
      }
    }

    const divv = grad(0)(0) + grad(1)(1) + grad(2)(2);
    for param c in 0..2 {
      var tau = muv*(grad(c)(dir) + grad(dir)(c));
      if c == dir then tau -= (2.0/3.0)*muv*divv;
      const vf = 0.5*(Vv[cm](IVX1+c) + Vv[idx](IVX1+c));
      f(IMX1 + c) -= tau;
      f(IENG)     -= vf*tau;
    }
  }

  /* radial tau_{R,phi} flux for cylindrical/polar (see above) */
  proc addViscFluxRPhi(const ref Vv, idx: 3*int, ivp: int, Rf: real,
                       dc: real, muv: real, ref f: StateVec) {
    const (i, j, k) = idx;
    if abs(Rf) > 1.0e-14 {
      const vm = Vv[i-1,j,k](ivp), vc = Vv[i,j,k](ivp);
      const vf = 0.5*(vm + vc);
      const tau = muv*((vc - vm)/dc - vf/Rf);
      f(IMX1 + (ivp - IVX1)) -= tau;
      f(IENG) -= vf*tau;
    }
  }

  inline proc centerCoord(dir: int, q: int): real {
    if dir == 0 then return x1c(q);
    if dir == 1 then return x2c(q);
    return x3c(q);
  }

  /* explicit thermal conduction: energy flux -kappa dT/dl across the
     face, with T = p/rho and dl the physical centre-to-centre spacing */
  proc addConductionFlux(const ref Vv, idx: 3*int, dir: int, dl: real,
                         kap: real, ref f: StateVec) {
    const wm = Vv[idx - dirOffset(dir)], wc = Vv[idx];
    const Tm = wm(IPRS)/wm(IRHO), Tc = wc(IPRS)/wc(IRHO);
    f(IENG) -= kap*(Tc - Tm)/dl;
  }

  /* passive-scalar diffusion: F = -rho_face * D * ds/dl */
  proc addScalarDiffFlux(const ref Vv, idx: 3*int, dir: int, dl: real,
                         sD: real, ref f: StateVec) {
    const wm = Vv[idx - dirOffset(dir)], wc = Vv[idx];
    const rhoF = 0.5*(wm(IRHO) + wc(IRHO));
    for param sl in ISC..NTOT-1 do
      f(sl) -= rhoF*sD*(wc(sl) - wm(sl))/dl;
  }

  /* CFL time step: advective everywhere + explicit-diffusion limits */
  proc computeDt(): real {
    if gpuEnabled {
      var swd: Time.stopwatch;
      if Gpu.gpuTime then swd.start();
      var dtl: [Gpu.LocD] real = 1.0e30;
      coforall loc in Locales do on loc {
        var bdt: [0..#Gpu.locBlocks[here.id].size] real = 1.0e30;
        coforall (b, ib) in zip(Gpu.locBlocks[here.id], 0..) {
          var v = 1.0e30;
          on b.dev {
            const fd = Gpu.mkFlat(b.DI);
            forall q in 0..#fd.size {
              const idx = Gpu.unflat(fd, q);
              b.DTB[idx] = cellDt(b.V, b.mask, idx);
            }
            v = min reduce b.DTB;
          }
          bdt[ib] = v;
        }
        dtl[here.id] = min reduce bdt;
      }
      if Gpu.gpuTime { swd.stop(); Gpu.tDt += swd.elapsed(); }
      return min(min reduce dtl, dtMax);
    }
    const dtmin = min reduce ([idx in DInt] cellDt(V, solveMask, idx));
    return min(dtmin, dtMax);
  }

  inline proc cellDt(const ref Vv, const ref Mm, idx: 3*int): real {
    if !Mm[idx] then return 1.0e30;
    const (i, j, k) = idx;
    const w = Vv[idx];
    const cs = soundSpeed(w);
    var dta = 1.0e30, dlmin = 1.0e30;
    if act1 {
      dta = min(dta, dl1(i)/(abs(w(IVX1)) + cs));
      dlmin = min(dlmin, dl1(i));
    }
    if act2 {
      // FARGO: only the residual azimuthal velocity limits the step
      const v2 = if useFargo then w(IVX2) - wBg(i) else w(IVX2);
      dta = min(dta, dl2(i, j)/(abs(v2) + cs));
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
