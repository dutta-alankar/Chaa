/* Recon.chpl — face reconstruction of the primitive variables.
 *
 *   constant  first-order donor cell
 *   linear    piecewise-linear MUSCL (minmod / vanleer / mc limiters)
 *   limo3     third-order LimO3 limiter (Cada & Torrilhon 2009)
 *   ppm       piecewise-parabolic face values with the extremum-
 *             preserving limiter of Colella & Sekora 2008 / Peterson &
 *             Hammett 2013 (needs NG >= 3 ghost layers)
 *
 * limo3 and ppm follow the implementations in Idefix
 * (src/fluid/RiemannSolver/slopeLimiter.hpp, Lesur et al. 2023).
 */
module Recon {
  use Params;
  use Math;

  inline proc minmod2(a: real, b: real): real {
    if a*b <= 0.0 then return 0.0;
    return if abs(a) < abs(b) then a else b;
  }

  inline proc limitedSlope(qm: real, q0: real, qp: real): real {
    const dm = q0 - qm, dp = qp - q0;
    select limCode {
      when LIM_MINMOD do return minmod2(dm, dp);
      when LIM_VANLEER {
        if dm*dp <= 0.0 then return 0.0;
        return 2.0*dm*dp/(dm + dp);
      }
      when LIM_MC {
        if dm*dp <= 0.0 then return 0.0;
        const c = 0.5*(dm + dp);
        return min(abs(c), 2.0*abs(dm), 2.0*abs(dp))*sgn(c);
      }
      otherwise do return 0.0;
    }
    return 0.0;
  }

  /* LimO3 limiter value (multiplies dvp); Cada & Torrilhon 2009 with
     the radius-of-curvature switch, as implemented in Idefix.
     Not inline (pure function): inlining the big reconstruction
     bodies into the GPU flux kernel overflows the 32 KB
     kernel-parameter limit (see Riemann.chpl). */
  proc limO3Lim(dvp: real, dvm: real, dx: real): real {
    param rad = 0.1, eps = 1.0e-12;
    const th = dvm/(dvp + 1.0e-16);
    const q = (2.0 + th)/3.0;

    var a = min(1.5, 2.0*th);
    a = min(q, a);
    const b = max(-0.5*th, a);
    const c = min(q, b);
    const psi0 = max(0.0, c);

    var eta = rad*dx;
    eta = (dvm*dvm + dvp*dvp)/(eta*eta);
    if eta <= 1.0 - eps then return q;
    if eta >= 1.0 + eps then return psi0;
    const psi = (1.0 - (eta - 1.0)/eps)*q + (1.0 + (eta - 1.0)/eps)*psi0;
    return 0.5*psi;
  }

  /* face L/R states for constant, linear and limo3 reconstruction;
     the face separates cell "m" (idx-1 along the sweep) from cell "c" */
  inline proc faceStates(wmm: StateVec, wm: StateVec,
                         wc:  StateVec, wp: StateVec, dx: real,
                         out wL: StateVec, out wR: StateVec) {
    select reconCode {
      when RECON_CONST {
        wL = wm;
        wR = wc;
      }
      when RECON_LINEAR {
        for param v in 0..NTOT-1 {
          wL(v) = wm(v) + 0.5*limitedSlope(wmm(v), wm(v), wc(v));
          wR(v) = wc(v) - 0.5*limitedSlope(wm(v),  wc(v), wp(v));
        }
        if wL(IRHO) <= 0.0 || wL(IPRS) <= 0.0 then wL = wm;
        if wR(IRHO) <= 0.0 || wR(IPRS) <= 0.0 then wR = wc;
      }
      when RECON_LIMO3 {
        for param v in 0..NTOT-1 {
          var dvm = wm(v) - wmm(v);
          var dvp = wc(v) - wm(v);
          wL(v) = wm(v) + 0.5*dvp*limO3Lim(dvp, dvm, dx);
          // positivity fallback to minmod (as in Idefix)
          if (v == IRHO || v == IPRS) && wL(v) <= 0.0 then
            wL(v) = wm(v) + 0.5*minmod2(dvm, dvp);
          dvm = dvp;
          dvp = wp(v) - wc(v);
          wR(v) = wc(v) - 0.5*dvm*limO3Lim(dvm, dvp, dx);
          if (v == IRHO || v == IPRS) && wR(v) <= 0.0 then
            wR(v) = wc(v) - 0.5*minmod2(dvm, dvp);
        }
      }
      otherwise { wL = wm; wR = wc; }
    }
  }

  /* ----------------------------- PPM -------------------------------- */

  /* limit an interpolated face value vph between v0 and vp1
     (CD11 sect. 4.3.1 / CS08 eq. 18 with FS18 corrections) */
  inline proc limitPPMFace(vm1: real, v0: real, vp1: real, vp2: real,
                           ref vph: real) {
    if (vp1 - vph)*(vph - v0) < 0.0 {
      const deltaL = vm1 - 2.0*v0 + vp1;
      const deltaC = 3.0*(v0 - 2.0*vph + vp1);
      const deltaR = v0 - 2.0*vp1 + vp2;
      var delta = 0.0;
      if sgn(deltaL) == sgn(deltaC) && sgn(deltaR) == sgn(deltaC) {
        param C = 1.25;
        delta = C*min(abs(deltaL), abs(deltaR));
        delta = sgn(deltaC)*min(delta, abs(deltaC));
      }
      vph = 0.5*(v0 + vp1) - delta/6.0;
    }
  }

  /* parabolic face values of cell "0" with extremum-preserving
     limiting (PH13 3.26-3.32); vl = left face, vr = right face.
     Not inline (pure function — see limO3Lim). */
  proc ppmStates(vm2: real, vm1: real, v0: real,
                 vp1: real, vp2: real,
                 out vl: real, out vr: real) {
    param n = 2;

    vr = 7.0/12.0*(v0 + vp1) - 1.0/12.0*(vm1 + vp2);
    vl = 7.0/12.0*(vm1 + v0) - 1.0/12.0*(vm2 + vp1);

    limitPPMFace(vm2, vm1, v0, vp1, vl);
    limitPPMFace(vm1, v0, vp1, vp2, vr);

    const d2qf  = 6.0*(vl + vr - 2.0*v0);
    const d2qc0 = vm1 + vp1 - 2.0*v0;
    const d2qcp = v0 + vp2 - 2.0*vp1;
    const d2qcm = vm2 + v0 - 2.0*vm1;

    var d2q = 0.0;
    if sgn(d2qf) == sgn(d2qc0) && sgn(d2qf) == sgn(d2qcp) &&
       sgn(d2qf) == sgn(d2qcm) {
      param C = 1.25;
      d2q = min(abs(d2qc0), abs(d2qcp));
      d2q = C*min(abs(d2qcm), d2q);
      d2q = sgn(d2qf)*min(abs(d2qf), d2q);
    }

    const qmax = max(abs(vm1), abs(v0), abs(vp1));
    var rho = 0.0;
    if abs(d2qf) > 1.0e-12*qmax then rho = d2q/d2qf;

    if (vr - v0)*(v0 - vl) <= 0.0 || (vm1 - v0)*(v0 - vp1) <= 0.0 {
      if rho <= 1.0 - 1.0e-12 {
        vl = v0 - rho*(v0 - vl);
        vr = v0 + rho*(vr - v0);
      }
    } else {
      if abs(vr - v0) >= n*abs(v0 - vl) then vr = v0 + n*(v0 - vl);
      if abs(vl - v0) >= n*abs(v0 - vr) then vl = v0 + n*(v0 - vr);
    }
  }

  /* ----------------------------- WENO-Z ----------------------------- */

  /* fifth-order WENO-Z (Borges et al. 2008) right-face value of the
     centre cell of the 5-point stencil.
     Not inline (pure function — see limO3Lim). */
  proc wenozFace(vm2: real, vm1: real, v0: real,
                 vp1: real, vp2: real): real {
    param eps = 1.0e-40;
    const b0 = 13.0/12.0*(vm2 - 2.0*vm1 + v0)**2
             + 0.25*(vm2 - 4.0*vm1 + 3.0*v0)**2;
    const b1 = 13.0/12.0*(vm1 - 2.0*v0 + vp1)**2 + 0.25*(vm1 - vp1)**2;
    const b2 = 13.0/12.0*(v0 - 2.0*vp1 + vp2)**2
             + 0.25*(3.0*v0 - 4.0*vp1 + vp2)**2;
    const tau5 = abs(b0 - b2);
    const a0 = 0.1*(1.0 + tau5/(b0 + eps));
    const a1 = 0.6*(1.0 + tau5/(b1 + eps));
    const a2 = 0.3*(1.0 + tau5/(b2 + eps));
    const q0 = (2.0*vm2 - 7.0*vm1 + 11.0*v0)/6.0;
    const q1 = (-vm1 + 5.0*v0 + 2.0*vp1)/6.0;
    const q2 = (2.0*v0 + 5.0*vp1 - vp2)/6.0;
    return (a0*q0 + a1*q1 + a2*q2)/(a0 + a1 + a2);
  }

  /* face L/R states for the 5-cell schemes (ppm, wenoz); needs cells
     idx-3..idx+2 along the sweep (NG >= 3) */
  inline proc faceStates6(wm3: StateVec, wmm: StateVec, wm: StateVec,
                          wc: StateVec, wp: StateVec, wpp: StateVec,
                          out wL: StateVec, out wR: StateVec) {
    if reconCode == RECON_PPM {
      for param v in 0..NTOT-1 {
        var dum, vlm, vrm: real;
        ppmStates(wm3(v), wmm(v), wm(v), wc(v), wp(v), vlm, vrm);
        wL(v) = vrm;                       // right face of cell m
        ppmStates(wmm(v), wm(v), wc(v), wp(v), wpp(v), vlm, dum);
        wR(v) = vlm;                       // left face of cell c
      }
    } else {                               // wenoz
      for param v in 0..NTOT-1 {
        wL(v) = wenozFace(wm3(v), wmm(v), wm(v), wc(v), wp(v));
        // mirrored stencil for the right state
        wR(v) = wenozFace(wpp(v), wp(v), wc(v), wm(v), wmm(v));
      }
    }
    if wL(IRHO) <= 0.0 || wL(IPRS) <= 0.0 then wL = wm;
    if wR(IRHO) <= 0.0 || wR(IPRS) <= 0.0 then wR = wc;
  }

  /* ---- compile-time-specialized variants (GPU flux path) -----------
     The GPU pipeline dispatches the runtime scheme code to a param so
     each device kernel only carries one reconstruction (register
     pressure: the all-scheme kernel spills ~10 KB of locals per
     thread).  The CPU path keeps the runtime-select versions above. */
  inline proc faceStatesP(param rc: int, wmm: StateVec, wm: StateVec,
                          wc: StateVec, wp: StateVec, dx: real,
                          out wL: StateVec, out wR: StateVec) {
    if rc == RECON_LINEAR {
      for param v in 0..NTOT-1 {
        wL(v) = wm(v) + 0.5*limitedSlope(wmm(v), wm(v), wc(v));
        wR(v) = wc(v) - 0.5*limitedSlope(wm(v),  wc(v), wp(v));
      }
      if wL(IRHO) <= 0.0 || wL(IPRS) <= 0.0 then wL = wm;
      if wR(IRHO) <= 0.0 || wR(IPRS) <= 0.0 then wR = wc;
    } else if rc == RECON_LIMO3 {
      for param v in 0..NTOT-1 {
        var dvm = wm(v) - wmm(v);
        var dvp = wc(v) - wm(v);
        wL(v) = wm(v) + 0.5*dvp*limO3Lim(dvp, dvm, dx);
        if (v == IRHO || v == IPRS) && wL(v) <= 0.0 then
          wL(v) = wm(v) + 0.5*minmod2(dvm, dvp);
        dvm = dvp;
        dvp = wp(v) - wc(v);
        wR(v) = wc(v) - 0.5*dvm*limO3Lim(dvm, dvp, dx);
        if (v == IRHO || v == IPRS) && wR(v) <= 0.0 then
          wR(v) = wc(v) - 0.5*minmod2(dvm, dvp);
      }
    } else {                       // constant / donor cell
      wL = wm;
      wR = wc;
    }
  }

  inline proc faceStates6P(param rc: int, wm3: StateVec, wmm: StateVec,
                           wm: StateVec, wc: StateVec, wp: StateVec,
                           wpp: StateVec,
                           out wL: StateVec, out wR: StateVec) {
    if rc == RECON_PPM {
      for param v in 0..NTOT-1 {
        var dum, vlm, vrm: real;
        ppmStates(wm3(v), wmm(v), wm(v), wc(v), wp(v), vlm, vrm);
        wL(v) = vrm;                       // right face of cell m
        ppmStates(wmm(v), wm(v), wc(v), wp(v), wpp(v), vlm, dum);
        wR(v) = vlm;                       // left face of cell c
      }
    } else {                               // wenoz
      for param v in 0..NTOT-1 {
        wL(v) = wenozFace(wm3(v), wmm(v), wm(v), wc(v), wp(v));
        wR(v) = wenozFace(wpp(v), wp(v), wc(v), wm(v), wmm(v));
      }
    }
    if wL(IRHO) <= 0.0 || wL(IPRS) <= 0.0 then wL = wm;
    if wR(IRHO) <= 0.0 || wR(IPRS) <= 0.0 then wR = wc;
  }
}
