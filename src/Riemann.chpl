/* Riemann.chpl — Riemann solvers: LLF (Rusanov), HLL, HLLC, and the
 * exact (Godunov) solver of Toro (2009, ch. 4). */
module Riemann {
  use Params, Eos;
  use Math;

  inline proc riemannFlux(wL: StateVec, wR: StateVec, dir: int): StateVec {
    var f: StateVec;
    select rsCode {
      when RS_LLF   do f = llf(wL, wR, dir);
      when RS_HLL   do f = hll(wL, wR, dir);
      when RS_HLLC  do f = hllc(wL, wR, dir);
      when RS_EXACT do f = exact(wL, wR, dir);
      otherwise     do f = hllc(wL, wR, dir);
    }
    // passive tracers: upwind concentration carried by the mass flux
    for param s in ISC..NTOT-1 do
      f(s) = (if f(IRHO) >= 0.0 then wL(s) else wR(s))*f(IRHO);
    return f;
  }

  /* compile-time-specialized variant (GPU flux path): each device
     kernel carries a single solver — the all-solver kernel spills
     several KB of locals per thread (see Recon.faceStatesP) */
  inline proc riemannFluxP(param rs: int, wL: StateVec, wR: StateVec,
                           dir: int): StateVec {
    var f: StateVec;
    if rs == RS_LLF        then f = llf(wL, wR, dir);
    else if rs == RS_HLL   then f = hll(wL, wR, dir);
    else if rs == RS_EXACT then f = exact(wL, wR, dir);
    else                        f = hllc(wL, wR, dir);
    for param s in ISC..NTOT-1 do
      f(s) = (if f(IRHO) >= 0.0 then wL(s) else wR(s))*f(IRHO);
    return f;
  }

  inline proc llf(wL: StateVec, wR: StateVec, dir: int): StateVec {
    const uL = prim2cons(wL), uR = prim2cons(wR);
    const fL = physFlux(wL, dir), fR = physFlux(wR, dir);
    const s = max(abs(wL(IVX1+dir)) + soundSpeed(wL),
                  abs(wR(IVX1+dir)) + soundSpeed(wR));
    var f: StateVec;
    for param v in 0..NVAR-1 do
      f(v) = 0.5*(fL(v) + fR(v)) - 0.5*s*(uR(v) - uL(v));
    return f;
  }

  /* Davis wave-speed estimates shared by HLL and HLLC */
  inline proc waveSpeeds(wL: StateVec, wR: StateVec, dir: int,
                         out sL: real, out sR: real) {
    const unL = wL(IVX1+dir), unR = wR(IVX1+dir);
    const cL = soundSpeed(wL), cR = soundSpeed(wR);
    sL = min(unL - cL, unR - cR);
    sR = max(unL + cL, unR + cR);
  }

  inline proc hll(wL: StateVec, wR: StateVec, dir: int): StateVec {
    var sL, sR: real;
    waveSpeeds(wL, wR, dir, sL, sR);
    const fL = physFlux(wL, dir), fR = physFlux(wR, dir);
    if sL >= 0.0 then return fL;
    if sR <= 0.0 then return fR;
    const uL = prim2cons(wL), uR = prim2cons(wR);
    var f: StateVec;
    for param v in 0..NVAR-1 do
      f(v) = (sR*fL(v) - sL*fR(v) + sL*sR*(uR(v) - uL(v)))/(sR - sL);
    return f;
  }

  inline proc hllc(wL: StateVec, wR: StateVec, dir: int): StateVec {
    var sL, sR: real;
    waveSpeeds(wL, wR, dir, sL, sR);
    const fL = physFlux(wL, dir), fR = physFlux(wR, dir);
    if sL >= 0.0 then return fL;
    if sR <= 0.0 then return fR;

    const uL = prim2cons(wL), uR = prim2cons(wR);
    const rhoL = wL(IRHO), rhoR = wR(IRHO);
    const unL = wL(IVX1+dir), unR = wR(IVX1+dir);
    const pL = wL(IPRS), pR = wR(IPRS);

    const sM = (pR - pL + rhoL*unL*(sL - unL) - rhoR*unR*(sR - unR))
             / (rhoL*(sL - unL) - rhoR*(sR - unR));
    const pStar = pL + rhoL*(sL - unL)*(sM - unL);

    var f: StateVec;
    if sM >= 0.0 {
      const fac = rhoL*(sL - unL)/(sL - sM);
      var uS: StateVec;
      uS(IRHO) = fac;
      uS(IMX1) = fac*wL(IVX1);
      uS(IMX2) = fac*wL(IVX2);
      uS(IMX3) = fac*wL(IVX3);
      uS(IMX1+dir) = fac*sM;
      uS(IENG) = fac*(uL(IENG)/rhoL + (sM - unL)*(sM + pL/(rhoL*(sL - unL))));
      for param v in 0..NVAR-1 do
        f(v) = fL(v) + sL*(uS(v) - uL(v));
    } else {
      const fac = rhoR*(sR - unR)/(sR - sM);
      var uS: StateVec;
      uS(IRHO) = fac;
      uS(IMX1) = fac*wR(IVX1);
      uS(IMX2) = fac*wR(IVX2);
      uS(IMX3) = fac*wR(IVX3);
      uS(IMX1+dir) = fac*sM;
      uS(IENG) = fac*(uR(IENG)/rhoR + (sM - unR)*(sM + pR/(rhoR*(sR - unR))));
      for param v in 0..NVAR-1 do
        f(v) = fR(v) + sR*(uS(v) - uR(v));
    }
    return f;
  }

  /* ---------------- exact (Godunov) Riemann solver -------------------
   *
   * The helpers below are deliberately NOT inline and take every
   * runtime constant (gam, prsFloor, the EOS mode) as an explicit
   * argument: Chapel outlines one GPU-kernel parameter per inlined
   * reference to an outer constant, and inlining the whole iterative
   * solver into the flux kernel overflows the 32 KB kernel-parameter
   * limit.  As self-contained functions they compile to plain device
   * calls and contribute nothing to the kernel's parameter list.   */

  private proc soundSpeedG(w: StateVec, g: real, iso: bool): real {
    if iso then return sqrt(w(IPRS)/w(IRHO));
    return sqrt(g*w(IPRS)/w(IRHO));
  }

  /* Toro's pressure function f_K(p) and its derivative for state K */
  private proc pressFn(p: real, rho: real, pk: real, ck: real, g: real,
                       out f: real, out df: real) {
    if p > pk {                                   // shock
      const A = 2.0/((g + 1.0)*rho);
      const B = (g - 1.0)/(g + 1.0)*pk;
      const sq = sqrt(A/(p + B));
      f  = (p - pk)*sq;
      df = sq*(1.0 - 0.5*(p - pk)/(B + p));
    } else {                                      // rarefaction
      f  = 2.0*ck/(g - 1.0)*((p/pk)**((g - 1.0)/(2.0*g)) - 1.0);
      df = (p/pk)**(-(g + 1.0)/(2.0*g))/(rho*ck);
    }
  }

  /* star-region pressure and velocity (Newton-Raphson) */
  private proc starState(wL: StateVec, wR: StateVec, dir: int,
                         g: real, pf: real, iso: bool,
                         out ps: real, out us: real) {
    const rhoL = wL(IRHO), rhoR = wR(IRHO);
    const uL = wL(IVX1+dir), uR = wR(IVX1+dir);
    const pL = wL(IPRS), pR = wR(IPRS);
    const cL = soundSpeedG(wL, g, iso), cR = soundSpeedG(wR, g, iso);

    var p = max(0.5*(pL + pR), pf);
    var fL, dfL, fR, dfR: real;
    for 1..60 {
      pressFn(p, rhoL, pL, cL, g, fL, dfL);
      pressFn(p, rhoR, pR, cR, g, fR, dfR);
      const dp = (fL + fR + (uR - uL))/(dfL + dfR);
      p = max(p - dp, 1.0e-14);
      if abs(dp) < 1.0e-12*p then break;
    }
    ps = p;
    pressFn(p, rhoL, pL, cL, g, fL, dfL);
    pressFn(p, rhoR, pR, cR, g, fR, dfR);
    us = 0.5*(uL + uR) + 0.5*(fR - fL);
  }

  /* sample the self-similar solution at xi = x/t = 0 (the cell face) and
     return (rho, un, p); transverse velocities are upwinded across the
     contact by the caller */
  private proc sampleFace(wL: StateVec, wR: StateVec, dir: int,
                          g: real, iso: bool, ps: real, us: real,
                          out rho: real, out un: real, out p: real) {
    const g1 = (g - 1.0)/(g + 1.0);
    const g2 = (g - 1.0)/(2.0*g);

    if us >= 0.0 {                                // face lies left of contact
      const rhoK = wL(IRHO), uK = wL(IVX1+dir), pK = wL(IPRS);
      const cK = soundSpeedG(wL, g, iso);
      if ps > pK {                                // left shock
        const sK = uK - cK*sqrt((g + 1.0)/(2.0*g)*ps/pK + g2);
        if sK >= 0.0 then { rho = rhoK; un = uK; p = pK; return; }
        rho = rhoK*((ps/pK + g1)/(g1*ps/pK + 1.0));
        un = us; p = ps; return;
      }
      // left rarefaction
      const sh = uK - cK;
      if sh >= 0.0 then { rho = rhoK; un = uK; p = pK; return; }
      const cs = cK*(ps/pK)**g2;
      if us - cs <= 0.0 {                         // left of/at the contact
        rho = rhoK*(ps/pK)**(1.0/g); un = us; p = ps; return;
      }
      // inside the fan
      un = 2.0/(g + 1.0)*(cK + (g - 1.0)/2.0*uK);
      const c = 2.0/(g + 1.0)*(cK + (g - 1.0)/2.0*uK);
      rho = rhoK*(c/cK)**(2.0/(g - 1.0));
      p = pK*(c/cK)**(2.0*g/(g - 1.0));
      return;
    }
    // face lies right of the contact (mirror)
    const rhoK = wR(IRHO), uK = wR(IVX1+dir), pK = wR(IPRS);
    const cK = soundSpeedG(wR, g, iso);
    if ps > pK {                                  // right shock
      const sK = uK + cK*sqrt((g + 1.0)/(2.0*g)*ps/pK + g2);
      if sK <= 0.0 then { rho = rhoK; un = uK; p = pK; return; }
      rho = rhoK*((ps/pK + g1)/(g1*ps/pK + 1.0));
      un = us; p = ps; return;
    }
    // right rarefaction
    const sh = uK + cK;
    if sh <= 0.0 then { rho = rhoK; un = uK; p = pK; return; }
    const cs = cK*(ps/pK)**g2;
    if us + cs >= 0.0 {
      rho = rhoK*(ps/pK)**(1.0/g); un = us; p = ps; return;
    }
    un = 2.0/(g + 1.0)*(-cK + (g - 1.0)/2.0*uK);
    const c = 2.0/(g + 1.0)*(cK - (g - 1.0)/2.0*uK);
    rho = rhoK*(c/cK)**(2.0/(g - 1.0));
    p = pK*(c/cK)**(2.0*g/(g - 1.0));
  }

  inline proc exact(wL: StateVec, wR: StateVec, dir: int): StateVec {
    var ps, us, rho, un, p: real;
    const iso = eosCode == EOS_ISO;
    starState(wL, wR, dir, gam, prsFloor, iso, ps, us);
    sampleFace(wL, wR, dir, gam, iso, ps, us, rho, un, p);

    // assemble the face primitive state; transverse components ride
    // with the flow across the contact
    var w = if us >= 0.0 then wL else wR;
    w(IRHO) = rho;
    w(IVX1+dir) = un;
    w(IPRS) = p;
    return physFlux(w, dir);
  }
}
