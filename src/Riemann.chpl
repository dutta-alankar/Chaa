/* Riemann.chpl — approximate Riemann solvers: LLF (Rusanov), HLL, HLLC. */
module Riemann {
  use Params, Eos;
  use Math;

  inline proc riemannFlux(wL: StateVec, wR: StateVec, dir: int): StateVec {
    select rsCode {
      when RS_LLF  do return llf(wL, wR, dir);
      when RS_HLL  do return hll(wL, wR, dir);
      when RS_HLLC do return hllc(wL, wR, dir);
      otherwise    do return hllc(wL, wR, dir);
    }
    return hllc(wL, wR, dir);
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
}
