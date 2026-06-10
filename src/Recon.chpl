/* Recon.chpl — slope-limited piecewise reconstruction (on primitives). */
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

  /* Left/right primitive states at the face that separates cell "m"
     (at idx-1 along the sweep) from cell "c" (at idx).  Needs the two
     cells on each side of each of those: wmm, wm, wc, wp.            */
  inline proc faceStates(wmm: StateVec, wm: StateVec,
                         wc:  StateVec, wp: StateVec,
                         out wL: StateVec, out wR: StateVec) {
    if reconCode == RECON_CONST {
      wL = wm;
      wR = wc;
    } else {
      for param v in 0..NVAR-1 {
        wL(v) = wm(v) + 0.5*limitedSlope(wmm(v), wm(v), wc(v));
        wR(v) = wc(v) - 0.5*limitedSlope(wm(v),  wc(v), wp(v));
      }
      // guard against reconstruction-induced negativity
      if wL(IRHO) <= 0.0 || wL(IPRS) <= 0.0 then wL = wm;
      if wR(IRHO) <= 0.0 || wR(IPRS) <= 0.0 then wR = wc;
    }
  }
}
