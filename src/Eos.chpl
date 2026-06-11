/* Eos.chpl — ideal-gas (gamma-law) equation of state and
 * primitive <-> conservative conversions with positivity floors.
 */
module Eos {
  use Params;
  use Math;

  /* build a primitive state from the five hydro fields
     (tracer slots default to zero) */
  inline proc mkPrim(rho: real, v1: real, v2: real, v3: real,
                     p: real): StateVec {
    var w: StateVec;
    w(IRHO) = rho;
    w(IVX1) = v1;
    w(IVX2) = v2;
    w(IVX3) = v3;
    w(IPRS) = p;
    return w;
  }

  inline proc prim2cons(w: StateVec): StateVec {
    var u: StateVec;
    u(IRHO) = w(IRHO);
    u(IMX1) = w(IRHO)*w(IVX1);
    u(IMX2) = w(IRHO)*w(IVX2);
    u(IMX3) = w(IRHO)*w(IVX3);
    u(IENG) = w(IPRS)/(gam - 1.0)
            + 0.5*w(IRHO)*(w(IVX1)**2 + w(IVX2)**2 + w(IVX3)**2);
    for param s in ISC..NTOT-1 do      // tracers: rho * concentration
      u(s) = w(IRHO)*w(s);
    return u;
  }

  inline proc cons2prim(u: StateVec): StateVec {
    var w: StateVec;
    const rho = max(u(IRHO), rhoFloor);
    w(IRHO) = rho;
    w(IVX1) = u(IMX1)/rho;
    w(IVX2) = u(IMX2)/rho;
    w(IVX3) = u(IMX3)/rho;
    const ek = 0.5*rho*(w(IVX1)**2 + w(IVX2)**2 + w(IVX3)**2);
    w(IPRS) = max((gam - 1.0)*(u(IENG) - ek), prsFloor);
    for param s in ISC..NTOT-1 do
      w(s) = u(s)/rho;
    return w;
  }

  /* for the isothermal EOS the state's pressure is kept at rho*cs^2,
     so cs = sqrt(p/rho) without the adiabatic-index factor */
  inline proc soundSpeed(w: StateVec): real {
    if eosCode == EOS_ISO then return sqrt(w(IPRS)/w(IRHO));
    return sqrt(gam*w(IPRS)/w(IRHO));
  }

  /* physical (advective + pressure) flux of the Euler equations along
     coordinate direction dir (0,1,2), from a primitive state */
  inline proc physFlux(w: StateVec, dir: int): StateVec {
    const un = w(IVX1 + dir);
    const E  = w(IPRS)/(gam - 1.0)
             + 0.5*w(IRHO)*(w(IVX1)**2 + w(IVX2)**2 + w(IVX3)**2);
    var f: StateVec;
    f(IRHO) = w(IRHO)*un;
    f(IMX1) = w(IRHO)*w(IVX1)*un;
    f(IMX2) = w(IRHO)*w(IVX2)*un;
    f(IMX3) = w(IRHO)*w(IVX3)*un;
    f(IMX1 + dir) += w(IPRS);
    f(IENG) = (E + w(IPRS))*un;
    for param s in ISC..NTOT-1 do      // tracers ride on the mass flux
      f(s) = w(s)*f(IRHO);
    return f;
  }
}
