/* DiskCavity.chpl — rotationally supported, locally isothermal disk
 * with a density cavity (Idefix test/HD/RWI-cavity profile):
 *
 *   rho(R)  = (1/R) * 0.5 * (1.01 + tanh((R - diskJumpR)/diskJumpW))
 *   cs(R)   = diskH0 / sqrt(R)            (--csIso=diskH0 --csSlope=-0.5)
 *   vphi(R) from exact centrifugal balance with the pressure gradient.
 *
 * Run in polar geometry with --gravCentral=1 --eos=isothermal.
 */
module DiskCavity {
  use Params, Grid, State, Eos;
  use Math;

  inline proc rhoProf(R: real): real do
    return 0.5*(1.01 + tanh((R - diskJumpR)/diskJumpW))/R;

  inline proc prsProf(R: real): real do
    return rhoProf(R)*(diskH0*diskH0/R);    // p = rho * cs^2(R)

  proc setup() {
    forall (i, j, k) in DInt {
      const R = x1c(i);
      const rho = rhoProf(R);
      // vphi^2 = GM/R + (R/rho) dp/dR  (radial equilibrium)
      const dR = 1.0e-6*R;
      const dpdR = (prsProf(R + dR) - prsProf(R - dR))/(2.0*dR);
      const vphi2 = gravCentral/R + R*dpdR/rho;
      const vphi = sqrt(max(vphi2, 0.0));
      V[i,j,k] = mkPrim(rho, 0.0, vphi, 0.0, rho*diskH0*diskH0/R);
    }
  }
}
