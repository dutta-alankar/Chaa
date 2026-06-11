/* RayleighTaylor.chpl — heavy-over-light hydrostatic layer with a
 * single-mode velocity seed; needs constant gravity, e.g. --grav2=-0.1
 * on [-0.25,0.25]x[-0.75,0.75] (domain centred on the interface y=0).
 * Parameters: rtRho{Top,Bot}, rtPrs0, rtPert.
 */
module RayleighTaylor {
  use Params, Grid, State, Eos;
  use Math;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const rho = if y > 0.0 then rtRhoTop else rtRhoBot;
      const p = rtPrs0 + rho*grav2*y;    // piecewise hydrostatic
      // single-mode perturbation, zero at the walls
      const Lx = x1max - x1min, Ly = x2max - x2min;
      const vy = rtPert*0.25*(1.0 + cos(2.0*pi*x/Lx))
                          *(1.0 + cos(2.0*pi*y/Ly));
      V[i,j,k] = mkPrim(rho, 0.0, vy, 0.0, p);
    }
  }
}
