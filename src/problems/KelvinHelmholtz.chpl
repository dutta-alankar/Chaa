/* KelvinHelmholtz.chpl — shear instability on a doubly periodic unit
 * box; two shear layers at y = 0.25 and 0.75 with a localised seed.
 * Parameters: khRho{In,Out}, khV0, khPert, khPrs.
 */
module KelvinHelmholtz {
  use Params, Grid, State;
  use Math;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const inner = abs(y - 0.5) < 0.25;
      const rho = if inner then khRhoIn else khRhoOut;
      const vx  = if inner then khV0 else -khV0;
      const vy = khPert*sin(4.0*pi*x)
               * (exp(-(y-0.25)**2/0.005) + exp(-(y-0.75)**2/0.005));
      V[i,j,k] = (rho, vx, vy, 0.0, khPrs);
    }
  }
}
