/* TwoBlast.chpl — Woodward & Colella (1984) interacting blast waves.
 *
 * 1D Cartesian on [0,1] with reflecting walls; run to t = 0.038.
 */
module TwoBlast {
  use Params, Grid, State, Eos;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i);
      var p = 0.01;
      if x < 0.1 then p = 1000.0;
      else if x > 0.9 then p = 100.0;
      V[i,j,k] = mkPrim(1.0, 0.0, 0.0, 0.0, p);
    }
  }
}
