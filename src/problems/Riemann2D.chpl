/* Riemann2D.chpl — Lax & Liu four-quadrant 2D Riemann problem,
 * configuration 3, on [0,1]^2 split at (0.5, 0.5); run to t = 0.3.
 */
module Riemann2D {
  use Params, Grid, State;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      if x >= 0.5 && y >= 0.5 then
        V[i,j,k] = (1.5, 0.0, 0.0, 0.0, 1.5);
      else if x < 0.5 && y >= 0.5 then
        V[i,j,k] = (0.5323, 1.206, 0.0, 0.0, 0.3);
      else if x < 0.5 && y < 0.5 then
        V[i,j,k] = (0.138, 1.206, 1.206, 0.0, 0.029);
      else
        V[i,j,k] = (0.5323, 0.0, 1.206, 0.0, 0.3);
    }
  }
}
