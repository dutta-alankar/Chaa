/* Sod.chpl — shock tube along x1.
 *
 * Works in any geometry: planar in Cartesian, a radial shock tube in
 * cylindrical/spherical (use --bcX1min=axis there).
 * Parameters: sodX0, sodRho{L,R}, sodVx{L,R}, sodPrs{L,R}.
 */
module Sod {
  use Params, Grid, State;

  proc setup() {
    forall (i, j, k) in DInt {
      if x1c(i) < sodX0 then
        V[i,j,k] = (sodRhoL, sodVxL, 0.0, 0.0, sodPrsL);
      else
        V[i,j,k] = (sodRhoR, sodVxR, 0.0, 0.0, sodPrsR);
    }
  }
}
