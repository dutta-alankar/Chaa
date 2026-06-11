/* Sod.chpl — shock tube along x1.
 *
 * Works in any geometry: planar in Cartesian, a radial shock tube in
 * cylindrical/spherical (use --bcX1min=axis there).
 * Parameters: sodX0, sodRho{L,R}, sodVx{L,R}, sodPrs{L,R}.
 */
module Sod {
  use Params, Grid, State, Eos;

  proc setup() {
    forall (i, j, k) in DInt {
      var w = if x1c(i) < sodX0
              then mkPrim(sodRhoL, sodVxL, 0.0, 0.0, sodPrsL)
              else mkPrim(sodRhoR, sodVxR, 0.0, 0.0, sodPrsR);
      if ISC < NTOT then              // dye the left state: the tracer
        w(ISC) = if x1c(i) < sodX0 then 1.0 else 0.0;  // marks the contact
      V[i,j,k] = w;
    }
  }
}
