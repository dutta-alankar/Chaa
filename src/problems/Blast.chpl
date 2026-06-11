/* Blast.chpl — over-pressured circular/spherical region (PLUTO-like
 * blast wave): density/pressure blastRho{In,Out}, blastP{in,out} inside
 * and outside radius blastR0 around (cen1,cen2,cen3).
 */
module Blast {
  use Params, Grid, State, Eos;
  use ProblemUtils;

  proc setup() {
    forall (i, j, k) in DInt {
      const inside = distFromCentre(i,j,k) < blastR0;
      V[i,j,k] = mkPrim(if inside then blastRhoIn else blastRhoOut,
                  0.0, 0.0, 0.0,
                  if inside then blastPin else blastPout);
    }
  }
}
