/* Sedov.chpl — Sedov-Taylor point explosion.
 *
 * Total energy sedovE0 is deposited as thermal pressure inside radius
 * sedovR0 around (cen1,cen2,cen3); the deposit volume is *measured* on
 * the actual mesh (inactive angular dimensions contribute their full
 * measure), so the similarity solution applies unchanged in every
 * geometry and dimensionality.
 */
module Sedov {
  use Params, Grid, State, Eos;
  use ProblemUtils;

  proc setup() {
    const vol = + reduce ([(i,j,k) in DInt]
                  (if distFromCentre(i,j,k) < sedovR0
                   then cellVol(i,j,k) else 0.0));
    if vol <= 0.0 then halt("sedov: no cells inside deposit radius");
    const pIn = (gam - 1.0)*sedovE0/vol;
    forall (i, j, k) in DInt {
      const p = if distFromCentre(i,j,k) < sedovR0 then pIn else sedovPrsAmb;
      V[i,j,k] = mkPrim(sedovRhoAmb, 0.0, 0.0, 0.0, p);
    }
  }
}
