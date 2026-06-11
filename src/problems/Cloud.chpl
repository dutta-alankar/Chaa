/* Cloud.chpl — wind-cloud interaction ("cloud-in-wind", after the
 * AthenaPK cloud problem): a dense spherical cloud of density contrast
 * cloudChi in pressure equilibrium with a uniform wind (inRho, inVx1,
 * inPrs).  The first tracer field dyes the cloud material, so mixing
 * can be followed quantitatively.
 *
 * Suggested BCs: --bcX1min=inflow, outflow-diode elsewhere.
 */
module Cloud {
  use Params, Grid, State, Eos;
  use ProblemUtils;

  proc setup() {
    forall (i, j, k) in DInt {
      const inside = distFromCentre(i,j,k) < cloudRad;
      var w = mkPrim(if inside then cloudChi*inRho else inRho,
                     if inside then 0.0 else inVx1,
                     0.0, 0.0, inPrs);
      if ISC < NTOT then
        w(ISC) = if inside then 1.0 else 0.0;
      V[i,j,k] = w;
    }
  }
}
