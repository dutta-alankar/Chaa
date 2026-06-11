/* Turbulence.chpl — driven turbulence in a periodic box (after the
 * AthenaPK turbulence problem): gas starts at rest and is stirred by
 * the solenoidal Ornstein-Uhlenbeck forcing module (--forceAmp > 0,
 * see src/Forcing.chpl).  Run isothermal (--eos=isothermal) for the
 * classic supersonic-turbulence setup, or ideal-gas for decaying heat.
 * The first tracer field is initialised as a half-box dye to visualise
 * and quantify mixing.
 */
module Turbulence {
  use Params, Grid, State, Eos;

  proc setup() {
    forall (i, j, k) in DInt {
      var w = mkPrim(inRho, 0.0, 0.0, 0.0, inPrs);
      if ISC < NTOT then
        w(ISC) = if x1c(i) < 0.5*(x1min + x1max) then 1.0 else 0.0;
      V[i,j,k] = w;
    }
  }
}
