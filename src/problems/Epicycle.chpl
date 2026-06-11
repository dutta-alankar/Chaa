/* Epicycle.chpl — epicyclic oscillation in the shearing box: the
 * equilibrium shear flow v_y = -q Omega x plus a uniform radial kick
 * v_x = waveAmp.  The kick oscillates at the epicyclic frequency
 *   kappa = sqrt(2 (2 - q)) Omega        (kappa = Omega for q = 3/2),
 * giving <v_x>(t) = waveAmp cos(kappa t) — a clean quantitative test
 * of the Coriolis/tidal source terms and shear-periodic boundaries.
 *
 * Run with --omegaRot>0, --bcX1min=shear-periodic,
 * --bcX1max=shear-periodic, periodic x2 (works with --fargo=on too).
 */
module Epicycle {
  use Params, Grid, State, Eos;

  proc setup() {
    forall (i, j, k) in DInt {
      V[i,j,k] = mkPrim(inRho,
                        waveAmp,
                        -shearQ*omegaRot*x1c(i),
                        0.0, inPrs);
    }
  }
}
