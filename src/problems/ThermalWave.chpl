/* ThermalWave.chpl — decaying entropy mode for the thermal-conduction
 * module (Idefix test/HD/thermalDiffusion): uniform pressure with a
 * small sinusoidal density (hence temperature) perturbation,
 *
 *   rho = 1 - twAmp * sin(2 pi x / Lx),   p = 1,  v = 0.
 *
 * At constant pressure the mode decays at the analytic rate
 *   Gamma = kappa * (gam-1)/gam * k^2,   k = 2 pi / Lx.
 * Run with --kappa>0 on a periodic domain.
 */
module ThermalWave {
  use Params, Grid, State;
  use Math;

  proc setup() {
    const Lx = x1max - x1min;
    forall (i, j, k) in DInt {
      const x = x1c(i);
      V[i,j,k] = (1.0 - twAmp*sin(2.0*pi*(x - x1min)/Lx),
                  0.0, 0.0, 0.0, 1.0);
    }
  }
}
