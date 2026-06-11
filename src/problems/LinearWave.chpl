/* LinearWave.chpl — travelling acoustic eigenmode (the AthenaPK/Athena
 * linear-wave convergence test, hydro sound wave):
 *
 *   rho = 1 + A sin(2 pi x/L),  vx = A cs sin(...),  p = p0 + A cs^2 sin(...)
 *
 * with p0 = 1/gam so cs = 1.  On a periodic box the exact solution
 * after t = L/cs is the initial condition; the L1 error measures the
 * accuracy of the reconstruction + time integrator combination.
 * Parameter: waveAmp.
 */
module LinearWave {
  use Params, Grid, State, Eos;
  use Math;

  proc setup() {
    const L = x1max - x1min;
    const p0 = 1.0/gam;                 // cs = sqrt(gam p0 / rho0) = 1
    forall (i, j, k) in DInt {
      const ph = sin(2.0*pi*(x1c(i) - x1min)/L);
      V[i,j,k] = mkPrim(1.0 + waveAmp*ph,
                        waveAmp*ph,
                        0.0, 0.0,
                        p0 + waveAmp*ph);
    }
  }
}
