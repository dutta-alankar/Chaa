/* KhIdefix.chpl — Kelvin-Helmholtz instability, Idefix test/HD/KHI:
 * isothermal shear flow (+-1) across a perturbed interface
 *   y_int = 0.5*(1 + 0.05*(sin(pi x / 2) + cos(4 pi x)))
 * on [0,4]x[0,1], periodic in x, outflow in y.
 * Run with --eos=isothermal --csIso=10.
 */
module KhIdefix {
  use Params, Grid, State;
  use Math;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const yInt = 0.5*(1.0 + 0.05*(sin(0.5*pi*x) + cos(4.0*pi*x)));
      V[i,j,k] = (1.0, if y > yInt then 1.0 else -1.0, 0.0, 0.0, 1.0);
    }
  }
}
