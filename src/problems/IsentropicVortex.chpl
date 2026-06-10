/* IsentropicVortex.chpl — smooth isentropic vortex advected diagonally
 * with unit velocity; after t = (domain length) the exact solution is
 * the initial condition, giving a clean accuracy/convergence test.
 * Parameters: vortexBeta, centre (cen1,cen2).
 */
module IsentropicVortex {
  use Params, Grid, State;
  use Math;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i) - cen1, y = x2c(j) - cen2;
      const r2 = x*x + y*y;
      const ex = exp(0.5*(1.0 - r2));
      const dT = -(gam - 1.0)*vortexBeta**2/(8.0*gam*pi*pi)*exp(1.0 - r2);
      const T = 1.0 + dT;
      const rho = T**(1.0/(gam - 1.0));
      V[i,j,k] = (rho,
                  1.0 - vortexBeta/(2.0*pi)*ex*y,
                  1.0 + vortexBeta/(2.0*pi)*ex*x,
                  0.0,
                  rho*T);
    }
  }
}
