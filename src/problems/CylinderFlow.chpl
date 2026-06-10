/* CylinderFlow.chpl — viscous flow past an immersed solid cylinder of
 * radius cylRad centred at (cen1,cen2): uniform inflow (inRho, inVx1,
 * inPrs) with --bcX1min=inflow and --mu>0.  The solid is enforced as an
 * internal boundary (solveMask) after every stage.
 */
module CylinderFlow {
  use Params, Grid, State, Eos;

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const inside = (x - cen1)**2 + (y - cen2)**2 < cylRad**2;
      solveMask[i,j,k] = !inside;
      V[i,j,k] = (inRho, if inside then 0.0 else inVx1, 0.0, 0.0, inPrs);
    }
    solveMask.updateFluff();
  }

  /* re-impose the solid state inside the cylinder after every stage */
  proc internalBC(t: real) {
    forall idx in DInt {
      if !solveMask[idx] {
        V[idx] = (inRho, 0.0, 0.0, 0.0, inPrs);
        U[idx] = prim2cons(V[idx]);
      }
    }
  }
}
