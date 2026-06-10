/* TaylorCouette.chpl — viscous azimuthal flow between rotating
 * cylinders (1D cylindrical R in [x1min,x1max], requires --mu>0 and
 * --bcX1min=userdef --bcX1max=userdef).  The steady state is the
 * analytic Couette profile v_phi = aR + b/R.
 * Parameters: tcOmega{In,Out} (wall angular velocities), inRho, inPrs.
 */
module TaylorCouette {
  use Params, Grid, State;

  proc setup() {
    const R1 = x1min, R2 = x1max;
    forall (i, j, k) in DInt {
      const R = x1c(i);
      const om = tcOmegaIn + (tcOmegaOut - tcOmegaIn)*(R - R1)/(R2 - R1);
      V[i,j,k] = (inRho, 0.0, 0.0, om*R, inPrs);
    }
  }

  /* no-slip rotating walls: v_phi fixed at the wall value, v_R
     reflected, density/pressure mirrored */
  proc userBC(side: int, t: real) {
    const ivp = if geom == Geom.cylindrical then IVX3 else IVX2;
    const vWall = if side == 0 then tcOmegaIn*x1min else tcOmegaOut*x1max;
    const Dg = if side == 0 then DAll[1-ng1..0, 1..nx2, 1..nx3]
                            else DAll[nx1+1..nx1+ng1, 1..nx2, 1..nx3];
    forall (i, j, k) in Dg {
      const src = if side == 0 then (1-i, j, k) else (2*nx1+1-i, j, k);
      var w = V[src];
      w(IVX1) = -w(IVX1);
      w(ivp)  = 2.0*vWall - w(ivp);
      V[i,j,k] = w;
    }
  }
}
