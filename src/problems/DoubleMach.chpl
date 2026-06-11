/* DoubleMach.chpl — double Mach reflection of a Mach 10 shock
 * (Woodward & Colella 1984) on [0,4]x[0,1], shock at 60 degrees.
 *
 * Use --bcX1min=inflow with the post-shock state and
 * --bcX2min=userdef --bcX2max=userdef; run to t = 0.2.
 */
module DoubleMach {
  use Params, Grid, State, Eos;
  use Math;

  const post = mkPrim(8.0, 8.25*sin(pi/3.0), -8.25*cos(pi/3.0),
                      0.0, 116.5);
  const pre  = mkPrim(1.4, 0.0, 0.0, 0.0, 1.0);

  proc setup() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const xs = 1.0/6.0 + y/tan(pi/3.0);
      V[i,j,k] = if x < xs then post else pre;
    }
  }

  proc userBC(side: int, t: real) {
    if side == 2 {                       // bottom: post-shock for x<1/6,
      const Dg = DAll[1-ng1..nx1+ng1, 1-ng2..0, 1..1];
      forall (i, j, k) in Dg {           // reflecting wall beyond
        if x1c(i) < 1.0/6.0 {
          V[i,j,k] = post;
        } else {
          var w = V[i, 1-j, k];
          w(IVX2) = -w(IVX2);
          V[i,j,k] = w;
        }
      }
    } else if side == 3 {                // top: exact moving-shock state
      const Dg = DAll[1-ng1..nx1+ng1, nx2+1..nx2+ng2, 1..1];
      const xs = 1.0/6.0 + (1.0 + 20.0*t)/tan(pi/3.0);
      forall (i, j, k) in Dg do
        V[i,j,k] = if x1c(i) < xs then post else pre;
    }
  }
}
