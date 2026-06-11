/* Fargo.chpl — FARGO orbital advection (Masset 2000; Mignone et al.
 * 2012 flavour).  The azimuthal velocity is split into a static
 * background w(R) plus a residual:
 *
 *   - the azimuthal sweep solves the Riemann problem in the frame
 *     comoving with w (states are shifted, the returned flux is
 *     transformed back to total momentum/energy — see Hydro.sweepDir);
 *   - the azimuthal CFL condition uses only the residual velocity;
 *   - after each full step every (R[,z]) row is rigidly shifted by
 *     w dt with a conservative, slope-limited remap (periodic).
 *
 * Backgrounds: Keplerian sqrt(GM/R) in polar geometry, the linear
 * shear -q Omega x in the Cartesian shearing box.  Enable with
 * --fargo=on (requires a periodic/shear-periodic, uniform x2).
 */
module Fargo {
  use Params, Grid, State, Eos;
  use Math;

  /* background orbital (linear) velocity at radial index i */
  inline proc wBg(i: int): real {
    if geom == Geom.cartesian then return -shearQ*omegaRot*x1c(i);
    return sqrt(gravCentral/x1c(i));        // polar, Keplerian
  }

  /* conservative remap of every azimuthal row by w(i)*dt */
  proc fargoShift(dt: real) {
    if !useFargo then return;
    forall (i, k) in {1..nx1, 1..nx3} {
      // shift in cells (x2 is uniform; angular for polar)
      const s = if geom == Geom.cartesian
                then wBg(i)*dt/dx2At(1)
                else (wBg(i)/x1c(i))*dt/dx2At(1);
      const n0 = floor(s): int;
      const f = s - n0;                     // f in [0,1)

      // gather the row (local on one locale per row in 1D/2D decomp)
      var q, qs: [0..#nx2] StateVec;
      for j in 1..nx2 do q[j-1] = U[i,j,k];

      // integer circular shift: content moves +x2 by n0 cells
      for j in 0..#nx2 do
        qs[j] = q[mod(j - n0, nx2)];

      // fractional advection by f (donor cell + minmod slope)
      for j in 1..nx2 {
        var w: StateVec;
        for param v in 0..NTOT-1 {
          const jm = mod(j-1 - 1, nx2), j0 = j-1, jp = mod(j-1 + 1, nx2);
          const jmm = mod(j-1 - 2, nx2);
          // flux through the right face of cell j0 and of cell jm
          const slR = minmod2(qs[j0](v) - qs[jm](v), qs[jp](v) - qs[j0](v));
          const slL = minmod2(qs[jm](v) - qs[jmm](v), qs[j0](v) - qs[jm](v));
          const fR = f*(qs[j0](v) + 0.5*(1.0 - f)*slR);
          const fL = f*(qs[jm](v) + 0.5*(1.0 - f)*slL);
          w(v) = qs[j0](v) - fR + fL;
        }
        U[i,j,k] = w;
      }
    }
  }

  inline proc minmod2(a: real, b: real): real {
    if a*b <= 0.0 then return 0.0;
    return if abs(a) < abs(b) then a else b;
  }
}
