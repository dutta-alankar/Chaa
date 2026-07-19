/* SelfGravity.chpl — Poisson solver for the gas's own gravity:
 *
 *     lap(Phi) = sgFourPiG * rho        (sgFourPiG = 4 pi G)
 *
 * solved with a (Jacobi-preconditioned-free) conjugate-gradient
 * iteration on the same finite-volume metric operators as the flux
 * divergence, so it works in every geometry and on every grid law.
 * Fully periodic domains use the Jeans swindle (the volume-weighted
 * mean density is subtracted); non-periodic sides impose Phi = 0 on
 * the ghosts (crude isolated boundaries).
 *
 * The potential is solved once per step on the evolving density (warm
 * started from the previous solution) and enters the momentum/energy
 * sources as g = -grad(Phi).
 */
module SelfGravity {
  use Params, Grid, State;
  use Math;

  // allocate the work arrays only when enabled
  const SGD = if sgFourPiG > 0.0 then DAll else DAll[1..0, 1..1, 1..1];
  var PHI: [SGD] real;          // gravitational potential (kept warm)
  var BB, RR, PP, AP: [SGD] real;

  const SGI = SGD[1..nx1, 1..nx2, 1..nx3];

  const allPeriodic =
    (bcCode(0) == BC_PERIODIC || bcCode(0) == BC_SHEAR) &&
    (!act2 || bcCode(2) == BC_PERIODIC) &&
    (!act3 || bcCode(4) == BC_PERIODIC);

  /* ghost fill for a scalar field: periodic wrap, else Dirichlet 0 */
  proc fieldBC(ref A: [SGD] real) {
    if act1 {
      const per = bcCode(0) == BC_PERIODIC || bcCode(0) == BC_SHEAR;
      forall (i, j, k) in SGD[1-ng1..0, 1..nx2, 1..nx3] do
        A[i,j,k] = if per then A[i+nx1, j, k] else 0.0;
      forall (i, j, k) in SGD[nx1+1..nx1+ng1, 1..nx2, 1..nx3] do
        A[i,j,k] = if per then A[i-nx1, j, k] else 0.0;
    }
    if act2 {
      const per = bcCode(2) == BC_PERIODIC;
      forall (i, j, k) in SGD[1-ng1..nx1+ng1, 1-ng2..0, 1..nx3] do
        A[i,j,k] = if per then A[i, j+nx2, k] else 0.0;
      forall (i, j, k) in SGD[1-ng1..nx1+ng1, nx2+1..nx2+ng2, 1..nx3] do
        A[i,j,k] = if per then A[i, j-nx2, k] else 0.0;
    }
    if act3 {
      const per = bcCode(4) == BC_PERIODIC;
      forall (i, j, k) in SGD[1-ng1..nx1+ng1, 1-ng2..nx2+ng2, 1-ng3..0] do
        A[i,j,k] = if per then A[i, j, k+nx3] else 0.0;
      forall (i, j, k) in SGD[1-ng1..nx1+ng1, 1-ng2..nx2+ng2,
                              nx3+1..nx3+ng3] do
        A[i,j,k] = if per then A[i, j, k-nx3] else 0.0;
    }
    A.updateFluff();
  }

  /* finite-volume Laplacian with the same metric factors as the
     conduction operator (valid in every geometry / grid law) */
  inline proc lapAt(const ref A: [SGD] real, i: int, j: int, k: int): real {
    var s = 0.0;
    if act1 {
      const gp = (A[i+1,j,k] - A[i,j,k])/(x1c(i+1) - x1c(i));
      const gm = (A[i,j,k] - A[i-1,j,k])/(x1c(i) - x1c(i-1));
      s += (fA1(i+1)*gp - fA1(i)*gm)*invV1(i);
    }
    if act2 {
      const m = if geom == Geom.polar || geom == Geom.spherical
                then x1c(i) else 1.0;
      const gp = (A[i,j+1,k] - A[i,j,k])/(m*(x2c(j+1) - x2c(j)));
      const gm = (A[i,j,k] - A[i,j-1,k])/(m*(x2c(j) - x2c(j-1)));
      s += (fA2(j+1)*gp - fA2(j)*gm)*invV2(j)*g2(i);
    }
    if act3 {
      const m = if geom == Geom.spherical then x1c(i)*sin(x2c(j)) else 1.0;
      const gp = (A[i,j,k+1] - A[i,j,k])/(m*(x3c(k+1) - x3c(k)));
      const gm = (A[i,j,k] - A[i,j,k-1])/(m*(x3c(k) - x3c(k-1)));
      s += (gp - gm)*invV3(k)*g3(i, j);
    }
    return s;
  }

  /* conjugate-gradient solve of lap(PHI) = b */
  proc solveGravity() {
    if sgFourPiG <= 0.0 then return;

    // right-hand side (mean removed when fully periodic)
    var rhobar = 0.0;
    if allPeriodic {
      const mtot = + reduce ([(i,j,k) in SGI] V[i,j,k](IRHO)*cellVol(i,j,k));
      const vtot = + reduce ([(i,j,k) in SGI] cellVol(i,j,k));
      rhobar = mtot/vtot;
    }
    forall (i, j, k) in SGI do
      BB[i,j,k] = sgFourPiG*(V[i,j,k](IRHO) - rhobar);

    fieldBC(PHI);
    forall (i, j, k) in SGI {
      RR[i,j,k] = BB[i,j,k] - lapAt(PHI, i, j, k);
      PP[i,j,k] = RR[i,j,k];
    }
    var rr = + reduce ([idx in SGI] RR[idx]*RR[idx]);
    const bb = max(+ reduce ([idx in SGI] BB[idx]*BB[idx]), 1.0e-300);

    var nIt = 0;
    while rr/bb > sgTol*sgTol && nIt < sgMaxIter {
      fieldBC(PP);
      forall (i, j, k) in SGI do
        AP[i,j,k] = lapAt(PP, i, j, k);
      const pap = + reduce ([idx in SGI] PP[idx]*AP[idx]);
      if abs(pap) < 1.0e-300 then break;
      const alpha = rr/pap;
      forall idx in SGI {
        PHI[idx] += alpha*PP[idx];
        RR[idx]  -= alpha*AP[idx];
      }
      const rrNew = + reduce ([idx in SGI] RR[idx]*RR[idx]);
      const beta = rrNew/rr;
      rr = rrNew;
      forall idx in SGI do
        PP[idx] = RR[idx] + beta*PP[idx];
      nIt += 1;
    }
    if nIt >= sgMaxIter then
      writeln("WARNING: self-gravity CG hit sgMaxIter (residual ",
              sqrt(rr/bb), ")");
    fieldBC(PHI);
  }

  /* acceleration g = -grad(Phi), central differences, physical lengths;
     components along the local orthonormal coordinate directions.
     The potential is passed in so the same code runs on the module's
     host array and on a GPU block's device copy. */
  inline proc sgAccel(const ref PH, i: int, j: int, k: int): 3*real {
    var g: 3*real;
    if act1 then
      g(0) = -(PH[i+1,j,k] - PH[i-1,j,k])/(x1c(i+1) - x1c(i-1));
    if act2 {
      const m = if geom == Geom.polar || geom == Geom.spherical
                then x1c(i) else 1.0;
      g(1) = -(PH[i,j+1,k] - PH[i,j-1,k])/(m*(x2c(j+1) - x2c(j-1)));
    }
    if act3 {
      const m = if geom == Geom.spherical then x1c(i)*sin(x2c(j)) else 1.0;
      g(2) = -(PH[i,j,k+1] - PH[i,j,k-1])/(m*(x3c(k+1) - x3c(k-1)));
    }
    return g;
  }
}
