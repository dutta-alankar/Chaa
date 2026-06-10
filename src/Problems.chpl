/* Problems.chpl — initial conditions, user-defined boundary conditions
 * and internal (immersed) boundaries for all bundled test problems.
 *
 *   sod            shock tube along x1 (any geometry: planar/cyl/sph)
 *   twoblast       Woodward & Colella interacting blast waves (1D)
 *   sedov          point explosion (energy deposit, any geometry/dim)
 *   blast          over-pressured circular/spherical region
 *   riemann2d      Lax & Liu four-quadrant 2D Riemann problem (config 3)
 *   dmr            double Mach reflection of a Mach 10 shock
 *   kh             Kelvin-Helmholtz instability
 *   rt             Rayleigh-Taylor instability (uniform gravity)
 *   vortex         isentropic vortex advection (smooth accuracy test)
 *   taylorCouette  viscous azimuthal flow between rotating cylinders
 *   cylinderFlow   viscous flow past an immersed cylinder
 */
module Problems {
  use Params, Grid, State, Eos;
  use Math;

  proc problemInit() {
    select problem {
      when "sod"           do initSod();
      when "twoblast"      do initTwoblast();
      when "sedov"         do initSedov();
      when "blast"         do initBlast();
      when "riemann2d"     do initRiemann2d();
      when "dmr"           do initDmr();
      when "kh"            do initKh();
      when "rt"            do initRt();
      when "vortex"        do initVortex();
      when "taylorCouette" do initTaylorCouette();
      when "cylinderFlow"  do initCylinderFlow();
      otherwise do halt("unknown problem: " + problem);
    }
  }

  inline proc distFromCentre(i: int, j: int, k: int): real {
    const p = physPos(i, j, k);
    return sqrt((p(0)-cen1)**2 + (p(1)-cen2)**2 + (p(2)-cen3)**2);
  }

  /* ------------------------------------------------------------------ */
  proc initSod() {
    forall (i, j, k) in DInt {
      if x1c(i) < sodX0 then
        V[i,j,k] = (sodRhoL, sodVxL, 0.0, 0.0, sodPrsL);
      else
        V[i,j,k] = (sodRhoR, sodVxR, 0.0, 0.0, sodPrsR);
    }
  }

  proc initTwoblast() {
    forall (i, j, k) in DInt {
      const x = x1c(i);
      var p = 0.01;
      if x < 0.1 then p = 1000.0;
      else if x > 0.9 then p = 100.0;
      V[i,j,k] = (1.0, 0.0, 0.0, 0.0, p);
    }
  }

  proc initSedov() {
    // measure the actual deposit volume so the analytic similarity
    // solution (total energy E0) applies in every geometry/dimension
    const vol = + reduce ([(i,j,k) in DInt]
                  (if distFromCentre(i,j,k) < sedovR0
                   then cellVol(i,j,k) else 0.0));
    if vol <= 0.0 then halt("sedov: no cells inside deposit radius");
    const pIn = (gam - 1.0)*sedovE0/vol;
    forall (i, j, k) in DInt {
      const p = if distFromCentre(i,j,k) < sedovR0 then pIn else sedovPrsAmb;
      V[i,j,k] = (sedovRhoAmb, 0.0, 0.0, 0.0, p);
    }
  }

  proc initBlast() {
    forall (i, j, k) in DInt {
      const inside = distFromCentre(i,j,k) < blastR0;
      V[i,j,k] = (if inside then blastRhoIn else blastRhoOut,
                  0.0, 0.0, 0.0,
                  if inside then blastPin else blastPout);
    }
  }

  /* Lax & Liu configuration 3, on [0,1]^2 split at (0.5, 0.5) */
  proc initRiemann2d() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      if x >= 0.5 && y >= 0.5 then
        V[i,j,k] = (1.5, 0.0, 0.0, 0.0, 1.5);
      else if x < 0.5 && y >= 0.5 then
        V[i,j,k] = (0.5323, 1.206, 0.0, 0.0, 0.3);
      else if x < 0.5 && y < 0.5 then
        V[i,j,k] = (0.138, 1.206, 1.206, 0.0, 0.029);
      else
        V[i,j,k] = (0.5323, 0.0, 1.206, 0.0, 0.3);
    }
  }

  /* ------------------------- double Mach reflection ----------------- */
  const dmrPost: StateVec = (8.0, 8.25*sin(pi/3.0), -8.25*cos(pi/3.0),
                             0.0, 116.5);
  const dmrPre:  StateVec = (1.4, 0.0, 0.0, 0.0, 1.0);

  proc initDmr() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const xs = 1.0/6.0 + y/tan(pi/3.0);
      V[i,j,k] = if x < xs then dmrPost else dmrPre;
    }
  }

  proc dmrBC(side: int, t: real) {
    if side == 2 {                       // bottom: post-shock for x<1/6,
      const Dg = DAll[1-ng1..nx1+ng1, 1-ng2..0, 1..1];
      forall (i, j, k) in Dg {           // reflecting wall beyond
        if x1c(i) < 1.0/6.0 {
          V[i,j,k] = dmrPost;
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
        V[i,j,k] = if x1c(i) < xs then dmrPost else dmrPre;
    }
  }

  /* ------------------------------------------------------------------ */
  proc initKh() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const inner = abs(y - 0.5) < 0.25;
      const rho = if inner then khRhoIn else khRhoOut;
      const vx  = if inner then khV0 else -khV0;
      // perturbation localised at the two shear layers
      const vy = khPert*sin(4.0*pi*x)
               * (exp(-(y-0.25)**2/0.005) + exp(-(y-0.75)**2/0.005));
      V[i,j,k] = (rho, vx, vy, 0.0, khPrs);
    }
  }

  proc initRt() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const rho = if y > 0.0 then rtRhoTop else rtRhoBot;
      const p = rtPrs0 + rho*grav2*y;    // piecewise hydrostatic
      // single-mode perturbation, zero at the walls (domain centred on 0)
      const Lx = x1max - x1min, Ly = x2max - x2min;
      const vy = rtPert*0.25*(1.0 + cos(2.0*pi*x/Lx))
                          *(1.0 + cos(2.0*pi*y/Ly));
      V[i,j,k] = (rho, 0.0, vy, 0.0, p);
    }
  }

  /* isentropic vortex advected diagonally; exact solution = initial
     condition translated by (t, t) */
  proc initVortex() {
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

  /* ------------------------- Taylor-Couette -------------------------- */
  proc initTaylorCouette() {
    const R1 = x1min, R2 = x1max;
    forall (i, j, k) in DInt {
      const R = x1c(i);
      const om = tcOmegaIn + (tcOmegaOut - tcOmegaIn)*(R - R1)/(R2 - R1);
      V[i,j,k] = (inRho, 0.0, 0.0, om*R, inPrs);
    }
  }

  /* no-slip rotating walls: v_phi fixed at the wall value, v_R
     reflected, density/pressure mirrored */
  proc tcBC(side: int, t: real) {
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

  /* ------------------------- flow past a cylinder -------------------- */
  proc initCylinderFlow() {
    forall (i, j, k) in DInt {
      const x = x1c(i), y = x2c(j);
      const inside = (x - cen1)**2 + (y - cen2)**2 < cylRad**2;
      solveMask[i,j,k] = !inside;
      V[i,j,k] = (inRho, if inside then 0.0 else inVx1, 0.0, 0.0, inPrs);
    }
    solveMask.updateFluff();
  }

  /* re-impose the solid state inside the cylinder after every stage */
  proc problemInternalBC(t: real) {
    if problem == "cylinderFlow" {
      forall idx in DInt {
        if !solveMask[idx] {
          V[idx] = (inRho, 0.0, 0.0, 0.0, inPrs);
          U[idx] = prim2cons(V[idx]);
        }
      }
    }
  }

  proc problemUserBC(side: int, t: real) {
    select problem {
      when "dmr"           do dmrBC(side, t);
      when "taylorCouette" do tcBC(side, t);
      otherwise do halt("problem " + problem +
                        " does not define user boundary conditions");
    }
  }
}
