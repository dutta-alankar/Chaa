/* Grid.chpl — structured mesh and curvilinear geometry factors.
 *
 * Each coordinate direction follows one of four closed-form grid laws
 * (Idefix-style):
 *   uniform   constant spacing
 *   log       x_f(i) = xmin (xmax/xmin)^((i-1)/n)   — spacing grows ~x
 *   log-dec   mirrored log                          — spacing shrinks ~x
 *   stretch   spacings in geometric progression with ratio stretchX*
 *             (>1 grows, <1 shrinks)
 * Coordinates and all metric factors remain closed-form functions of
 * the index — no coordinate arrays are stored, so every geometric query
 * is communication-free on any locale (performance portable).
 *
 * Coordinate meaning per geometry (PLUTO conventions):
 *   cartesian   : x1=x,  x2=y,     x3=z
 *   cylindrical : x1=R,  x2=z,     x3=phi (axisymmetric; nx3 must be 1,
 *                                          v3 = v_phi evolves passively)
 *   polar       : x1=R,  x2=phi,   x3=z
 *   spherical   : x1=r,  x2=theta, x3=phi
 *
 * The finite-volume divergence in direction d is
 *     (A_d F)_right - (A_d F)_left) * invV_d   [ * g2(i) or g3(i,j) ]
 * and the geometric source terms below are evaluated with the *same*
 * centroid factors so that a uniform pressure field is balanced to
 * machine precision in every geometry.
 */
module Grid {
  use Params;
  use Math;
  import CompileParams.NG;   // ghost layers (compile_params.chpl)

  const ng1 = if act1 then NG else 0,
        ng2 = if act2 then NG else 0,
        ng3 = if act3 then NG else 0;

  // mean spacings (exact for uniform; used only as scales elsewhere)
  const dx1 = (x1max - x1min)/nx1,
        dx2 = (x2max - x2min)/nx2,
        dx3 = (x3max - x3min)/nx3;

  /* stretched law, forward orientation: a uniform block of nu cells of
     spacing h at the start, then spacings in geometric progression
     h*r, h*r^2, ... (continuous across the junction).  nu = 0 gives a
     pure geometric progression starting from h.  h follows from the
     total length. */
  inline proc stretchH(L: real, n: int, nu: int, r: real): real {
    const ns = n - nu;
    if nu == 0 then
      return L*(r - 1.0)/(r**(n:real) - 1.0);
    return L/(nu + r*(r**(ns:real) - 1.0)/(r - 1.0));
  }

  inline proc stretchFwd(xmin: real, L: real, n: int, nu: int,
                         r: real, i: int): real {
    const h = stretchH(L, n, nu, r);
    if nu == 0 then
      return xmin + h*(r**(i - 1.0) - 1.0)/(r - 1.0);
    if i <= nu + 1 then return xmin + (i - 1.0)*h;
    return xmin + nu*h + h*r*(r**(i - 1.0 - nu) - 1.0)/(r - 1.0);
  }

  inline proc stretchFwdIdx(xmin: real, L: real, n: int, nu: int,
                            r: real, x: real): real {
    const h = stretchH(L, n, nu, r);
    if nu == 0 then
      return 1.0 + log(1.0 + (x - xmin)*(r - 1.0)/h)/log(r);
    if x <= xmin + nu*h then return 1.0 + (x - xmin)/h;
    return nu + 1.0
           + log(1.0 + (x - xmin - nu*h)*(r - 1.0)/(h*r))/log(r);
  }

  /* generic face coordinate under a grid law; face i is the *left*
     face of cell i, i in 1..n+1 (ghost faces extrapolate the law).
     stretch: ratio r > 1 grows away from a uniform block of nu cells
     at the *start*; r < 1 with nu > 0 is the mirror image — spacing
     shrinks into a uniform block of nu cells at the *end*. */
  inline proc lawFace(code: int, xmin: real, xmax: real, n: int,
                      r: real, nu: int, i: int): real {
    select code {
      when GRID_UNIFORM do return xmin + (i-1.0)*(xmax - xmin)/n;
      when GRID_LOG     do return xmin*(xmax/xmin)**((i-1.0)/n);
      when GRID_LOGDEC  do
        return xmin + xmax - xmin*(xmax/xmin)**((n+1.0-i)/n);
      when GRID_STRETCH {
        if abs(r - 1.0) < 1.0e-14 then
          return xmin + (i-1.0)*(xmax - xmin)/n;
        if r > 1.0 || nu == 0 then
          return stretchFwd(xmin, xmax - xmin, n, nu, r, i);
        // r < 1 with a uniform block: mirror of the 1/r forward law
        return xmin + xmax
               - stretchFwd(xmin, xmax - xmin, n, nu, 1.0/r, n + 2 - i);
      }
      otherwise do return xmin + (i-1.0)*(xmax - xmin)/n;
    }
    return 0.0;
  }

  /* inverse of the face law: fractional face index of coordinate x
     (used by the tracer particles to locate themselves) */
  inline proc lawIndex(code: int, xmin: real, xmax: real, n: int,
                       r: real, nu: int, x: real): real {
    select code {
      when GRID_UNIFORM do return 1.0 + (x - xmin)*n/(xmax - xmin);
      when GRID_LOG     do
        return 1.0 + n*log(x/xmin)/log(xmax/xmin);
      when GRID_LOGDEC  do
        return (n + 1.0) - n*log((xmin + xmax - x)/xmin)/log(xmax/xmin);
      when GRID_STRETCH {
        if abs(r - 1.0) < 1.0e-14 then
          return 1.0 + (x - xmin)*n/(xmax - xmin);
        if r > 1.0 || nu == 0 then
          return stretchFwdIdx(xmin, xmax - xmin, n, nu, r, x);
        return (n + 2.0)
               - stretchFwdIdx(xmin, xmax - xmin, n, nu, 1.0/r,
                               xmin + xmax - x);
      }
      otherwise do return 1.0 + (x - xmin)*n/(xmax - xmin);
    }
    return 1.0;
  }

  inline proc x1f(i: int): real do
    return lawFace(gridCode(0), x1min, x1max, nx1, stretchX1,
                   stretchUniX1, i);
  inline proc x2f(j: int): real do
    return lawFace(gridCode(1), x2min, x2max, nx2, stretchX2,
                   stretchUniX2, j);
  inline proc x3f(k: int): real do
    return lawFace(gridCode(2), x3min, x3max, nx3, stretchX3,
                   stretchUniX3, k);

  inline proc x1c(i: int): real do return 0.5*(x1f(i) + x1f(i+1));
  inline proc x2c(j: int): real do return 0.5*(x2f(j) + x2f(j+1));
  inline proc x3c(k: int): real do return 0.5*(x3f(k) + x3f(k+1));

  /* local cell widths */
  inline proc dx1At(i: int): real do return x1f(i+1) - x1f(i);
  inline proc dx2At(j: int): real do return x2f(j+1) - x2f(j);
  inline proc dx3At(k: int): real do return x3f(k+1) - x3f(k);

  /* distance between adjacent cell centres across face q along dir */
  inline proc dcAt(dir: int, q: int): real {
    if dir == 0 then return x1c(q) - x1c(q-1);
    if dir == 1 then return x2c(q) - x2c(q-1);
    return x3c(q) - x3c(q-1);
  }

  /* ---- direction 1 (x | R | r) ---- */
  inline proc fA1(i: int): real {
    select geom {
      when Geom.cartesian                  do return 1.0;
      when Geom.cylindrical, Geom.polar    do return abs(x1f(i));
      when Geom.spherical                  do return x1f(i)**2;
    }
    return 1.0;
  }

  inline proc invV1(i: int): real {
    select geom {
      when Geom.cartesian do return 1.0/dx1At(i);
      when Geom.cylindrical, Geom.polar {
        const rm = x1f(i), rp = x1f(i+1);
        return 2.0/abs(rp**2 - rm**2);
      }
      when Geom.spherical {
        const rm = x1f(i), rp = x1f(i+1);
        return 3.0/(rp**3 - rm**3);
      }
    }
    return 1.0/dx1At(i);
  }

  /* radius used in geometric source terms; chosen so that uniform
     pressure exactly balances the area-weighted flux divergence */
  inline proc rGeo(i: int): real {
    select geom {
      when Geom.cylindrical, Geom.polar do return x1c(i);
      when Geom.spherical {
        const rm = x1f(i), rp = x1f(i+1);
        return (2.0/3.0)*(rp**3 - rm**3)/(rp**2 - rm**2);
      }
      otherwise do return 1.0;
    }
    return 1.0;
  }

  /* ---- direction 2 (y | z | phi | theta) ---- */
  inline proc fA2(j: int): real {
    if geom == Geom.spherical then return sin(x2f(j));
    return 1.0;
  }

  inline proc invV2(j: int): real {
    if geom == Geom.spherical then
      return 1.0/(cos(x2f(j)) - cos(x2f(j+1)));
    return 1.0/dx2At(j);
  }

  /* extra 1/r factor multiplying the direction-2 divergence */
  inline proc g2(i: int): real {
    select geom {
      when Geom.polar     do return 1.0/x1c(i);
      when Geom.spherical do return 1.0/rGeo(i);
      otherwise           do return 1.0;
    }
    return 1.0;
  }

  /* well-balanced cot(theta) at the cell centroid (spherical) */
  inline proc cotGeo(j: int): real {
    const tm = x2f(j), tp = x2f(j+1);
    return (sin(tp) - sin(tm))/(cos(tm) - cos(tp));
  }

  /* mean sin(theta) over the cell (spherical) */
  inline proc sinGeo(j: int): real {
    if geom == Geom.spherical then
      return (cos(x2f(j)) - cos(x2f(j+1)))/dx2At(j);
    return 1.0;
  }

  /* ---- direction 3 (z | phi) ---- */
  inline proc invV3(k: int): real do return 1.0/dx3At(k);

  /* extra metric factor multiplying the direction-3 divergence */
  inline proc g3(i: int, j: int): real {
    select geom {
      when Geom.spherical do return 1.0/(rGeo(i)*sinGeo(j));
      otherwise           do return 1.0;   // polar x3=z, cartesian x3=z
    }
    return 1.0;
  }

  /* ---- physical (linear) cell sizes, for the CFL condition ---- */
  inline proc dl1(i: int): real do return dx1At(i);

  inline proc dl2(i: int, j: int): real {
    select geom {
      when Geom.polar     do return x1c(i)*dx2At(j);
      when Geom.spherical do return x1c(i)*dx2At(j);
      otherwise           do return dx2At(j);
    }
    return dx2At(j);
  }

  inline proc dl3(i: int, j: int, k: int): real {
    select geom {
      when Geom.spherical do return x1c(i)*sinGeo(j)*dx3At(k);
      otherwise           do return dx3At(k);   // polar: x3 = z
    }
    return dx3At(k);
  }

  /* ---- true cell volume; inactive angular dimensions contribute their
     full measure (2 for cos(theta), 2*pi for phi) so that e.g. a 1D
     spherical Sedov deposit normalises against the real 4*pi volume ---- */
  inline proc cellVol(i: int, j: int, k: int): real {
    var f1, f2, f3: real;
    select geom {
      when Geom.cartesian {
        f1 = dx1At(i);
        f2 = if act2 then dx2At(j) else 1.0;
        f3 = if act3 then dx3At(k) else 1.0;
      }
      when Geom.cylindrical {            // x1=R, x2=z, x3=phi(axisym)
        f1 = 0.5*abs(x1f(i+1)**2 - x1f(i)**2);
        f2 = if act2 then dx2At(j) else 1.0;
        f3 = if act3 then dx3At(k) else 2.0*pi;
      }
      when Geom.polar {                  // x1=R, x2=phi, x3=z
        f1 = 0.5*abs(x1f(i+1)**2 - x1f(i)**2);
        f2 = if act2 then dx2At(j) else 2.0*pi;
        f3 = if act3 then dx3At(k) else 1.0;
      }
      when Geom.spherical {              // x1=r, x2=theta, x3=phi
        f1 = (x1f(i+1)**3 - x1f(i)**3)/3.0;
        f2 = if act2 then (cos(x2f(j)) - cos(x2f(j+1))) else 2.0;
        f3 = if act3 then dx3At(k) else 2.0*pi;
      }
    }
    return f1*f2*f3;
  }

  /* ---- physical position of a cell centre, mapped to Cartesian-like
     coordinates (used for initial conditions / distances) ---- */
  inline proc physPos(i: int, j: int, k: int): 3*real {
    select geom {
      when Geom.cartesian {
        return (x1c(i), if act2 then x2c(j) else 0.0,
                        if act3 then x3c(k) else 0.0);
      }
      when Geom.cylindrical {            // meridional (R, z) plane
        return (x1c(i), if act2 then x2c(j) else 0.0, 0.0);
      }
      when Geom.polar {
        const R = x1c(i);
        const ph = if act2 then x2c(j) else 0.0;
        return (R*cos(ph), R*sin(ph), if act3 then x3c(k) else 0.0);
      }
      when Geom.spherical {
        const r  = x1c(i);
        const th = if act2 then x2c(j) else 0.5*pi;
        const ph = if act3 then x3c(k) else 0.0;
        return (r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th));
      }
    }
    return (x1c(i), 0.0, 0.0);
  }

  /* cylindrical radius of a cell centre (distance from the rotation
     axis), used by the locally-isothermal sound-speed profile */
  inline proc cylRadiusAt(i: int, j: int, k: int): real {
    select geom {
      when Geom.cylindrical, Geom.polar do return abs(x1c(i));
      when Geom.spherical do
        return x1c(i)*(if act2 then sin(x2c(j)) else 1.0);
      otherwise do return 1.0;
    }
    return 1.0;
  }

  /* isothermal sound speed squared: cs^2 = csIso^2 * R_cyl^(2*csSlope)
     (csSlope = 0 gives a globally constant sound speed; csSlope = -0.5
     with central gravity GM=1 gives cs = csIso * v_K, i.e. a constant
     disk aspect ratio) */
  inline proc cs2At(i: int, j: int, k: int): real {
    if csSlope == 0.0 then return csIso*csIso;
    return csIso*csIso * cylRadiusAt(i, j, k)**(2.0*csSlope);
  }

  /* physical node (cell-corner) position mapped to output x/y/z space;
     used by the VTK / XDMF writers */
  proc nodePos(i: int, j: int, k: int): 3*real {
    select geom {
      when Geom.cartesian {
        return (x1f(i), if act2 then x2f(j) else 0.0,
                        if act3 then x3f(k) else 0.0);
      }
      when Geom.cylindrical {                 // meridional (R, z) plane
        return (x1f(i), if act2 then x2f(j) else 0.0, 0.0);
      }
      when Geom.polar {
        const R = x1f(i);
        const ph = if act2 then x2f(j) else 0.0;
        return (R*cos(ph), R*sin(ph), if act3 then x3f(k) else 0.0);
      }
      when Geom.spherical {
        const r = x1f(i);
        const th = if act2 then x2f(j) else 0.5*pi;
        if !act3 then                          // meridional plane
          return (r*sin(th), r*cos(th), 0.0);
        const ph = x3f(k);
        return (r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th));
      }
    }
    return (x1f(i), 0.0, 0.0);
  }
}
