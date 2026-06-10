/* Grid.chpl — uniform structured mesh and curvilinear geometry factors.
 *
 * The grid is uniform in each *coordinate* (x, R or r, z, theta, phi), so
 * coordinates and all metric factors are closed-form functions of the
 * index.  No coordinate arrays are stored, which makes every geometric
 * query communication-free on any locale (performance portable).
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

  param NG = 2;          // ghost layers (2 needed for linear reconstruction)

  const ng1 = if act1 then NG else 0,
        ng2 = if act2 then NG else 0,
        ng3 = if act3 then NG else 0;

  const dx1 = (x1max - x1min)/nx1,
        dx2 = (x2max - x2min)/nx2,
        dx3 = (x3max - x3min)/nx3;

  /* face and centre coordinates: face i is the *left* face of cell i */
  inline proc x1f(i: int): real do return x1min + (i-1)*dx1;
  inline proc x2f(j: int): real do return x2min + (j-1)*dx2;
  inline proc x3f(k: int): real do return x3min + (k-1)*dx3;
  inline proc x1c(i: int): real do return x1min + (i-0.5)*dx1;
  inline proc x2c(j: int): real do return x2min + (j-0.5)*dx2;
  inline proc x3c(k: int): real do return x3min + (k-0.5)*dx3;

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
      when Geom.cartesian do return 1.0/dx1;
      when Geom.cylindrical, Geom.polar {
        const rm = x1f(i), rp = x1f(i+1);
        return 2.0/abs(rp**2 - rm**2);
      }
      when Geom.spherical {
        const rm = x1f(i), rp = x1f(i+1);
        return 3.0/(rp**3 - rm**3);
      }
    }
    return 1.0/dx1;
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
    return 1.0/dx2;
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
      return (cos(x2f(j)) - cos(x2f(j+1)))/dx2;
    return 1.0;
  }

  /* ---- direction 3 (z | phi) ---- */
  inline proc invV3(): real do return 1.0/dx3;

  /* extra metric factor multiplying the direction-3 divergence */
  inline proc g3(i: int, j: int): real {
    select geom {
      when Geom.spherical do return 1.0/(rGeo(i)*sinGeo(j));
      otherwise           do return 1.0;   // polar x3=z, cartesian x3=z
    }
    return 1.0;
  }

  /* ---- physical (linear) cell sizes, for the CFL condition ---- */
  inline proc dl1(): real do return dx1;

  inline proc dl2(i: int): real {
    select geom {
      when Geom.polar     do return x1c(i)*dx2;
      when Geom.spherical do return x1c(i)*dx2;
      otherwise           do return dx2;
    }
    return dx2;
  }

  inline proc dl3(i: int, j: int): real {
    select geom {
      when Geom.spherical do return x1c(i)*sinGeo(j)*dx3;
      otherwise           do return dx3;   // polar: x3 = z
    }
    return dx3;
  }

  /* ---- true cell volume; inactive angular dimensions contribute their
     full measure (2 for cos(theta), 2*pi for phi) so that e.g. a 1D
     spherical Sedov deposit normalises against the real 4*pi volume ---- */
  inline proc cellVol(i: int, j: int, k: int): real {
    var f1, f2, f3: real;
    select geom {
      when Geom.cartesian {
        f1 = dx1;
        f2 = if act2 then dx2 else 1.0;
        f3 = if act3 then dx3 else 1.0;
      }
      when Geom.cylindrical {            // x1=R, x2=z, x3=phi(axisym)
        f1 = 0.5*abs(x1f(i+1)**2 - x1f(i)**2);
        f2 = if act2 then dx2 else 1.0;
        f3 = if act3 then dx3 else 2.0*pi;
      }
      when Geom.polar {                  // x1=R, x2=phi, x3=z
        f1 = 0.5*abs(x1f(i+1)**2 - x1f(i)**2);
        f2 = if act2 then dx2 else 2.0*pi;
        f3 = if act3 then dx3 else 1.0;
      }
      when Geom.spherical {              // x1=r, x2=theta, x3=phi
        f1 = (x1f(i+1)**3 - x1f(i)**3)/3.0;
        f2 = if act2 then (cos(x2f(j)) - cos(x2f(j+1))) else 2.0;
        f3 = if act3 then dx3 else 2.0*pi;
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
