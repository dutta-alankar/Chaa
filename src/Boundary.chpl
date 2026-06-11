/* Boundary.chpl — ghost-cell boundary conditions.
 *
 * Ghost cells are owned elements of the padded distributed domain, so a
 * boundary pass is an ordinary (distributed) forall over a thin slab.
 * Sides are applied x1, then x2, then x3; each pass spans the *full*
 * extent of the other directions so edge/corner ghosts end up valid.
 *
 * side ids: 0:x1min 1:x1max 2:x2min 3:x2max 4:x3min 5:x3max
 */
module Boundary {
  use Params, Grid, State, Eos;
  use Problems;

  inline proc isInterior(idx: 3*int): bool {
    const (i, j, k) = idx;
    return i >= 1 && i <= nx1 && j >= 1 && j <= nx2 && k >= 1 && k <= nx3;
  }

  proc applyBCs(t: real) {
    if act1 { applySide(0, t); applySide(1, t); }
    if act2 { applySide(2, t); applySide(3, t); }
    if act3 { applySide(4, t); applySide(5, t); }

    forall (i, j, k) in DAll {
      if !isInterior((i, j, k)) {
        if eosCode == EOS_ISO then
          V[i,j,k](IPRS) = V[i,j,k](IRHO)*cs2At(i, j, k);
        U[i,j,k] = prim2cons(V[i,j,k]);
      }
    }
    syncHalos();
  }

  proc applySide(side: int, t: real) {
    const code = bcCode(side);
    if code == BC_USERDEF {
      problemUserBC(side, t);
      return;
    }

    const Dg = ghostSlab(side);
    // azimuthal velocity slot flipped by the "axis" condition
    const ivAxis = if geom == Geom.polar then IVX2 else IVX3;
    const nv = IVX1 + side/2;            // normal velocity slot

    forall (i, j, k) in Dg {
      select code {
        when BC_ZEROGRAD {
          V[i,j,k] = V[clampInt(side, i, j, k)];
        }
        when BC_OUT_DIODE {
          // zero-gradient, but never let material flow back in
          var w = V[clampInt(side, i, j, k)];
          if side % 2 == 0 then w(nv) = min(w(nv), 0.0);
                           else w(nv) = max(w(nv), 0.0);
          V[i,j,k] = w;
        }
        when BC_IN_DIODE {
          // zero-gradient, but never let material flow out
          var w = V[clampInt(side, i, j, k)];
          if side % 2 == 0 then w(nv) = max(w(nv), 0.0);
                           else w(nv) = min(w(nv), 0.0);
          V[i,j,k] = w;
        }
        when BC_PERIODIC {
          V[i,j,k] = V[periodicSrc(side, i, j, k)];
        }
        when BC_REFLECT, BC_AXIS {
          var w = V[mirrorSrc(side, i, j, k)];
          w(nv) = -w(nv);
          if code == BC_AXIS then w(ivAxis) = -w(ivAxis);
          V[i,j,k] = w;
        }
        when BC_INFLOW {
          V[i,j,k] = mkPrim(inRho, inVx1, inVx2, inVx3, inPrs);
        }
        otherwise { }
      }
    }
  }

  proc ghostSlab(side: int) {
    const r1 = 1-ng1..nx1+ng1,
          r2 = 1-ng2..nx2+ng2,
          r3 = 1-ng3..nx3+ng3;
    select side {
      when 0 do return DAll[1-ng1..0,        r2, r3];
      when 1 do return DAll[nx1+1..nx1+ng1,  r2, r3];
      when 2 do return DAll[r1, 1-ng2..0,        r3];
      when 3 do return DAll[r1, nx2+1..nx2+ng2,  r3];
      when 4 do return DAll[r1, r2, 1-ng3..0      ];
      otherwise do return DAll[r1, r2, nx3+1..nx3+ng3];
    }
  }

  inline proc clampInt(side: int, i: int, j: int, k: int): 3*int {
    select side {
      when 0 do return (1,   j, k);
      when 1 do return (nx1, j, k);
      when 2 do return (i, 1,   k);
      when 3 do return (i, nx2, k);
      when 4 do return (i, j, 1  );
      otherwise do return (i, j, nx3);
    }
  }

  inline proc periodicSrc(side: int, i: int, j: int, k: int): 3*int {
    select side {
      when 0 do return (i+nx1, j, k);
      when 1 do return (i-nx1, j, k);
      when 2 do return (i, j+nx2, k);
      when 3 do return (i, j-nx2, k);
      when 4 do return (i, j, k+nx3);
      otherwise do return (i, j, k-nx3);
    }
  }

  inline proc mirrorSrc(side: int, i: int, j: int, k: int): 3*int {
    select side {
      when 0 do return (1-i,       j, k);
      when 1 do return (2*nx1+1-i, j, k);
      when 2 do return (i, 1-j,       k);
      when 3 do return (i, 2*nx2+1-j, k);
      when 4 do return (i, j, 1-k      );
      otherwise do return (i, j, 2*nx3+1-k);
    }
  }
}
