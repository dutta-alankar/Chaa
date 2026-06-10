/* State.chpl — distributed field arrays.
 *
 * All fields live on a single block-distributed domain with stencil
 * fluff (halo) regions.  Chapel's StencilDist keeps cached copies of
 * neighbouring elements on each locale; `updateFluff()` refreshes the
 * caches and is the *only* explicit communication in the whole code.
 * Physical ghost cells are ordinary owned elements of the padded domain.
 */
module State {
  use Params, Grid;
  use StencilDist;

  const fullSpace = {1-ng1..nx1+ng1, 1-ng2..nx2+ng2, 1-ng3..nx3+ng3};

  const DAll = fullSpace dmapped new stencilDist(fullSpace,
                                                 fluff=(ng1, ng2, ng3));

  /* interior (physical) cells */
  const DInt = DAll[1..nx1, 1..nx2, 1..nx3];

  var V:    [DAll] StateVec;   // primitive  (rho, v1, v2, v3, p)
  var U:    [DAll] StateVec;   // conservative (rho, m1, m2, m3, E)
  var U0:   [DAll] StateVec;   // RK stage buffer
  var RHS:  [DAll] StateVec;   // accumulated right-hand side
  var FLX:  [DAll] StateVec;   // per-sweep face fluxes (face = left face)

  /* false marks cells inside an immersed solid (e.g. the cylinder) */
  var solveMask: [DAll] bool = true;

  inline proc syncHalos() {
    V.updateFluff();
    U.updateFluff();
  }
}
