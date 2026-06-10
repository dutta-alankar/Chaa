/* ProblemUtils.chpl — helpers shared by problem initialisers. */
module ProblemUtils {
  use Params, Grid;
  use Math;

  /* distance of a cell centre from the configured centre (cen1,cen2,cen3),
     measured in the physically mapped coordinates of the geometry */
  inline proc distFromCentre(i: int, j: int, k: int): real {
    const p = physPos(i, j, k);
    return sqrt((p(0)-cen1)**2 + (p(1)-cen2)**2 + (p(2)-cen3)**2);
  }
}
