# Tracer particles

Chaa carries massless Lagrangian tracer particles advected with the
gas. They need **no problem code**: enable them on any run with

```sh
--nParticles=256 [--partSeed=4321]
```

## What they do

- particles are scattered uniformly over the domain at start-up
  (deterministic for a given `partSeed`);
- each step they advance with **RK2 (midpoint)** in the velocity field,
  tri-linearly interpolated from cell centres;
- the interpolation works on **every grid law** — the closed-form
  coordinate maps are analytically inverted to locate a particle in
  index space (uniform, log and stretched directions alike);
- periodic (and shear-periodic) sides wrap positions; all other sides
  clamp at the domain edge;
- positions are written beside every field dump as

```
<outDir>/<problem>.particles.NNNN.txt     # columns: id  x  y  z
```

```python
import numpy as np
p = np.loadtxt("output/cloud.particles.0005.txt")
ids, x, y = p[:, 0], p[:, 1], p[:, 2]
```

## Accuracy

The `vortex-particles` CI case seeds 64 particles in the isentropic
vortex and integrates one full advection period: particles return to
their initial positions (periodic distance ≪ domain size; mean ≈ 0.06 L
at 64², dominated by interpolation error near the vortex core).

## Notes

- particle positions live in the *physically mapped* coordinates of the
  geometry (`physPos`), so they are directly comparable with the
  VTK/XDMF meshes;
- in multi-locale runs particles are advanced from locale 0 with remote
  velocity reads — fine at typical particle counts (≲10⁵); fully
  distributed particles are future work;
- velocity is frozen over the step (the field is interpolated at the
  step's final state), consistent with the RK2 accuracy of the gas.
