# Tracer particles

Chaa carries massless Lagrangian tracer particles advected with the
gas. They need **no problem code**: enable them on any run with

```sh
--nParticles=256 [--partSeed=4321]
```

## What they do

- particles start from the positions given by the problem's
  `problemParticleInit` hook — the **default is a uniform random
  scatter over the whole domain** (deterministic for a given
  `partSeed`);
- each step they advance with **RK2 (midpoint)** in the velocity field,
  tri-linearly interpolated from cell centres;
- the interpolation works on **every grid law** — the closed-form
  coordinate maps are analytically inverted to locate a particle in
  index space (uniform, log and stretched directions alike);
- periodic and shear-periodic sides wrap positions (shear-periodic x1
  crossings get the azimuthal image offset of the shearing box); all
  other sides clamp at the domain edge;
- positions are written beside every field dump as

```
<outDir>/<problem>.particles.NNNN.txt     # columns: id  x  y  z
```

```python
import numpy as np
p = np.loadtxt("output/cloud.particles.0005.txt")
ids, x, y = p[:, 0], p[:, 1], p[:, 2]
```

or use the bundled plotting tools, which overplot particles on 2D field
maps automatically (see [Plotting & analysis](plotting.md)):

```sh
python tools/plot_fields.py output/
```

## Custom initial positions

Any problem can seed the particles itself by registering a hook in
`Problems.problemParticleInit` — fill the position array and return
`true` (return `false` to keep the default random scatter):

```chapel
// in your problem module
proc particleInit(ref pos: [?D] 3*real): bool {
  for p in D {
    const th = 2.0*pi*p/D.size;
    pos[p] = (cen1 + 1.5*cos(th), cen2 + 1.5*sin(th), 0.0);
  }
  return true;
}
```

```chapel
// in src/Problems.chpl
proc problemParticleInit(ref pos: [] 3*real): bool {
  select problem {
    when "myproblem" do return MyProblem.particleInit(pos);
    otherwise do return false;
  }
}
```

The bundled example is the isentropic vortex: with `--partRingR=1.5`
its hook seeds the particles on a circle around the vortex centre
(instead of the random scatter), and the `vortex-particles-ring` CI
case checks that the ring stays a ring and rotates at the analytic
rate. Positions handed back by the hook are wrapped/clamped into the
domain before use.

## Distributed implementation

Particles are **fully distributed** (owner-computes, the same strategy
as AthenaPK's swarms):

- every locale keeps the particles that live inside its block of the
  distributed grid in a local bag;
- the RK2 update reads only locale-local field data — interpolation
  stencils that poke across a block edge are served by the same
  StencilDist fluff caches as the hydro stencils;
- after each step, particles that crossed a block boundary are handed
  to the receiving locale (a per-locale outbox, exchanged once per
  step);
- output gathers positions by particle id, so dump files are identical
  regardless of the locale count.

A 4-locale run reproduces single-locale trajectories to 10⁻¹² (the
run-to-run reduction-order difference); the bookkeeping is checked at
every dump — a lost particle halts the run rather than silently
corrupting statistics.

## Accuracy

Two CI cases validate the particles quantitatively in the isentropic
vortex (period t = 10 at 64²):

- `vortex-particles` — 64 randomly scattered tracers return to their
  initial positions after one advection period (mean periodic distance
  ≈ 0.06 of the domain, dominated by interpolation error near the
  vortex core);
- `vortex-particles-ring` — 64 hook-seeded tracers on a ring of radius
  1.5 around the vortex: the ring radius is preserved to 0.4 % and the
  measured rotation (4.19 rad) matches the analytic angular rate
  ω(R) = β/2π·e^{(1−R²)/2} (4.26 rad) to ~2 %.

## Notes

- particle positions live in the *physically mapped* coordinates of the
  geometry (`physPos`), so they are directly comparable with the
  VTK/XDMF meshes;
- velocity is frozen over the step (the field is interpolated at the
  step's final state), consistent with the RK2 accuracy of the gas.
