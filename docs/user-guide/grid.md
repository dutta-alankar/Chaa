# The grid

Chaa runs on structured meshes whose coordinates follow a closed-form
**grid law** per direction — no coordinate arrays are stored, so every
metric query is a pure function of the index (communication-free on any
locale).

## Grid laws (`--gridX1`, `--gridX2`, `--gridX3`)

| law | faces | use when |
|---|---|---|
| `uniform` | \(x_i = x_{\min} + (i{-}1)\,L/n\) | default |
| `log` | \(x_i = x_{\min}\,(x_{\max}/x_{\min})^{(i-1)/n}\) | resolution ∝ 1/x — disks, winds, anything radial (needs \(x_{\min}>0\)) |
| `log-dec` | mirror of `log` | fine spacing at the *outer* edge |
| `stretch` | geometric progression, see below | concentrate resolution near one end with a controlled ratio |

## The stretched law in detail

`--gridX1=stretch --stretchX1=r --stretchUniX1=nu` builds a **uniform
anchor block plus a geometric progression**:

- a block of `nu` cells of uniform spacing \(h\) provides the starting
  value, and the remaining cells grow geometrically from it:
  \(h r, h r^2, h r^3, \dots\) (spacing is continuous across the
  junction);
- **r > 1** → the uniform block sits at the *beginning* (the fine end),
  cell sizes grow towards \(x_{\max}\);
- **r < 1** (with `nu > 0`) → the mirror image: cell sizes shrink from
  \(x_{\min}\) into a uniform block at the *end*;
- `nu = 0` (default) → a pure geometric progression across the whole
  direction (either orientation, depending on r).

\(h\) follows from the total length:
\(L = h\,[\,n_u + r(r^{n_s}-1)/(r-1)\,]\) with \(n_s = n - n_u\).

```sh
# 100 uniform cells near x=0, then 1%-per-cell growth to x=1:
./build/bin/chaa --problem=sod --nx1=400 --gridX1=stretch \
                 --stretchX1=1.01 --stretchUniX1=100
```

The `sod-stretch-anchor` CI case verifies the anchor block is uniform
to round-off, the ratio is honoured exactly, and the shock-tube
solution still matches the exact Riemann solution.

## What adapts automatically

All of these are evaluated with *local* spacings, so every grid law
works in every geometry with no further configuration:

- finite-volume areas/volumes and the well-balanced geometric sources,
- reconstruction stencils (PLM/LimO3/PPM/WENO-Z use the local cell
  width; formal high-order accuracy degrades gracefully on strongly
  stretched regions),
- the CFL condition (physical cell sizes, e.g. \(r\,\Delta\theta\)),
- viscous/conductive/diffusive fluxes,
- tracer-particle location (the closed-form laws are analytically
  invertible),
- output coordinates (txt/VTK/HDF5 carry the true cell faces/centres).

Validated hard: a 1D spherical Sedov blast on a `log` grid spanning a
**117× spacing range** lands on the analytic similarity radius to
0.28 % (`sedov-1d-log` in CI).

## Geometry

Coordinate meanings (PLUTO conventions) and the finite-volume metric
treatment are described in [Equations & geometry](equations.md):
`cartesian` (x,y,z), `cylindrical` (R,z; axisymmetric),
`polar` (R,φ,z), `spherical` (r,θ,φ) — any subset of directions active.
