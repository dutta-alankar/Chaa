# Boundary conditions

One condition per side: `--bcX1min`, `--bcX1max`, `--bcX2min`,
`--bcX2max`, `--bcX3min`, `--bcX3max`. Ghost cells are *owned* elements
of the padded distributed domain, so every boundary pass is an ordinary
distributed `forall` over a thin slab — there are no special cases at
partition edges.

## Available conditions

| value | what it does |
|---|---|
| `zero-gradient` | copies the nearest interior cell into the ghosts (free in/outflow). Default. `outflow` is accepted as a legacy alias. |
| `outflow-diode` | zero-gradient, but the ghost's normal velocity is clamped towards the exterior — material may leave, **nothing can flow back in**. Use on the open sides of wind tunnels and blast outflows (see `cloud-wind`). |
| `inflow-diode` | the mirror diode: zero-gradient with the normal velocity clamped towards the interior — **nothing can flow out**. |
| `periodic` | wraps the opposite side. |
| `reflect` | mirror image with the normal velocity flipped (solid wall / symmetry plane). |
| `axis` | `reflect` plus a flipped azimuthal velocity — the correct parity at coordinate axes (r=0 in cylindrical/polar/spherical, θ=0,π in spherical). Face areas vanish on the axis, so no flux leaks through by construction. |
| `inflow` | fixed ghost state from `--inRho --inVx1..3 --inPrs` (supersonic/forced inflow). |
| `shear-periodic` | shearing-box radial boundary (x1 sides, Cartesian, `--omegaRot>0`): the periodic image is sampled at the azimuthal offset ∓qΩL_x·t (linear interpolation) and the background velocity jump ±qΩL_x is added to v_y. |
| `userdef` | the problem's `userBC(side, t)` hook runs instead — arbitrary, time-dependent conditions ([write your own](../custom-problem.md)). |

Side ids for hooks: 0/1 = x1min/max, 2/3 = x2min/max, 4/5 = x3min/max.

## Order and corners

Sides are applied x1 → x2 → x3, each pass spanning the full ghost
extent of the other directions, which leaves edge/corner ghosts
consistent. After the primitive ghosts are set, ghost conservatives are
rebuilt (isothermal runs re-impose p = ρcs² first) and the halo caches
refresh — the only explicit communication in the code.

## Choosing between the outflow family

- `zero-gradient`: cheapest; can let ambient material drift in if the
  exterior pressure would push inward.
- `outflow-diode`: same cost, suppresses the spurious back-inflow —
  prefer it for blast waves leaving the box and wind-tunnel exits.
- `inflow` (fixed state): when the upstream state must stay pinned
  (wind tunnels, the post-shock feed of the double Mach reflection).
