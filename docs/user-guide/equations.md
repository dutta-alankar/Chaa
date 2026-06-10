# Equations & geometry

## Equations

Chaa solves the compressible Euler equations in conservative form,

$$
\partial_t \rho + \nabla\!\cdot(\rho\mathbf{v}) = 0,
$$

$$
\partial_t(\rho\mathbf{v}) + \nabla\!\cdot(\rho\mathbf{v}\mathbf{v} + p\,\mathbb{I})
  = \rho\,\mathbf{g} + \nabla\!\cdot\boldsymbol{\tau},
$$

$$
\partial_t E + \nabla\!\cdot\big[(E+p)\,\mathbf{v}\big]
  = \rho\,\mathbf{v}\!\cdot\!\mathbf{g}
  + \nabla\!\cdot(\boldsymbol{\tau}\cdot\mathbf{v})
  + \nabla\!\cdot(\kappa\nabla T),
$$

with \(E = p/(\gamma-1) + \tfrac12\rho v^2\) and \(T = p/\rho\).
Optional right-hand sides: constant gravity, central point-mass gravity
\(\mathbf g = -GM\,\hat r/(r^2+\epsilon^2)^{3/2}\), explicit viscosity
\(\boldsymbol\tau = \mu(\nabla\mathbf v + \nabla\mathbf v^{T}
- \tfrac23\nabla\!\cdot\!\mathbf v\,\mathbb I)\) and explicit thermal
conduction.

With `--eos=isothermal` the energy equation is dropped and
\(p = \rho\,c_s^2\) is enforced after every stage, with optionally a
locally isothermal profile \(c_s = c_{s,0}\,R_{\rm cyl}^{s}\)
(`csIso`, `csSlope`) — constant disk aspect ratio for \(s=-1/2\).

## Coordinates

Chaa follows PLUTO's axis conventions:

| `geometry` | x1 | x2 | x3 | notes |
|---|---|---|---|---|
| `cartesian` | x | y | z | |
| `cylindrical` | R | z | φ | axisymmetric: `nx3=1`, v₃=v_φ evolves passively |
| `polar` | R | φ | z | |
| `spherical` | r | θ | φ | θ ∈ [0, π] |

Any subset of dimensions can be active: 1D (x or r), 2D (x,y / R,z /
R,φ / r,θ) and 3D (x,y,z / R,φ,z / r,θ,φ).

## Finite-volume discretisation

The divergence in direction \(d\) is discretised with face areas
\(A_d\) and cell volumes per direction, e.g. radially in spherical
coordinates

$$
(\nabla\cdot F)_1 \approx
\frac{A_{i+1/2}F_{i+1/2} - A_{i-1/2}F_{i-1/2}}{\Delta V_i},
\qquad A = r^2,\;
\Delta V = \tfrac13(r_{i+1/2}^3 - r_{i-1/2}^3).
$$

### Well-balanced geometric sources

Curvilinear momentum equations pick up geometric source terms (e.g.
\((2p + \rho(v_\theta^2{+}v_\phi^2))/r\) radially in spherical
coordinates). Chaa evaluates them with **the same centroid factors as
the flux divergence** — the centroid radius
\(\tilde r = \tfrac23(r_+^3-r_-^3)/(r_+^2-r_-^2)\) and the centroid
cotangent \(\widetilde{\cot\theta} =
(\sin\theta_+-\sin\theta_-)/(\cos\theta_--\cos\theta_+)\) — so that a
uniform-pressure fluid is balanced to machine precision. In practice: a
2D (R,z) Sedov blast stays spherical to 0.01 %, and a 3D (r,θ,φ) blast
has a shock radius independent of angle to round-off.

The compensating \(p\,\cot\theta/r\) source belongs to the θ flux
operator and is therefore dropped automatically when the θ dimension is
inactive (equatorial symmetry).

## Boundary conditions

Ghost cells are real (owned) elements of the padded distributed domain;
each side applies one of `periodic`, `outflow`, `reflect`, `axis`,
`inflow` or a problem-defined `userdef` hook. At coordinate axes
(`axis`): the face area vanishes there, so no flux leaks through r = 0
or θ = 0, π by construction.

## Viscosity and conduction

- **Cartesian:** the full stress tensor at faces, with transverse
  derivatives from the halo stencil; viscous heating \(\tau\cdot v\)
  included in the energy flux.
- **Cylindrical/polar:** the dominant \(\tau_{R\phi} =
  \mu\,R\,\partial_R(v_\phi/R)\) term, which drives Taylor–Couette
  flow; validated to 0.35 % against the analytic Couette profile.
- **Thermal conduction:** \(-\kappa\,\partial T/\partial \ell\) along
  every active direction; the decaying entropy-mode test reproduces the
  analytic rate \(\Gamma = \kappa\,\tfrac{\gamma-1}{\gamma}k^2\) to ~3 %.

Both impose the explicit diffusion time-step limit
\(\Delta t \le \tfrac12\,\mathrm{cflVisc}\,\Delta\ell^2/(\chi\,n_{\rm dim})\).
