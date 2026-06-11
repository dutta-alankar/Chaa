# Cross-code validation

Chaa is validated not only against exact/similarity solutions but
**against two production codes run locally at matched configurations**:
[Idefix](https://github.com/idefix-code/idefix) (Kokkos, serial build)
and [AthenaPK](https://github.com/parthenon-hpc-lab/athenapk)
(Parthenon, serial build). Small reference profiles frozen from those
runs live in `tests/reference/` and are re-checked by the `ref-*` CI
cases on every push — CI never needs to build the external codes.

## Method

For each overlapping test the *same* grid, domain, γ/cs, CFL-equivalent
time step, reconstruction order and final time were used; remaining
differences (Riemann solver flavour, limiter details) are
scheme-level, so agreement is expected at or below each code's distance
to the exact solution.

## Results

| test | matched configuration | metric | result |
|---|---|---|---|
| Sod | 500 cells, plm/rk2, cfl 0.8 (Idefix: roe; Chaa: hllc) | L1(ρ) | **2.7×10⁻⁴** (vs 1.6×10⁻³ to the exact solution) |
| Sod, isothermal | 500 cells, cs=1 | L1(ρ) | **3.8×10⁻⁴** |
| Mach reflection (DMR) | 480×120, t=0.2 | rel. L1(ρ), full field | **0.31 %**; ρmax 22.00 vs 22.01 |
| Sedov 3D | 64³, γ=5/3, hll, periodic | rel. L1 radial ρ profile | **3.3×10⁻⁴**; peaks equal to 3 decimals |
| KHI (isothermal) | 256×64, cs=10, t=1 | ⟨v_y²⟩ | within **13 %** (nonlinear stage of an instability) |
| thermal diffusion | 500 cells, κ=0.1 | δT decay over t=0.2 | 0.755 vs 0.729 (**3 %**) |
| Sod (AthenaPK) | 256 cells, hlle/plm/vl2, t=0.4 | L1(ρ) | **5.6×10⁻⁴** |
| linear wave (AthenaPK) | 128 cells, 5 periods, vl2 | L1(δρ)/amp | 1.1×10⁻² vs 5.3×10⁻³ (same order; limiter detail) |

Two instructive non-discrepancies found during the campaign:

- **Idefix's thermalDiffusion test zeroes the velocity field every
  step** (an internal boundary labelled "cancel any motion") to isolate
  pure conduction. Chaa's run keeps the physical acoustic response, so
  the *density* fields differ while the conduction rates agree — the
  velocity Chaa develops matches the analytic acoustic estimate.
- **AthenaPK's rk1 linear-wave configuration is linearly unstable**
  (forward Euler + PLM); both codes blow up on it, and their
  instability amplitudes after 5 periods agree to 2 % — the codes even
  fail identically. With `vl2` both are stable and dissipate at the
  same order.

**No Chaa code changes were required**: every apparent mismatch traced
to configuration or diagnostic conventions, not solver differences.
One convention worth knowing: Idefix/AthenaPK define the CFL number
against the *sum* of directional signal speeds, Chaa against the
per-direction minimum — a 3D Idefix `CFL 0.9` corresponds to roughly
`--cfl=0.3` in Chaa.

## Reproducing the campaign

```sh
# Idefix (serial):
git clone --recurse-submodules https://github.com/idefix-code/idefix
cd idefix/test/HD/sod && IDEFIX_DIR=$PWD/../../.. cmake $IDEFIX_DIR && make -j
./idefix -i idefix.ini

# AthenaPK (serial):
git clone --recurse-submodules https://github.com/parthenon-hpc-lab/athenapk
cmake -Bbuild -DPARTHENON_DISABLE_MPI=ON -DKokkos_ENABLE_SERIAL=ON \
      -DPARTHENON_ENABLE_PYTHON_MODULE_CHECK=OFF && cmake --build build -j
build/bin/athenaPK -i inputs/sod.in
```

then run the matched Chaa configurations listed in `tests/cases.conf`
(the `ref-*` entries) and compare with the readers in
`tests/validate/validate.py`. The frozen profiles record the generating
configurations in their headers.
