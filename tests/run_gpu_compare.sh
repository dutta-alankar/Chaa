#!/usr/bin/env bash
# GPU integration test: run a matrix of small configurations twice —
# once with a GPU-enabled binary (CHAA_GPU_BIN) and once with the
# reference CPU binary (CHAA_BIN) — and require every field dump (and
# particle dump) to agree to round-off.  The matrix covers every code
# path that behaves differently in a GPU build: all reconstructions/
# integrators/Riemann solvers, curvilinear geometry, viscosity,
# conduction, cooling, isothermal EOS, user-defined BCs, OU forcing,
# tracer particles, FARGO, self-gravity and the restart machinery.
#
#   CHAA_BIN=build/bin/chaa CHAA_GPU_BIN=build-gpu/bin/chaa \
#     tests/run_gpu_compare.sh
#
# Works with the cpu-as-device runtime (CHPL_GPU=cpu) for functional
# checking and with real GPUs.  PY needs numpy + h5py.
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
BIN=${CHAA_BIN:-build/bin/chaa}
GBIN=${CHAA_GPU_BIN:-build-gpu/bin/chaa}
: "${PY:=python3}"
OUT=${TEST_OUT:-test-output}

fail=0

run_pair() { # run_pair <name> <flags...>
  local name="$1"; shift
  rm -rf "$OUT/gpucmp-$name" "$OUT/gpucmp-$name-ref"
  echo "=== $name: GPU run ==="
  "$GBIN" "$@" --outDir="$OUT/gpucmp-$name" --logEvery=1000000000 \
    | grep -E "gpu|done" || true
  echo "=== $name: CPU reference ==="
  "$BIN" "$@" --outDir="$OUT/gpucmp-$name-ref" --logEvery=1000000000 \
    > /dev/null
  echo "=== $name: compare ==="
  "$PY" tests/validate/gpu_check.py \
        "$OUT/gpucmp-$name" "$OUT/gpucmp-$name-ref" || fail=1
}

# 1. shock tube, default hllc/linear/rk2 (1D Cartesian)
run_pair sod --problem=sod --nx1=128 --tstop=0.2 --outFormats=txt

# 2. exact Riemann solver + PPM + RK3
run_pair sod-exact-ppm --problem=sod --riemann=exact --recon=ppm \
  --integrator=rk3 --nx1=96 --tstop=0.15 --outFormats=txt

# 3. isothermal EOS + WENO-Z + VL2 (+ llf solver)
run_pair sod-iso-wenoz --problem=sod --eos=isothermal --csIso=1.0 \
  --recon=wenoz --integrator=vl2 --riemann=llf --nx1=96 --tstop=0.1 \
  --outFormats=txt

# 4. limo3 + hll + a stretched grid
run_pair sod-limo3 --problem=sod --recon=limo3 --riemann=hll \
  --gridX1=stretch --stretchX1=1.02 --nx1=96 --tstop=0.15 \
  --outFormats=txt

# 5. 1D spherical Sedov (geometric source terms, reflect BC)
run_pair sedov-sph --problem=sedov --geometry=spherical --nx1=128 \
  --x1min=0 --x1max=1 --bcX1min=reflect --tstop=0.05 --outFormats=txt

# 6. thermal conduction + optically thin cooling
run_pair cond-cool --problem=sod --nx1=96 --kappa=0.005 \
  --coolLambda0=0.5 --coolAlpha=0.5 --tstop=0.1 --outFormats=txt

# 7. 2D Kelvin-Helmholtz with tracers + full viscous stress + scalar
#    diffusion (Cartesian, periodic)
run_pair kh-visc --problem=kh --geometry=cartesian --nx1=48 --nx2=48 \
  --x1min=0 --x1max=1 --x2min=0 --x2max=1 \
  --bcX1min=periodic --bcX1max=periodic \
  --bcX2min=periodic --bcX2max=periodic \
  --mu=0.001 --scDiff=0.001 --tstop=0.1 --outFormats=hdf5

# 8. 2D polar blast (curvilinear geometry + metric factors)
run_pair blast-polar --problem=blast --geometry=polar --nx1=48 --nx2=96 \
  --x1min=0.5 --x1max=1.5 --x2min=0 --x2max=6.28318530717959 \
  --bcX2min=periodic --bcX2max=periodic --cen1=-1 --cen2=0 \
  --blastR0=0.1 --tstop=0.05 --outFormats=hdf5

# 9. Taylor-Couette (user-defined BCs + cylindrical viscous source)
run_pair taylor-couette --problem=taylorCouette --geometry=cylindrical \
  --nx1=48 --x1min=1.0 --x1max=2.0 --bcX1min=userdef --bcX1max=userdef \
  --mu=0.05 --inRho=1.0 --inPrs=10.0 --tcOmegaIn=1.0 --tcOmegaOut=0.0 \
  --tstop=0.5 --outFormats=txt

# 10. double Mach reflection (time-dependent user BC, inflow BC)
run_pair dmr --problem=dmr --geometry=cartesian --nx1=96 --nx2=24 \
  --x1min=0 --x1max=4 --x2min=0 --x2max=1 --bcX1min=inflow \
  --bcX2min=userdef --bcX2max=userdef --inRho=8.0 \
  --inVx1=7.144709581221619 --inVx2=-4.125 --inPrs=116.5 \
  --tstop=0.05 --outFormats=hdf5

# 11. OU-driven turbulence (forcing tables on the device)
run_pair turbulence --problem=turbulence --geometry=cartesian \
  --nx1=32 --nx2=32 --x1min=0 --x1max=1 --x2min=0 --x2max=1 \
  --bcX1min=periodic --bcX1max=periodic \
  --bcX2min=periodic --bcX2max=periodic \
  --eos=isothermal --csIso=1 --forceAmp=2.0 --forceTcorr=0.5 \
  --forceKmin=1 --forceKmax=2 --forceSeed=1234 --tstop=0.2 \
  --outFormats=hdf5

# 12. Lagrangian tracer particles (host advection off device fields)
run_pair particles --problem=vortex --geometry=cartesian --nx1=48 \
  --nx2=48 --x1min=0 --x1max=10 --x2min=0 --x2max=10 --cen1=5 --cen2=5 \
  --bcX1min=periodic --bcX1max=periodic \
  --bcX2min=periodic --bcX2max=periodic \
  --nParticles=64 --partSeed=7 --tstop=1.0 --outFormats=hdf5

# 13. shearing box + FARGO (shear-periodic BCs, host-staged remap)
run_pair epicycle-fargo --problem=epicycle --nx1=32 --nx2=32 \
  --x1min=-0.5 --x1max=0.5 --x2min=-0.5 --x2max=0.5 \
  --bcX1min=shear-periodic --bcX1max=shear-periodic \
  --bcX2min=periodic --bcX2max=periodic --omegaRot=1.0 --shearQ=1.5 \
  --waveAmp=0.01 --inPrs=10 --fargo=on --tstop=0.5 --outFormats=hdf5

# 14. self-gravity (host CG solve + device potential)
run_pair selfgrav --problem=thermalWave --nx1=64 --x1min=-0.5 \
  --x1max=0.5 --bcX1min=periodic --bcX1max=periodic --twAmp=0.01 \
  --sgFourPiG=1.0 --dtMax=1e-3 --tstop=0.01 --outFormats=txt

# 15. restart machinery in the GPU build: a run checkpointed mid-way and
#     resumed must reproduce the uninterrupted GPU run exactly
name=restart
rm -rf "$OUT/gpucmp-$name" "$OUT/gpucmp-$name-ref"
echo "=== $name: uninterrupted GPU run ==="
"$GBIN" --problem=sod --nx1=128 --tstop=0.2 --outDt=0.05 \
  --outFormats=txt --outDir="$OUT/gpucmp-$name-ref" \
  --logEvery=1000000000 > /dev/null
echo "=== $name: stopped + resumed GPU run ==="
mkdir -p "$OUT/gpucmp-$name"
touch "$OUT/gpucmp-$name/stop"      # graceful stop after the 1st step
"$GBIN" --problem=sod --nx1=128 --tstop=0.2 --outDt=0.05 \
  --outFormats=txt --outDir="$OUT/gpucmp-$name" \
  --logEvery=1000000000 > /dev/null
"$GBIN" --problem=sod --nx1=128 --tstop=0.2 --outDt=0.05 \
  --outFormats=txt --outDir="$OUT/gpucmp-$name" --restart=true \
  --logEvery=1000000000 > /dev/null
echo "=== $name: compare ==="
"$PY" tests/validate/gpu_check.py \
      "$OUT/gpucmp-$name" "$OUT/gpucmp-$name-ref" --exact || fail=1

if [ "$fail" -ne 0 ]; then
  echo "GPU comparison test: FAILURES (see above)"
  exit 1
fi
echo "GPU comparison test: all comparisons passed"
