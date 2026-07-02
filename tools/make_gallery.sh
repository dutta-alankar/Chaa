#!/usr/bin/env bash
# tools/make_gallery.sh — regenerate every figure committed under
# docs/assets/plots/ from scratch: runs the test cases, then produces
# the analytic-comparison and field/slice plots with the bundled tools.
#
#   PY=/path/to/python tools/make_gallery.sh
#
# (python needs numpy, h5py, matplotlib; the chaa binary must be built.)
set -e
cd "$(dirname "$0")/.."
PY=${PY:-python3}
OUT=docs/assets/plots
mkdir -p "$OUT"

CASES="sod-1d-cart sod-1d-iso sedov-1d-sph sedov-2d-cyl sedov-3d-cart
       sedov-3d-sph blast-2d-polar blast-3d-polar riemann2d dmr kh rt
       vortex vortex-particles-ring taylor-couette cylinder-flow
       twoblast-1d thermal-diffusion cooling-box linear-wave
       turbulence-2d"
for c in $CASES; do
  echo "== running $c"
  PY=$PY tests/run_case.sh "$c" > /dev/null
done

# the epicycle oscillation needs a run with frequent dumps
echo "== running epicycle (frequent dumps)"
./build/bin/chaa --problem=epicycle --nx1=32 --nx2=32 \
  --x1min=-0.5 --x1max=0.5 --x2min=-0.5 --x2max=0.5 \
  --bcX1min=shear-periodic --bcX1max=shear-periodic \
  --bcX2min=periodic --bcX2max=periodic \
  --omegaRot=1.0 --shearQ=1.5 --waveAmp=0.01 --inPrs=10 \
  --tstop=12.566 --outDt=0.3 --outFormats=hdf5 \
  --outDir=test-output/epicycle-dumps > /dev/null

T=test-output
echo "== analytic comparisons"
$PY tools/plot_compare.py sod            $T/sod-1d-cart       --save $OUT/sod-vs-exact.png
$PY tools/plot_compare.py sod-iso        $T/sod-1d-iso        --save $OUT/sodiso-vs-exact.png
$PY tools/plot_compare.py sedov          $T/sedov-1d-sph      --save $OUT/sedov1d-radius.png
$PY tools/plot_compare.py sedov          $T/sedov-3d-cart     --save $OUT/sedov3d-radius.png
$PY tools/plot_compare.py taylor-couette $T/taylor-couette    --save $OUT/taylor-couette-profile.png
$PY tools/plot_compare.py thermal-wave   $T/thermal-diffusion --save $OUT/thermal-decay.png
$PY tools/plot_compare.py cooling        $T/cooling-box       --save $OUT/cooling-townsend.png
$PY tools/plot_compare.py linear-wave    $T/linear-wave       --save $OUT/linear-wave-return.png
$PY tools/plot_compare.py vortex         $T/vortex            --save $OUT/vortex-return.png
$PY tools/plot_compare.py epicycle       $T/epicycle-dumps    --save $OUT/epicycle-oscillation.png

echo "== field maps & slices"
$PY tools/plot_fields.py $T/sedov-2d-cyl          --log rho        --save $OUT/sedov2d-cyl-fields.png
$PY tools/plot_fields.py $T/blast-2d-polar        --fields rho,prs --save $OUT/blast-polar-fields.png
$PY tools/plot_fields.py $T/dmr                   --fields rho     --save $OUT/dmr-fields.png
$PY tools/plot_fields.py $T/kh                    --fields rho,sc0 --save $OUT/kh-fields.png
$PY tools/plot_fields.py $T/riemann2d             --fields rho     --save $OUT/riemann2d-fields.png
$PY tools/plot_fields.py $T/rt                    --fields rho     --save $OUT/rt-fields.png
$PY tools/plot_fields.py $T/cylinder-flow         --fields vx1     --save $OUT/cylinder-wake.png
$PY tools/plot_fields.py $T/twoblast-1d           --fields rho,prs --save $OUT/twoblast-fields.png
$PY tools/plot_fields.py $T/turbulence-2d         --fields sc0,vx1 --save $OUT/turbulence-fields.png
$PY tools/plot_fields.py $T/vortex-particles-ring --fields rho     --save $OUT/ring-particles.png
# 3D runs: slice plots (mid-plane and an off-centre cut)
$PY tools/plot_fields.py $T/sedov-3d-cart  --fields rho --slice x3,0.5  --save $OUT/sedov3d-slice.png
$PY tools/plot_fields.py $T/sedov-3d-cart  --fields rho --slice x1,0.75 --save $OUT/sedov3d-slice-x1.png
$PY tools/plot_fields.py $T/blast-3d-polar --fields prs --slice x3,0.5  --save $OUT/blast3d-polar-slice.png
$PY tools/plot_fields.py $T/sedov-3d-sph   --fields rho --slice x3,0.5  --save $OUT/sedov3d-sph-slice.png

echo "gallery written to $OUT"
