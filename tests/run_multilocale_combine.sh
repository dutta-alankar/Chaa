#!/usr/bin/env bash
# Multi-locale integration test: run Chaa on several locales (gasnet),
# combine the per-locale piece output into single global files with
# tools/combine_pieces.py, and require the result to match a
# single-locale reference run.
#
# Needs a single-locale binary (CHAA_BIN, default build/bin/chaa) and a
# multi-locale one (CHAA_ML_BIN, default build-gasnet/bin/chaa), plus
# python3 with numpy + h5py (PY).
#   tests/run_multilocale_combine.sh [num-locales]
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
NL=${1:-3}

# Chapel refuses more co-locales than physical cores, so clamp NL to
# what the machine actually has (CI runners can be as small as 2 cores).
if command -v lscpu > /dev/null 2>&1; then
  CORES=$(lscpu -p=Core,Socket 2>/dev/null | grep -v '^#' | sort -u | wc -l)
elif command -v sysctl > /dev/null 2>&1; then
  CORES=$(sysctl -n hw.physicalcpu 2>/dev/null || echo 0)
else
  CORES=0
fi
if [ "${CORES:-0}" -ge 2 ] && [ "$NL" -gt "$CORES" ]; then
  echo "note: only $CORES physical cores available; using $CORES locales instead of $NL"
  NL=$CORES
fi

BIN=${CHAA_BIN:-build/bin/chaa}
ML=${CHAA_ML_BIN:-build-gasnet/bin/chaa}
: "${PY:=python3}"
OUT=${TEST_OUT:-test-output}

run_pair() { # run_pair <name> <flags...>
  local name="$1"; shift
  rm -rf "$OUT/mlcombine-$name" "$OUT/mlcombine-$name-ref"
  echo "=== $name: $NL-locale run ==="
  "$ML" -nl "$NL" "$@" --outDir="$OUT/mlcombine-$name" \
        --logEvery=1000000000
  echo "=== $name: single-locale reference ==="
  "$BIN" "$@" --outDir="$OUT/mlcombine-$name-ref" --logEvery=1000000000
  echo "=== $name: combine pieces ==="
  "$PY" tools/combine_pieces.py "$OUT/mlcombine-$name" --clean
  echo "=== $name: compare with reference ==="
  "$PY" tests/validate/combine_check.py \
        "$OUT/mlcombine-$name" "$OUT/mlcombine-$name-ref"
}

# Cartesian (rectilinear VTK) with a tracer field
run_pair kh --problem=kh --geometry=cartesian --nx1=48 --nx2=48 \
  --x1min=0 --x1max=1 --x2min=0 --x2max=1 \
  --bcX1min=periodic --bcX1max=periodic \
  --bcX2min=periodic --bcX2max=periodic \
  --tstop=0.2 --outFormats=hdf5,vtk

# curvilinear polar (mapped structured-grid VTK)
run_pair blast --problem=blast --geometry=polar --nx1=48 --nx2=96 \
  --x1min=0.5 --x1max=1.5 --x2min=0 --x2max=6.28318530717959 \
  --bcX2min=periodic --bcX2max=periodic \
  --cen1=-1 --cen2=0 --blastR0=0.1 --tstop=0.05 --outFormats=hdf5,vtk

echo "multi-locale combine test: all comparisons passed"
