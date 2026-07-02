#!/usr/bin/env bash
# tools/bench.sh — Chaa scaling benchmarks (strong + weak, single- and
# multi-locale).
#
# Usage:
#   tools/bench.sh                # single-locale (thread) scaling with
#                                 # build/bin/chaa
#   tools/bench.sh build-gasnet/bin/chaa 4    # add multi-locale scaling
#                                 # with a CHPL_COMM=gasnet binary up to
#                                 # the given locale count
#
# The benchmark problem is the 3D Cartesian Sedov blast (no I/O), and
# the figure of merit is the Mcell-updates/s the driver prints.
set -u
cd "$(dirname "$0")/.."

BIN=${BIN:-build/bin/chaa}
GASNET_BIN=${1:-}
MAXNL=${2:-4}
STEPS=${STEPS:-50}
MAXTHREADS=$(getconf _NPROCESSORS_ONLN)

run() { # run <binary> <nl> <threads/locale> <nx1> <nx2> <nx3>
  local bin=$1 nl=$2 thr=$3 n1=$4 n2=$5 n3=$6 launch=()
  if [ "$nl" -gt 1 ] || [ "$bin" != "$BIN" ]; then launch=(-nl "$nl"); fi
  local out
  out=$(CHPL_RT_NUM_THREADS_PER_LOCALE=$thr "$bin" \
        ${launch[@]+"${launch[@]}"} \
        --problem=sedov --geometry=cartesian \
        --nx1="$n1" --nx2="$n2" --nx3="$n3" \
        --x1min=-1.2 --x1max=1.2 --x2min=-1.2 --x2max=1.2 \
        --x3min=-1.2 --x3max=1.2 --sedovR0=0.12 \
        --tstop=1e9 --maxSteps="$STEPS" --logEvery=1000000 \
        --outFormats=none 2>&1 | grep "Mcell-updates/s")
  local rate
  rate=$(echo "$out" | sed -E 's/.*\(([0-9.]+) Mcell-updates.*/\1/')
  printf "| %2s | %2s | %4sx%4sx%4s | %8s |\n" "$nl" "$thr" "$n1" "$n2" "$n3" "$rate"
}

echo "## single-locale strong scaling (fixed 128^3, varying threads)"
echo "| locales | threads | grid | Mcell/s |"
echo "|---|---|---|---|"
for t in 1 2 4 "$MAXTHREADS"; do
  [ "$t" -le "$MAXTHREADS" ] && run "$BIN" 1 "$t" 128 128 128
done

echo
echo "## single-locale weak scaling (64^3 cells per thread)"
echo "| locales | threads | grid | Mcell/s |"
echo "|---|---|---|---|"
run "$BIN" 1 1 64 64 64
run "$BIN" 1 2 128 64 64
run "$BIN" 1 4 128 128 64
[ "$MAXTHREADS" -ge 8 ] && run "$BIN" 1 8 128 128 128

if [ -n "$GASNET_BIN" ]; then
  echo
  echo "## multi-locale strong scaling (fixed 128^3, fixed total threads)"
  echo "| locales | threads | grid | Mcell/s |"
  echo "|---|---|---|---|"
  thr_tot=$MAXTHREADS
  nl=1
  while [ "$nl" -le "$MAXNL" ]; do
    run "$GASNET_BIN" "$nl" $((thr_tot / nl)) 128 128 128
    nl=$((nl * 2))
  done

  echo
  echo "## multi-locale weak scaling (128x64x64 cells per locale)"
  echo "| locales | threads | grid | Mcell/s |"
  echo "|---|---|---|---|"
  thr=$((MAXTHREADS / MAXNL))
  run "$GASNET_BIN" 1 "$thr" 128 64 64
  [ "$MAXNL" -ge 2 ] && run "$GASNET_BIN" 2 "$thr" 128 128 64
  [ "$MAXNL" -ge 4 ] && run "$GASNET_BIN" 4 "$thr" 128 128 128
fi
