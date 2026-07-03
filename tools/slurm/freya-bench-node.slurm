#!/bin/bash -l
# tools/slurm/freya-bench-node.slurm — within-node thread scaling
# (strong: fixed 128^3; weak: 64x64x64 cells per 5 threads) with the
# single-locale (CHPL_COMM=none) binary on one Freya node.
#
#   sbatch tools/slurm/freya-bench-node.slurm
#SBATCH -p p.24h
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH --mem=100G
#SBATCH -t 03:00:00
#SBATCH -J chaa-bench-node
#SBATCH -o chaa-bench-node-%j.out
set -e
cd "$SLURM_SUBMIT_DIR"
source tools/slurm/freya-env.sh
BIN=${CHAA_BIN:-$PWD/build/bin/chaa}
STEPS=${STEPS:-50}

run() { # run <threads> <nx1> <nx2> <nx3>
  local rate
  rate=$(CHPL_RT_NUM_THREADS_PER_LOCALE=$1 "$BIN" \
      --problem=sedov --geometry=cartesian \
      --nx1=$2 --nx2=$3 --nx3=$4 \
      --x1min=-1.2 --x1max=1.2 --x2min=-1.2 --x2max=1.2 \
      --x3min=-1.2 --x3max=1.2 --sedovR0=0.12 \
      --tstop=1e9 --maxSteps=$STEPS --logEvery=1000000 \
      --outFormats=none 2>&1 | sed -nE 's/.*\(([0-9.]+) Mcell-updates.*/\1/p')
  printf "NODE %2d threads  %4dx%4dx%4d  %s Mcell/s\n" "$1" "$2" "$3" "$4" "$rate"
}

echo "== strong scaling: fixed 128^3, 1..40 threads"
for t in 1 2 5 10 20 40; do run $t 128 128 128; done

echo "== weak scaling: 64^3 cells per 5 threads"
run 5  64  64  64
run 10 128 64  64
run 20 128 128 64
run 40 128 128 128
echo "NODE-BENCH-DONE"
