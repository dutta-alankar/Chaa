#!/bin/bash -l
# tools/slurm/freya-bench-multi.slurm — across-node scaling with the
# multi-locale (gasnet/ibv) binary: one locale per node, 40 threads per
# locale.  Submit with the maximum node count; every power-of-two
# locale count is measured inside the allocation.
#
# GASNet is bootstrapped by srun via PMIx (the compiled-in pmi
# spawner), launching chaa_real directly — no Chapel launcher, no
# node-to-node ssh needed.  CHPL_RT_MAX_HEAP_SIZE caps the registered
# ibv segment below the job's cgroup memory limit.
#
#   sbatch -N 8 tools/slurm/freya-bench-multi.slurm
#SBATCH -p p.24h
#SBATCH --ntasks-per-node=1
#SBATCH -c 40
#SBATCH --mem=180G
#SBATCH -t 04:00:00
#SBATCH -J chaa-bench-multi
#SBATCH -o chaa-bench-multi-%j.out
set -e
cd "$SLURM_SUBMIT_DIR"
source tools/slurm/freya-env.sh
export CHPL_RT_NUM_THREADS_PER_LOCALE=40
export CHPL_RT_MAX_HEAP_SIZE=16g
BIN=${CHAA_BIN:-$PWD/build-gasnet/bin/chaa_real}
STEPS=${STEPS:-50}

run() { # run <nl> <nx1> <nx2> <nx3>
  local rate
  rate=$(srun --mpi=pmix --ntasks-per-node=1 -c 40 -n $1 "$BIN" -nl $1 \
      --problem=sedov --geometry=cartesian \
      --nx1=$2 --nx2=$3 --nx3=$4 \
      --x1min=-1.2 --x1max=1.2 --x2min=-1.2 --x2max=1.2 \
      --x3min=-1.2 --x3max=1.2 --sedovR0=0.12 \
      --tstop=1e9 --maxSteps=$STEPS --logEvery=1000000 \
      --outFormats=none --outDir=/tmp/chaa-bench \
      2>&1 | sed -nE 's/.*\(([0-9.]+) Mcell-updates.*/\1/p')
  printf "MULTI %2d locales  %4dx%4dx%4d  %s Mcell/s\n" "$1" "$2" "$3" "$4" "$rate"
}

echo "== strong scaling: fixed 256^3, 1..$SLURM_NNODES locales"
nl=1
while [ $nl -le "$SLURM_NNODES" ]; do
  run $nl 256 256 256
  nl=$((nl * 2))
done

echo "== weak scaling: 256^3 cells per locale"
run 1 256 256 256
[ "$SLURM_NNODES" -ge 2 ] && run 2 512 256 256
[ "$SLURM_NNODES" -ge 4 ] && run 4 512 512 256
[ "$SLURM_NNODES" -ge 8 ] && run 8 512 512 512
echo "MULTI-BENCH-DONE"
