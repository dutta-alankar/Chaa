#!/bin/bash -l
# freya-gpu-multi.slurm — multi-node GPU benchmark on Freya:
# 2 nodes x 2 A100 = 4 GPUs, ONE Chapel locale per node driving two
# GPUs each (the same srun/PMIx pattern as the CPU multi-locale runs).
#
# Why not one co-locale per GPU: launching 4 gasnet+GPU tasks per node
# trips slurmstepd OOM/launch failures on Freya, and a single locale
# driving all 4 GPUs of a node hits a Chapel 2.8 runtime fault (CUDA
# 400 "invalid resource handle" — 2 concurrent GPUs per process are
# fine).  One locale x 2 GPUs per node avoids both.
#
# Needs the gasnet GPU binary:
#   source tools/slurm/freya-gpu-env.sh
#   export CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv
#   cmake -B build-gpu-gasnet -DCHAA_GPU=ON && cmake --build build-gpu-gasnet
#
#   cd ~/ptmp/Chaa && sbatch tools/slurm/freya-gpu-multi.slurm

#SBATCH -J chaa-gpu-multi
#SBATCH -o chaa-gpu-multi-%j.out
#SBATCH -e chaa-gpu-multi-%j.err
#SBATCH --partition=p.gpu.ampere
#SBATCH --constraint="gpu"
#SBATCH --nodes=2
#SBATCH --gres=gpu:a100:2    # 2 GPUs per node, one locale per node
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --time=01:00:00

set -e
cd ~/ptmp/Chaa
source tools/slurm/freya-gpu-env.sh
export CHPL_RT_MAX_HEAP_SIZE=64g
# the GPU gasnet runtime probes max pinnable memory by pinning until
# failure — inside the job cgroup that is an instant OOM kill; cap it
export GASNET_PHYSMEM_MAX="96 GB"
export GASNET_PHYSMEM_PROBE=0

BIN=build-gpu-gasnet/bin/chaa_real
NL=$SLURM_NNODES
NGPU=$((2*NL))

run() { # run <label> <nx1> <n23>
  local label="$1" nx1="$2" n23="$3"
  local out
  # pre-create the output dir: Lustre metadata caching can make the
  # binary's exists() see a directory a previous job removed, so its
  # own mkdir is skipped and the restart write fails
  mkdir -p bench-gpu-multi-out
  out=$(srun --mpi=pmix --ntasks-per-node=1 -c 24 -n $NL \
        $BIN -nl $NL --problem=sedov --geometry=cartesian \
        --nx1=$nx1 --nx2=$n23 --nx3=$n23 \
        --x1min=-0.5 --x1max=0.5 --x2min=-0.5 --x2max=0.5 \
        --x3min=-0.5 --x3max=0.5 --tstop=0.02 --outFormats=txt \
        --outDir=bench-gpu-multi-out --logEvery=1000000000 2>&1) || true
  rm -rf bench-gpu-multi-out
  local rate
  rate=$(echo "$out" | grep -oE '\([0-9.]+ Mcell-updates/s\)' \
         | tr -d '( )' | cut -d' ' -f1 | sed 's/Mcell-updates\/s//' || true)
  if [ -z "$rate" ]; then
    echo "RUN-FAILED (no rate); last output lines:"
    echo "$out" | tail -25
    return 0
  fi
  echo "GPUBENCH $label n=$nx1 gpus=$NGPU rate=$rate"
  echo "$out" | tail -3
}

echo "== strong scaling point: 256^3 total on $NGPU GPUs ($NL nodes) =="
run strong 256 256

echo "== weak scaling point: 256^3 per GPU on $NGPU GPUs ($NL nodes) =="
run weak $((256*NGPU)) 256
