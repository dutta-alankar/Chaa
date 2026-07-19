#!/bin/bash -l
# freya-gpu-bench.slurm — GPU benchmarks on a Freya A100 node:
# 3D Sedov (no output during the loop), evolution-loop Mcell/s.
#
#   - single-GPU throughput vs problem size (128^3 .. 384^3)
#   - single-node GPU scaling: 1/2/4 A100s (one process, blocks split
#     over the visible devices via CUDA_VISIBLE_DEVICES)
#   - CPU-node reference at 256^3 for the speedup figure
#
#   cd ~/ptmp/Chaa && sbatch tools/slurm/freya-gpu-bench.slurm

#SBATCH -J chaa-gpu-bench
#SBATCH -o chaa-gpu-bench.%j.out
#SBATCH -e chaa-gpu-bench.%j.err
#SBATCH --partition=p.gpu.ampere
#SBATCH --constraint="gpu"
#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH --time=03:00:00

set -e
cd ~/ptmp/Chaa
source tools/slurm/freya-gpu-env.sh

BIN=build-gpu/bin/chaa
run_bench() { # run_bench <label> <n> <gpus>
  local label="$1" n="$2" gpus="$3"
  local out
  out=$(CUDA_VISIBLE_DEVICES=$gpus $BIN --problem=sedov \
        --geometry=cartesian --nx1=$n --nx2=$n --nx3=$n \
        --x1min=-0.5 --x1max=0.5 --x2min=-0.5 --x2max=0.5 \
        --x3min=-0.5 --x3max=0.5 --tstop=0.02 --outFormats=txt \
        --outDir=bench-gpu-out --logEvery=1000000000 2>&1) || true
  rm -rf bench-gpu-out
  local rate
  rate=$(echo "$out" | grep -oE '\([0-9.]+ Mcell-updates/s\)' \
         | tr -d '( )' | cut -d' ' -f1 | sed 's/Mcell-updates\/s//' || true)
  if [ -z "$rate" ]; then
    echo "RUN-FAILED (no rate); last output lines:"
    echo "$out" | tail -8
    return 0
  fi
  echo "GPUBENCH $label n=$n gpus=$(echo $gpus | tr ',' ' ' | wc -w) rate=$rate"
  echo "$out" | tail -2
}

echo "== size scan, 1 GPU =="
for n in 128 192 256 320 384; do
  run_bench size $n 0
done

echo "== GPU scaling at 256^3 (strong) =="
run_bench strong 256 0
run_bench strong 256 0,1
run_bench strong 256 0,1,2,3

echo "== weak scaling, 256^3 per GPU (x1 extent grows) =="
weak() { # weak <nx1> <gpus>
  local nx1="$1" gpus="$2"
  local out
  out=$(CUDA_VISIBLE_DEVICES=$gpus $BIN --problem=sedov \
        --geometry=cartesian --nx1=$nx1 --nx2=256 --nx3=256 \
        --x1min=-0.5 --x1max=0.5 --x2min=-0.5 --x2max=0.5 \
        --x3min=-0.5 --x3max=0.5 --tstop=0.02 --outFormats=txt \
        --outDir=bench-gpu-out --logEvery=1000000000 2>&1) || true
  rm -rf bench-gpu-out
  local rate
  rate=$(echo "$out" | grep -oE '\([0-9.]+ Mcell-updates/s\)' \
         | tr -d '( )' | cut -d' ' -f1 | sed 's/Mcell-updates\/s//' || true)
  if [ -z "$rate" ]; then
    echo "RUN-FAILED (no rate); last output lines:"
    echo "$out" | tail -8
    return 0
  fi
  echo "GPUBENCH weak n=$nx1 gpus=$(echo $gpus | tr ',' ' ' | wc -w) rate=$rate"
  echo "$out" | tail -2
}
weak 256 0
weak 512 0,1
weak 1024 0,1,2,3
