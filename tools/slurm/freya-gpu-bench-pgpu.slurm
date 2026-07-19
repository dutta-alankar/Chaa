#!/bin/bash -l
# freya-gpu-bench-pgpu.slurm — GPU benchmarks on Freya's p.gpu
# partition (2x V100 per node; request them explicitly with the gres
# below).  P100s are NOT supported: Chapel passes kernel arguments by
# value and Chaa's physics kernels exceed the 4 KB kernel-parameter
# limit of pre-Volta GPUs (sm_70+ allows 32 KB) — ptxas rejects the
# sm_60 build.
#
#   sbatch tools/slurm/freya-gpu-bench-pgpu.slurm

#SBATCH -J chaa-gpu-pgpu
#SBATCH -o chaa-gpu-pgpu-%j.out
#SBATCH -e chaa-gpu-pgpu-%j.err
#SBATCH --partition=p.gpu
#SBATCH --constraint="gpu"
#SBATCH --nodes=1
#SBATCH --gres=gpu:v100:2    # V100 only (see above)
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
#SBATCH --time=02:00:00

set -e
cd ~/ptmp/Chaa
source tools/slurm/freya-gpu-env.sh

gpuname=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
case "$gpuname" in
  *V100*) BIN=build-gpu-v100/bin/chaa ;;
  *P100*) echo "P100 (sm_60) is not supported — pre-Volta 4 KB kernel-"
          echo "parameter limit; request V100s (--gres=gpu:v100:2)"
          exit 1 ;;
  *)      echo "unexpected GPU: $gpuname"; exit 1 ;;
esac
echo "GPU: $gpuname -> $BIN"

run_bench() { # run_bench <label> <n> <gpus>
  local label="$1" n="$2" gpus="$3"
  local out
  out=$(CUDA_VISIBLE_DEVICES=$gpus $BIN --problem=sedov \
        --geometry=cartesian --nx1=$n --nx2=$n --nx3=$n \
        --x1min=-0.5 --x1max=0.5 --x2min=-0.5 --x2max=0.5 \
        --x3min=-0.5 --x3max=0.5 --tstop=0.02 --outFormats=txt \
        --outDir=bench-pgpu-out --logEvery=1000000000 2>&1) || true
  rm -rf bench-pgpu-out
  local rate
  rate=$(echo "$out" | grep -oE '\([0-9.]+ Mcell-updates/s\)' \
         | tr -d '( )' | cut -d' ' -f1 | sed 's/Mcell-updates\/s//' || true)
  if [ -z "$rate" ]; then
    echo "RUN-FAILED (no rate); last output lines:"
    echo "$out" | tail -8
    return 0
  fi
  echo "GPUBENCH $label n=$n gpus=$(echo $gpus | tr ',' ' ' | wc -w) rate=$rate gpu=$gpuname"
  echo "$out" | tail -2
}

echo "== size scan, 1 GPU =="
for n in 128 192 256; do
  run_bench size $n 0
done

echo "== 2-GPU strong scaling at 256^3 =="
run_bench strong 256 0
run_bench strong 256 0,1

echo "== quick validation: GPU-vs-CPU on two cases =="
CHAA_BIN=build/bin/chaa CHAA_GPU_BIN=$BIN PY=python \
  tests/run_gpu_compare.sh 2>&1 | tail -3 || echo "COMPARE FAILED"
