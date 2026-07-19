#!/bin/bash -l
# freya-gpu-tests.slurm — validate the GPU build on a Freya A100 node:
#   1. the GPU-vs-CPU comparison matrix (tests/run_gpu_compare.sh),
#   2. the full validated test suite driven with the GPU binary.
#
#   cd ~/ptmp/Chaa && sbatch tools/slurm/freya-gpu-tests.slurm
#
# Expects build/bin/chaa (CPU) and build-gpu/bin/chaa (GPU) to exist
# (see tools/slurm/README.md).

#SBATCH -J chaa-gpu-tests
#SBATCH -o chaa-gpu-tests.%j.out
#SBATCH -e chaa-gpu-tests.%j.err
#SBATCH --partition=p.gpu.ampere
#SBATCH --constraint="gpu"
#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=120G
#SBATCH --time=04:00:00

set -e
cd ~/ptmp/Chaa
source tools/slurm/freya-gpu-env.sh

echo "== GPU-vs-CPU comparison matrix =="
CHAA_BIN=build/bin/chaa CHAA_GPU_BIN=build-gpu/bin/chaa \
  PY=python tests/run_gpu_compare.sh

echo "== full validated suite with the GPU binary =="
pass=0; fail=0; failed=""
while IFS='|' read -r name flags; do
  name=$(echo "$name" | xargs)
  case "$name" in
    ""|\#*) continue ;;
    cylinder-flow) echo "SKIP $name (internal BCs are CPU-only)";
                   continue ;;
  esac
  if CHAA_BIN=build-gpu/bin/chaa PY=python \
       tests/run_case.sh "$name" > gpu-case-$name.log 2>&1; then
    echo "PASS $name"; pass=$((pass+1))
  else
    echo "FAIL $name (see gpu-case-$name.log)"; fail=$((fail+1))
    failed="$failed $name"
  fi
done < tests/cases.conf
echo "GPU suite: $pass passed, $fail failed${failed:+ ($failed)}"
[ "$fail" -eq 0 ]
