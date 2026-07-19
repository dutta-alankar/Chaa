#!/bin/bash -l
# tools/slurm/build-chapel-gpu.slurm — build a GPU-enabled Chapel 2.8
# on a Freya compute node.
#
# This is a second, separate Chapel installation: GPU code generation
# needs Chapel's LLVM backend (CHPL_LLVM=bundled), while the CPU
# install on Freya uses the C backend (CHPL_LLVM=none).  Building the
# bundled LLVM takes ~1.5-2 h on 40 cores.  No GPU is needed at build
# time — only the CUDA toolkit module.
#
#   mkdir -p ~/ptmp/Chaa/chapel-gpu && cd ~/ptmp/Chaa/chapel-gpu
#   wget https://github.com/chapel-lang/chapel/releases/download/2.8.0/chapel-2.8.0.tar.gz
#   tar xzf chapel-2.8.0.tar.gz
#   cd ~/ptmp/Chaa && sbatch tools/slurm/build-chapel-gpu.slurm
#SBATCH -p p.24h
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH --mem=90G
#SBATCH -t 05:00:00
#SBATCH -J chpl-gpu-build
#SBATCH -o chpl-gpu-build-%j.out
set -e
module purge
module load gcc/14 cuda/12.8 cmake/3.28
cd ~/ptmp/Chaa/chapel-gpu/chapel-2.8.0
export CHPL_HOME=$PWD
export CHPL_LLVM=bundled
export CHPL_LOCALE_MODEL=gpu
export CHPL_GPU=nvidia
export CHPL_CUDA_PATH="$CUDA_HOME"
export CHPL_GPU_ARCH=sm_80          # A100; rebuild the runtime with
                                    # sm_70/sm_60 for V100/P100 (p.gpu)
echo "=== compiler (bundled LLVM) + comm=none GPU runtime ==="
make -j32
echo "=== multi-node GPU runtime (GASNet over ibv) ==="
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv \
  CHPL_LAUNCHER=slurm-gasnetrun_ibv make -j32 runtime
echo "BUILD-DONE"
bin/linux64-x86_64/chpl --version | head -2
