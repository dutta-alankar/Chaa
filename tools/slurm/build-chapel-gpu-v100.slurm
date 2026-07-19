#!/bin/bash -l
# one-off: separate GPU Chapel tree for the V100 nodes.  The GPU
# runtime bakes device code (reduction helpers etc.) for CHPL_GPU_ARCH
# at runtime-build time, so an sm_80 runtime cannot serve sm_70
# binaries — V100 needs its own tree.
#SBATCH -p p.24h
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH --mem=90G
#SBATCH -t 04:00:00
#SBATCH -J chpl-gpu-v100
#SBATCH -o chpl-gpu-v100-%j.out
set -e
module purge
module load gcc/14 cuda/12.8 cmake/3.28
mkdir -p ~/ptmp/tmp && export TMPDIR=$HOME/ptmp/tmp
mkdir -p ~/ptmp/Chaa/chapel-gpu-v100
cd ~/ptmp/Chaa/chapel-gpu-v100
if [ ! -d chapel-2.8.0 ]; then
  tar xzf ~/ptmp/Chaa/chapel-gpu/chapel-2.8.0.tar.gz -C . 2>/dev/null \
    || tar xzf ~/ptmp/Chaa/chapel-2.8.0.tar.gz -C .
fi
cd chapel-2.8.0
export CHPL_HOME=$PWD
export CHPL_LLVM=bundled
export CHPL_LOCALE_MODEL=gpu
export CHPL_GPU=nvidia
export CHPL_CUDA_PATH="$CUDA_HOME"
export CHPL_GPU_ARCH=sm_70
make -j32
echo "BUILD-V100CHAPEL-OK"
