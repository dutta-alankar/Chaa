#!/bin/bash -l
# tools/slurm/build-chapel.slurm — build Chapel 2.8 on a Freya compute
# node: the single-locale (CHPL_COMM=none) compiler+runtime plus the
# multi-locale GASNet/ibv runtime used by Chaa.
#
#   sbatch tools/slurm/build-chapel.slurm
#
# Expects the Chapel source tree at ~/ptmp/Chaa/chapel-2.8.0 (see
# tools/slurm/README.md).
#SBATCH -p p.24h
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH --mem=60G
#SBATCH -t 02:00:00
#SBATCH -J chpl-build
#SBATCH -o chpl-build-%j.out
set -e
module purge
module load gcc/14 cmake/3.28
cd ~/ptmp/Chaa/chapel-2.8.0
export CHPL_HOME=$PWD
export CHPL_LLVM=none
echo "=== base build (single-locale: compiler + comm=none runtime) ==="
make -j32
echo "=== multi-locale runtime (GASNet over ibv + slurm launcher) ==="
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv \
  CHPL_LAUNCHER=slurm-gasnetrun_ibv make -j32 runtime
echo "BUILD-DONE"
