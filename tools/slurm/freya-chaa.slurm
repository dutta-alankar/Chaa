#!/bin/bash -l
# tools/slurm/freya-chaa.slurm — run a multi-locale Chaa simulation on
# MPCDF Freya.  One Chapel locale per node, all 40 cores per locale;
# GASNet is bootstrapped by srun via PMIx (chaa_real is launched
# directly — no Chapel launcher, no node-to-node ssh required).
#
#   sbatch -N 4 tools/slurm/freya-chaa.slurm --problem=sedov \
#          --geometry=cartesian --nx1=256 --nx2=256 --nx3=256 ...
#
# Everything after the script name is passed to chaa verbatim.
#SBATCH -p p.24h
#SBATCH --ntasks-per-node=1
#SBATCH -c 40
#SBATCH --mem=180G
#SBATCH -t 04:00:00
#SBATCH -J chaa
#SBATCH -o chaa-%j.out
set -e
cd "$SLURM_SUBMIT_DIR"
source tools/slurm/freya-env.sh
export CHPL_RT_NUM_THREADS_PER_LOCALE=40
# cap the registered ibv segment below the job's memory limit
export CHPL_RT_MAX_HEAP_SIZE=${CHPL_RT_MAX_HEAP_SIZE:-64g}

BIN=${CHAA_BIN:-$PWD/build-gasnet/bin/chaa_real}
echo "== $SLURM_NNODES locale(s) on: $SLURM_JOB_NODELIST"
srun --mpi=pmix --ntasks-per-node=1 -c 40 -n "$SLURM_NNODES" \
     "$BIN" -nl "$SLURM_NNODES" "$@"
