#!/bin/bash -l
# tools/slurm/freya-tests.slurm — run the full validated Chaa test
# suite on one Freya compute node (single locale, all cores).
#
#   sbatch tools/slurm/freya-tests.slurm
#SBATCH -p p.24h
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH --mem=60G
#SBATCH -t 02:00:00
#SBATCH -J chaa-tests
#SBATCH -o chaa-tests-%j.out
set -e
cd "$SLURM_SUBMIT_DIR"
source tools/slurm/freya-env.sh
export CHAA_BIN=${CHAA_BIN:-$PWD/build/bin/chaa}
PY=$(which python3) tests/run_case.sh all
