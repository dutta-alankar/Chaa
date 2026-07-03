# tools/slurm/freya-env.sh — environment for building/running Chaa on
# MPCDF's Freya cluster (2x Xeon Gold 6138, 40 cores/node, Omni-Path).
# Source this from an interactive shell or a SLURM script:
#     source tools/slurm/freya-env.sh
# Adjust CHPL_HOME if your Chapel tree lives elsewhere; on other
# clusters change the modules and the GASNet substrate to match the
# fabric (ibv for InfiniBand/Omni-Path verbs, ofi/udp otherwise).
module purge
# hdf5-serial is exposed under gcc/14 in freya's module hierarchy;
# serial HDF5 is the right choice for Chaa (multi-locale output writes
# independent per-locale piece files — no MPI-IO code path)
module load gcc/14 cmake/3.28 hdf5-serial/1.14.1

export CHPL_HOME=${CHPL_HOME:-$HOME/ptmp/Chaa/chapel-2.8.0}
export CHPL_LLVM=none
export PATH=$CHPL_HOME/bin/linux64-x86_64:$PATH

# the module sets HDF5_ROOT (which Chaa's CMake honours)
export LD_LIBRARY_PATH=$HDF5_ROOT/lib:$LD_LIBRARY_PATH

# python (validation, plotting): uv-managed venv
#   uv venv ~/ptmp/Chaa/venv --python 3.12
#   uv pip install --python ~/ptmp/Chaa/venv/bin/python numpy h5py matplotlib
source $HOME/ptmp/Chaa/venv/bin/activate

# multi-locale (GASNet over the Omni-Path verbs interface + slurm):
#   export CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv
#   export CHPL_LAUNCHER=slurm-gasnetrun_ibv
# single-locale (node-local, no comm layer): leave CHPL_COMM unset.
