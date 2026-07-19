# freya-gpu-env.sh — environment for building/running GPU Chaa on Freya
# (MPCDF).  Source this before configuring, compiling or running the
# GPU build:   source tools/slurm/freya-gpu-env.sh
#
# The GPU toolchain is separate from the CPU one (freya-env.sh): Chapel
# needs its LLVM backend for GPU code generation, so a second Chapel
# installation is built with the bundled LLVM (see
# build-chapel-gpu.slurm) against the CUDA module.

# (Chapel's GPU backend uses its own LLVM for device code generation,
# not nvcc, so the CUDA/gcc pairing is unconstrained: use the newest
# supported CUDA and the same gcc/hdf5 stack as the CPU build.)
module purge
module load gcc/14 cuda/12.8 cmake/3.28
module load hdf5-serial/1.14.1
# on the GPU nodes the module sets HDF5_HOME but not the runtime path
export LD_LIBRARY_PATH="$HDF5_HOME/lib:$LD_LIBRARY_PATH"

# GPU-enabled Chapel (built by build-chapel-gpu.slurm)
export CHPL_HOME=$HOME/ptmp/Chaa/chapel-gpu/chapel-2.8.0
export PATH="$CHPL_HOME/bin/linux64-x86_64:$PATH"

export CHPL_LLVM=bundled           # LLVM backend (required for GPUs)
export CHPL_LOCALE_MODEL=gpu
export CHPL_GPU=nvidia
export CHPL_CUDA_PATH="$CUDA_HOME"
export CHPL_GPU_ARCH=sm_80         # A100 (p.gpu.ampere); use sm_60 for
                                   # P100 / sm_70 for V100 (p.gpu)
export CHPL_RT_NUM_THREADS_PER_LOCALE=MAX_LOGICAL

# chpl's GPU codegen writes GB-scale temporaries; compute-node /tmp is
# RAM-backed and counted against the job cgroup (compiles hang at ~0%
# CPU when it fills) — keep compiler temporaries on Lustre
mkdir -p "$HOME/ptmp/tmp"
export TMPDIR="$HOME/ptmp/tmp"

# python venv (numpy/h5py for the validation scripts)
[ -f "$HOME/ptmp/Chaa/venv/bin/activate" ] && \
  source "$HOME/ptmp/Chaa/venv/bin/activate"
