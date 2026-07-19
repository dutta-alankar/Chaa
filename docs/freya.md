# Running on Freya (MPCDF)

A complete, tested walkthrough for MPCDF's **Freya** cluster: building
Chapel three ways — single-locale, multi-locale and GPU — building
Chaa against each, and submitting the matching SLURM jobs.  All
scripts live in [`tools/slurm/`](https://github.com/dutta-alankar/Chaa/tree/main/tools/slurm);
the working directory throughout is `~/ptmp/Chaa` (a clone of the
repository on the large ptmp filesystem).

Freya's relevant hardware:

| partition | nodes | per node |
|---|---|---|
| `p.24h` / `p.test` | CPU | 2× Xeon Gold 6138 (40 cores), ~190 GB usable, Omni-Path |
| `p.gpu.ampere` | 11 | 4× NVIDIA **A100**-PCIE-40GB, 48 cores, 380 GB |
| `p.gpu` | 7 | 2× V100 (supported) or 2× P100 (pre-Volta, **not** supported), 40 cores |

## 0. One-time setup

```sh
ssh <user>@freya01.bc.mpcdf.mpg.de
cd ~/ptmp
git clone https://github.com/dutta-alankar/Chaa.git Chaa

# python for the validators / plotting (uv manages its own CPython)
curl -LsSf https://astral.sh/uv/install.sh | sh
uv venv ~/ptmp/Chaa/venv --python 3.12
~/ptmp/Chaa/venv/bin/python -m ensurepip 2>/dev/null || true
uv pip install --python ~/ptmp/Chaa/venv/bin/python numpy h5py matplotlib

# Chapel source (used by both the CPU and the GPU installations)
cd ~/ptmp/Chaa
wget https://github.com/chapel-lang/chapel/releases/download/2.8.0/chapel-2.8.0.tar.gz
tar xzf chapel-2.8.0.tar.gz                      # -> CPU Chapel tree
mkdir -p chapel-gpu && tar xzf chapel-2.8.0.tar.gz -C chapel-gpu   # -> GPU tree
```

!!! warning "Compiling on compute nodes: set `TMPDIR`"
    Compute-node `/tmp` is RAM-backed and counted against the job's
    memory cgroup.  Chapel's GPU code generation writes multi-GB
    temporaries there and the compile **hangs at ~0 % CPU** when it
    fills.  Both environment scripts export
    `TMPDIR=$HOME/ptmp/tmp` — keep that if you write your own jobs.

## 1. Chapel for single-locale programs (CPU)

Freya has no LLVM, so the CPU Chapel uses the C backend
(`CHPL_LLVM=none`) — that is all a single-locale (one node, all
cores) binary needs.  Build it once, in a batch job:

```sh
cd ~/ptmp/Chaa && sbatch tools/slurm/build-chapel.slurm
```

which does (on a `p.24h` node, ~20 min):

```sh
module purge && module load gcc/14 cmake/3.28
cd ~/ptmp/Chaa/chapel-2.8.0
export CHPL_HOME=$PWD  CHPL_LLVM=none
make -j32                          # compiler + CHPL_COMM=none runtime
```

Afterwards `source tools/slurm/freya-env.sh` provides the compiler
(plus `gcc/14`, `cmake/3.28`, `hdf5-serial/1.14.1` — note the hdf5
module is only visible **after** `gcc/14` is loaded).

**Build and run Chaa (single-locale).**

```sh
source tools/slurm/freya-env.sh
cmake -B build && cmake --build build          # -> build/bin/chaa
```

A single-locale job asks for one full node and simply runs the
binary; all 40 cores are used by Chapel's tasking layer:

```sh
#!/bin/bash -l
#SBATCH -p p.24h -N 1 -n 1 -c 40 --mem=100G -t 04:00:00
source ~/ptmp/Chaa/tools/slurm/freya-env.sh
cd ~/ptmp/Chaa
./build/bin/chaa --problem=sedov --geometry=cartesian \
    --nx1=256 --nx2=256 --nx3=256 --tstop=0.05 --outFormats=hdf5
```

(`tools/slurm/freya-tests.slurm` runs the full validated test suite
this way; `freya-bench-node.slurm` the within-node thread-scaling
benchmark.)

## 2. Chapel for multi-locale programs (CPU, GASNet)

Multi-locale needs a second **runtime** configuration in the same
Chapel tree: GASNet over InfiniBand/Omni-Path verbs.  The build job
above already does it; the manual step is:

```sh
cd ~/ptmp/Chaa/chapel-2.8.0 && export CHPL_HOME=$PWD CHPL_LLVM=none
module purge && module load gcc/14 cmake/3.28
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv \
  CHPL_LAUNCHER=slurm-gasnetrun_ibv make -j32 runtime
```

**Build Chaa (multi-locale).**  The comm setting is taken from the
environment at configure/compile time:

```sh
source tools/slurm/freya-env.sh
export CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv
cmake -B build-gasnet && cmake --build build-gasnet
```

**Submitting multi-locale jobs.**  Three Freya-specific rules:

1. **Never invoke the `chaa` launcher wrapper inside a batch job** —
   `slurm-gasnetrun_ibv` submits *its own* job (deadlock), and its
   ssh spawner is blocked from login to compute nodes.  Launch the
   real binary `chaa_real` directly with `srun --mpi=pmix` (GASNet
   was built with PMIx support).
2. **Cap the heap**: set `CHPL_RT_MAX_HEAP_SIZE` below the job's
   `--mem`, or GASNet's segment registration is OOM-killed at start.
3. `srun` steps do **not** inherit the batch `-c` — repeat it.

```sh
#!/bin/bash -l
#SBATCH -p p.24h -N 4 --ntasks-per-node=1 -c 40 --mem=180G -t 04:00:00
source ~/ptmp/Chaa/tools/slurm/freya-env.sh
export CHPL_RT_MAX_HEAP_SIZE=64g
cd ~/ptmp/Chaa
srun --mpi=pmix --ntasks-per-node=1 -c 40 -n $SLURM_NNODES \
     ./build-gasnet/bin/chaa_real -nl $SLURM_NNODES \
     --problem=sedov --geometry=cartesian \
     --nx1=512 --nx2=512 --nx3=512 --tstop=0.05
```

(one Chapel locale per node, 40 threads each;
`freya-bench-multi.slurm` runs the across-node strong/weak-scaling
campaign this way).

## 3. Chapel for GPU programs

GPU code generation requires Chapel's **LLVM backend**, which the
C-backend installation above does not have — so the GPU Chapel is a
**second, separate installation** built with the bundled LLVM against
the CUDA toolkit (no GPU is needed at build time):

```sh
cd ~/ptmp/Chaa && sbatch tools/slurm/build-chapel-gpu.slurm
```

which does (on a `p.24h` node, ~20 min for the compiler + both GPU
runtimes; LLVM itself is built once):

```sh
module purge && module load gcc/14 cuda/12.8 cmake/3.28
cd ~/ptmp/Chaa/chapel-gpu/chapel-2.8.0
export CHPL_HOME=$PWD
export CHPL_LLVM=bundled          # LLVM backend: required for GPUs
export CHPL_LOCALE_MODEL=gpu
export CHPL_GPU=nvidia
export CHPL_CUDA_PATH=$CUDA_HOME
export CHPL_GPU_ARCH=sm_80        # A100; sm_70 for the V100 nodes
make -j32                         # compiler + comm=none GPU runtime

# multi-node GPU runtime (GASNet/ibv), same pattern as the CPU one:
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv \
  CHPL_LAUNCHER=slurm-gasnetrun_ibv make -j32 runtime
```

Afterwards `source tools/slurm/freya-gpu-env.sh` selects this
installation (it also exports `LD_LIBRARY_PATH=$HDF5_HOME/lib:…` —
on the GPU nodes the hdf5 module sets `HDF5_HOME` but not the runtime
library path — and the `TMPDIR` fix above).

**Build Chaa (GPU).**  One CMake flag; compiling works on the login
node:

```sh
source tools/slurm/freya-gpu-env.sh
cmake -B build-gpu -DCHAA_GPU=ON && cmake --build build-gpu

# V100 variant (p.gpu partition) — needs the sm_70 Chapel tree: the
# GPU runtime bakes its device helpers for one architecture, so an
# sm_80 runtime cannot serve sm_70 binaries (build the tree once with
# build-chapel-gpu-v100.slurm, i.e. the same recipe with
# CHPL_GPU_ARCH=sm_70 in ~/ptmp/Chaa/chapel-gpu-v100):
export CHPL_HOME=$HOME/ptmp/Chaa/chapel-gpu-v100/chapel-2.8.0
export PATH="$CHPL_HOME/bin/linux64-x86_64:$PATH"
CHPL_GPU_ARCH=sm_70 cmake -B build-gpu-v100 -DCHAA_GPU=ON
cmake --build build-gpu-v100

# multi-node GPU binary (gasnet runtime):
export CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv
cmake -B build-gpu-gasnet -DCHAA_GPU=ON && cmake --build build-gpu-gasnet
```

**Submitting single-node GPU jobs** (`p.gpu.ampere`, 4× A100): one
process drives all GPUs visible to it — the grid splits across them
along x1 — so a single-node job needs one task and a `--gres`
request; `CUDA_VISIBLE_DEVICES` selects a subset:

```sh
#!/bin/bash -l
#SBATCH --partition=p.gpu.ampere
#SBATCH --constraint="gpu"
#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:4    # 4 GPUs of this node
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH -t 04:00:00
source ~/ptmp/Chaa/tools/slurm/freya-gpu-env.sh
cd ~/ptmp/Chaa
./build-gpu/bin/chaa --problem=sedov ... # all 4 GPUs
CUDA_VISIBLE_DEVICES=0 ./build-gpu/bin/chaa ...   # exactly one GPU
```

For the Volta nodes use `--partition=p.gpu --gres=gpu:v100:2` and the
`build-gpu-v100` binary (the P100 nodes in the same partition are
pre-Volta and rejected by `ptxas` — see
[Running on GPUs](user-guide/gpu.md#current-limitations)).

**Submitting multi-node GPU jobs**: one task (co-locale) per GPU with
the gasnet binary — the CPU multi-locale rules (PMIx launch, heap
cap, repeated `-c`) apply unchanged, plus
`CHPL_RT_LOCALES_PER_NODE=<gpus per node>` so each co-locale binds
one GPU:

```sh
#!/bin/bash -l
#SBATCH --partition=p.gpu.ampere
#SBATCH --constraint="gpu"
#SBATCH --nodes=2            # 1 or more full GPU nodes
#SBATCH --gres=gpu:a100:4    # 4 GPUs per node
#SBATCH --ntasks-per-node=4  # one co-locale per GPU
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH -t 01:00:00
source ~/ptmp/Chaa/tools/slurm/freya-gpu-env.sh
export CHPL_RT_LOCALES_PER_NODE=4
export CHPL_RT_MAX_HEAP_SIZE=32g       # per co-locale, under --mem/4
cd ~/ptmp/Chaa
NL=$((4*SLURM_NNODES))
srun --mpi=pmix --ntasks-per-node=4 -c 12 -n $NL \
     ./build-gpu-gasnet/bin/chaa_real -nl $NL --problem=sedov ...
```

Ready-made jobs: `freya-gpu-tests.slurm` (GPU-vs-CPU comparison
matrix + the full validated suite on an A100), `freya-gpu-bench.slurm`
(A100 size scan + 1/2/4-GPU scaling), `freya-gpu-bench-pgpu.slurm`
(V100), `freya-gpu-multi.slurm` (2 nodes × 4 GPUs).

## Quick reference

| goal | Chapel build | Chaa build | job pattern |
|---|---|---|---|
| 1 node, CPU | `CHPL_LLVM=none make` | `cmake -B build` | `-N 1 -c 40`, run `build/bin/chaa` |
| N nodes, CPU | + `CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv make runtime` | `CHPL_COMM=gasnet … cmake -B build-gasnet` | `srun --mpi=pmix -n N chaa_real -nl N` |
| 1 node, GPUs | separate tree: `CHPL_LLVM=bundled CHPL_LOCALE_MODEL=gpu CHPL_GPU=nvidia make` | `cmake -B build-gpu -DCHAA_GPU=ON` | `--gres=gpu:a100:4 -n 1`, run `build-gpu/bin/chaa` |
| N nodes, GPUs | + gasnet GPU runtime | `CHPL_COMM=gasnet … -DCHAA_GPU=ON` | `--ntasks-per-node=4` + `CHPL_RT_LOCALES_PER_NODE=4`, `srun --mpi=pmix chaa_real -nl 4N` |
