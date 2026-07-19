# Chaa on MPCDF Freya — setup & running guide

> The same material, rendered with the rest of the documentation:
> [Running on Freya (MPCDF)](https://dutta-alankar.github.io/Chaa/freya/)
> — including the GPU Chapel build and GPU job patterns.

Step-by-step instructions for building Chapel (single- and
multi-locale) on [MPCDF Freya](https://docs.mpcdf.mpg.de) and compiling
and running Chaa there. Everything below was validated on Freya
(2× Intel Xeon Gold 6138, 40 cores/node, Omni-Path fabric, SLURM);
for another cluster, adapt the modules, the GASNet substrate and the
partition names.

All paths assume the layout used by the scripts in this directory:

```
~/ptmp/Chaa/                  work directory (on the ptmp filesystem)
├── chapel-2.8.0/             Chapel source tree (built in place)
├── venv/                     uv-managed python (validation, plotting)
├── build/                    Chaa, single-locale (CHPL_COMM=none)
├── build-gasnet/             Chaa, multi-locale (CHPL_COMM=gasnet)
└── ...                       the Chaa git checkout itself
```

## 1. One-time setup

```sh
mkdir -p ~/ptmp/Chaa && cd ~/ptmp/Chaa
git clone https://github.com/dutta-alankar/Chaa.git .        # or your fork

# Chapel source
curl -LO https://github.com/chapel-lang/chapel/releases/download/2.8.0/chapel-2.8.0.tar.gz
tar xzf chapel-2.8.0.tar.gz

# python for the test-suite validators and the plotting tools
# (uv downloads a managed CPython; no module needed)
curl -LsSf https://astral.sh/uv/install.sh | sh
uv venv venv --python 3.12
uv pip install --python venv/bin/python numpy h5py matplotlib
```

The module environment used everywhere (this is what
[`freya-env.sh`](freya-env.sh) loads):

```sh
module purge
module load gcc/14 cmake/3.28 hdf5-serial/1.14.1
export CHPL_HOME=~/ptmp/Chaa/chapel-2.8.0
export CHPL_LLVM=none          # C backend: no LLVM dependency on the cluster
export PATH=$CHPL_HOME/bin/linux64-x86_64:$PATH
```

Notes:

- `hdf5-serial` is only visible after loading `gcc/14` (module
  hierarchy). Serial HDF5 is the right choice: Chaa's multi-locale
  output writes independent per-locale piece files, so `hdf5-parallel`
  would only add an unused MPI dependency.
- `CHPL_LLVM=none` selects Chapel's C backend — it avoids building or
  finding LLVM on the cluster and compiles Chaa with the loaded gcc.

## 2. Building Chapel

Build on a compute node (it matches the target CPU and keeps the login
node free) — [`build-chapel.slurm`](build-chapel.slurm) does exactly
this:

```sh
sbatch tools/slurm/build-chapel.slurm     # ~15 min on a 40-core node
```

or by hand:

```sh
cd $CHPL_HOME
# single-locale compiler + runtime (CHPL_COMM=none):
make -j32
# multi-locale runtime: GASNet over the Omni-Path verbs interface,
# SLURM-aware launcher
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv \
  CHPL_LAUNCHER=slurm-gasnetrun_ibv make -j32 runtime
```

Both runtimes coexist in one tree; the `CHPL_COMM` environment at
*application compile time* selects which one a binary uses.

## 3. Building Chaa

```sh
cd ~/ptmp/Chaa
source tools/slurm/freya-env.sh

# single-locale binary (node-local runs, the test suite):
cmake -B build && cmake --build build

# multi-locale binary:
export CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv CHPL_LAUNCHER=slurm-gasnetrun_ibv
cmake -B build-gasnet && cmake --build build-gasnet
```

The multi-locale build produces `build-gasnet/bin/chaa` (a launcher)
and `build-gasnet/bin/chaa_real` (the actual program — this is what
the SLURM scripts run).

## 4. Running single-locale

One node, all cores; ordinary sbatch. To run the full validated test
suite, `sbatch tools/slurm/freya-tests.slurm`; for a simulation:

```sh
#!/bin/bash -l
#SBATCH -p p.24h -N 1 -n 1 -c 40 --mem=100G -t 04:00:00
cd ~/ptmp/Chaa
source tools/slurm/freya-env.sh
export CHPL_RT_NUM_THREADS_PER_LOCALE=40
./build/bin/chaa --problem=sedov --geometry=spherical --nx1=1024 ...
```

## 5. Running multi-locale

Use [`freya-chaa.slurm`](freya-chaa.slurm) — one Chapel locale per
node, 40 threads per locale, arguments passed through to chaa:

```sh
sbatch -N 4 tools/slurm/freya-chaa.slurm \
       --problem=sedov --geometry=cartesian \
       --nx1=512 --nx2=512 --nx3=512 --tstop=0.5 --outFormats=hdf5
```

The essential pattern inside (worth knowing when writing your own
scripts):

```sh
export CHPL_RT_NUM_THREADS_PER_LOCALE=40
export CHPL_RT_MAX_HEAP_SIZE=64g
srun --mpi=pmix --ntasks-per-node=1 -c 40 -n $SLURM_NNODES \
     ./build-gasnet/bin/chaa_real -nl $SLURM_NNODES [chaa args...]
```

Three things make this work — all learned the hard way:

1. **Launch `chaa_real` directly with `srun --mpi=pmix`.** GASNet was
   compiled with PMI(x) spawner support, so srun bootstraps the
   locales inside the allocation. Do *not* run the `chaa` launcher
   binary inside an sbatch script: Chapel's `slurm-gasnetrun_ibv`
   launcher submits its own allocation (deadlocking yours), and its
   default ssh spawner needs node-to-node ssh that Freya doesn't
   allow.
2. **Cap the heap with `CHPL_RT_MAX_HEAP_SIZE`.** The ibv conduit
   registers ("pins") the Chapel heap at startup; the default tries to
   take most of the node's physical memory, which exceeds the job's
   cgroup limit and gets the job OOM-killed at startup. Pick a value
   comfortably under `--mem` and large enough for your grid (fields
   need ≈ 240 B/cell/locale plus halo).
3. **Repeat `-c 40` on the srun line.** Recent SLURM versions do not
   propagate `--cpus-per-task` from sbatch to srun steps; without it
   each locale may be confined to one core.

Multi-locale output arrives as per-locale HDF5/VTK piece files plus an
XDMF collection; the python tools (`tools/chaa_io.py`,
`tools/plot_fields.py`) reassemble pieces transparently.

## 6. Benchmarks

The scaling campaigns from the
[benchmarks documentation](https://dutta-alankar.github.io/Chaa/benchmarks/)
are reproduced with:

```sh
sbatch      tools/slurm/freya-bench-node.slurm    # within-node thread scaling
sbatch -N 8 tools/slurm/freya-bench-multi.slurm   # across-node strong + weak
python tools/plot_bench.py chaa-bench-*.out --save scaling.png
```

## 7. GPUs (p.gpu.ampere: 4x A100 per node)

GPU code generation needs Chapel's LLVM backend, so a **second**
Chapel installation is built (the CPU one above uses the C backend).
Everything is scripted; details and design notes live in the
[GPU documentation](https://dutta-alankar.github.io/Chaa/user-guide/gpu/).

```sh
# one-time: GPU Chapel (bundled LLVM + CUDA + gasnet runtime, ~2 h)
mkdir -p ~/ptmp/Chaa/chapel-gpu && cd ~/ptmp/Chaa/chapel-gpu
wget https://github.com/chapel-lang/chapel/releases/download/2.8.0/chapel-2.8.0.tar.gz
tar xzf chapel-2.8.0.tar.gz
cd ~/ptmp/Chaa && sbatch tools/slurm/build-chapel-gpu.slurm

# Chaa GPU binary (compiles on the login node; no GPU needed)
source tools/slurm/freya-gpu-env.sh
cmake -B build-gpu -DCHAA_GPU=ON && cmake --build build-gpu

# validation on an A100 node: GPU-vs-CPU comparison matrix + full suite
sbatch tools/slurm/freya-gpu-tests.slurm

# benchmarks: size scan + 1/2/4-GPU strong & weak scaling
sbatch tools/slurm/freya-gpu-bench.slurm
```

Single-node runs drive all visible GPUs from one process
(`CUDA_VISIBLE_DEVICES` selects a subset).  Multi-node runs launch one
co-locale per GPU with the gasnet binary:

```sh
srun --mpi=pmix --ntasks-per-node=4 -n $((4*SLURM_NNODES)) \
     build-gpu-gasnet/bin/chaa_real -nl $((4*SLURM_NNODES))
```

The `CHPL_RT_MAX_HEAP_SIZE`/PMIx caveats from section 5 apply
unchanged; `CHPL_GPU_ARCH=sm_80` targets the A100s and `sm_70` the
V100s in `p.gpu` (`sbatch tools/slurm/freya-gpu-bench-pgpu.slurm`).
The P100 nodes are pre-Volta and **not supported**: Chaa's kernels
exceed their 4 KB kernel-parameter limit (`ptxas` rejects sm_60).
