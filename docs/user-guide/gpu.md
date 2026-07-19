# Running on GPUs

Chaa can hand the hot loops — reconstruction, Riemann solves, flux
divergence, source terms, the RK update, floors, cooling and the CFL
reduction — to NVIDIA GPUs.  GPU execution is a **compile-time**
choice: the same source builds either the classic CPU binary (default,
completely unchanged) or a GPU binary, selected by one CMake flag.

```sh
cmake -B build-gpu -DCHAA_GPU=ON      # needs a GPU-enabled chpl (below)
cmake --build build-gpu
./build-gpu/bin/chaa --problem=sedov --nx1=256 --nx2=256 --nx3=256 ...
```

(The GPU build compiles with `--gpu-block-size=128` —
`-DCHAA_GPU_BLOCK_SIZE=…` to change it: the flux kernel is
register-heavy and Chapel's default of 512 threads per block exceeds
the per-block register file, failing at launch with CUDA error 701.)

## How the GPU port works

The numerical kernels are written once, as generic procedures shared
by both builds; only where the data lives differs.

- Each locale carves its slab of the (x1-split) grid into one **padded
  block per visible GPU**.  A block owns device-resident copies of the
  state (`V`, `U`, `U0`, `RHS`, `FLX`, the face states, the solve
  mask, forcing tables and, if enabled, the gravitational potential),
  so a whole RK stage runs as a sequence of device kernels with **no
  host round-trips for the interior**.
- Every kernel iterates a **flattened 1D index space** (Chapel 2.8
  parallelises a GPU forall only over one dimension — a
  multidimensional forall runs ~50–100× below memory bandwidth), and
  the flux evaluation is split into **reconstruction and Riemann
  kernels specialised at compile time to the active scheme** (a fused
  all-scheme kernel needs 255 registers plus ~16 KB of spilled locals
  per thread).
- Face fluxes are recomputed on both sides of a block cut (one extra
  face plane per block), so blocks never exchange fluxes.
- **Standard x2/x3 boundary conditions run on the device** (those
  directions are never split across blocks): periodic, zero-gradient,
  reflect, axis, inflow and the diodes are filled by a ghost-slab
  kernel.  Only the thin **x1 planes** are staged through the host
  `StencilDist` arrays each stage — x1 physical BCs and inter-block /
  inter-locale halo exchange run on the host, as do `userdef` and
  `shear-periodic` sides of any direction (automatic fallback), so
  every BC type behaves exactly as on the CPU.
- Output, restart files, Lagrangian particles, FARGO's row remap and
  the self-gravity CG solve keep their host implementations; the
  interior is pulled off the devices on demand (at dump cadence for
  output/restart, once per step for particles/FARGO/self-gravity).
  Restart dumps remain **byte-identical across a stop/resume** of the
  same GPU binary and configuration.

At startup the GPU build runs a kernel-launch self-test and reports

```
  gpus       : 4 visible across 1 locale(s)
  gpu        : 4 device(s), self-test ok (... kernel launches)
```

halting if the loops silently fell back to host execution (which would
be orders of magnitude slower, not wrong).  `--gpuDiag=true` prints
kernel-launch and host↔device transfer counts per locale at the end of
the run — useful to confirm that per-step traffic is only the ghost
shell.

### Current limitations

- **Volta (sm_70) or newer is required.**  Chapel passes kernel
  arguments by value, and Chaa's physics kernels exceed the 4 KB
  kernel-parameter limit of pre-Volta GPUs (`ptxas` rejects an
  `sm_60` build); Volta and later allow 32 KB.  V100 and A100 are
  validated; P100 (Pascal) is not supported.
- **At most two concurrently driven GPUs per process.**  Chapel
  2.8's runtime faults (CUDA "invalid resource handle") when one
  locale drives more than two devices; two scale well (weak
  efficiency 94 %).  For more GPUs, scale across nodes/locales with
  the GASNet binary — one locale per node or per GPU pair.
- **One GPU architecture per Chapel runtime.**  The GPU runtime bakes
  its built-in device helpers (reductions etc.) for `CHPL_GPU_ARCH`
  at runtime-build time — a binary compiled for `sm_70` against an
  `sm_80` runtime crashes at launch.  Mixed clusters need one Chapel
  tree per architecture (on Freya: `chapel-gpu` for the A100s,
  `chapel-gpu-v100` for the V100s).
- `cylinderFlow` (per-stage internal/immersed boundaries) and problems
  using the `problemBodyForce` hook are rejected at startup — these
  hooks run per-stage on the host interior; use the CPU build.
- Only `CHPL_GPU=nvidia` has been tested (AMD may work — Chapel
  supports it — but is unvalidated).
- Chapel 2.8's `stencilDist` cannot itself target GPU sublocales;
  that is why Chaa manages the device blocks explicitly.

## Building a GPU-enabled Chapel

GPU code generation needs Chapel's **LLVM backend** — a Chapel built
with `CHPL_LLVM=none` (e.g. the C-backend cluster build used for the
CPU binary) cannot compile GPU code, so a second Chapel installation
is required:

```sh
wget https://github.com/chapel-lang/chapel/releases/download/2.8.0/chapel-2.8.0.tar.gz
tar xzf chapel-2.8.0.tar.gz && cd chapel-2.8.0
export CHPL_HOME=$PWD
export CHPL_LLVM=bundled          # builds LLVM too: ~1.5-2 h, once
export CHPL_LOCALE_MODEL=gpu
export CHPL_GPU=nvidia
export CHPL_CUDA_PATH=$CUDA_HOME  # CUDA toolkit prefix (>= 11.7)
export CHPL_GPU_ARCH=sm_80        # A100; sm_70 V100, sm_60 P100,
                                  # sm_90 H100
make -j
export PATH=$CHPL_HOME/bin/linux64-x86_64:$PATH
```

No GPU is needed at build time, only the CUDA toolkit.  For
**multi-node** GPU runs additionally build the GASNet runtime
configuration (same as for CPU multi-locale):

```sh
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=ibv \
  CHPL_LAUNCHER=slurm-gasnetrun_ibv make -j runtime
```

There is also a **cpu-as-device** mode (`CHPL_GPU=cpu` instead of
`nvidia`, no CUDA needed) that compiles the GPU code path but executes
it on the CPU — slow, but ideal for functional testing on machines
without a GPU; the GPU comparison test below runs in this mode in CI.

## Running

### Single node

One process drives all visible GPUs — the grid is split across them
along x1:

```sh
./build-gpu/bin/chaa --problem=sedov --nx1=256 --nx2=256 --nx3=256 ...
CUDA_VISIBLE_DEVICES=0     ./build-gpu/bin/chaa ...   # exactly 1 GPU
CUDA_VISIBLE_DEVICES=0,1   ./build-gpu/bin/chaa ...   # 2 GPUs
```

### Multiple nodes

Build Chaa against the gasnet+GPU runtime (`CHPL_COMM=gasnet
CHPL_COMM_SUBSTRATE=ibv` in the environment at `cmake`/build time,
build directory e.g. `build-gpu-gasnet`) and launch **one co-locale
per GPU**, exactly like the CPU multi-locale runs but with
`--ntasks-per-node=<gpus per node>`:

```sh
srun --mpi=pmix --ntasks-per-node=4 -n $((4*SLURM_NNODES)) \
     ./build-gpu-gasnet/bin/chaa_real -nl $((4*SLURM_NNODES))
```

Chapel binds each co-locale to one GPU automatically; Chaa then runs
one device block per locale.  All the CPU-cluster caveats from
[Running in parallel](parallel.md) (heap size, PMIx launch) apply
unchanged, plus `GASNET_PHYSMEM_MAX` (cap the GPU runtime's
pinnable-memory probe inside a job cgroup).

!!! warning "Multi-node GPU: correct, not yet fast"
    Across nodes every stage's x1 ghost exchange travels
    device → host → network → host → device, which dominates the step
    at GPU speeds (measured ~4 Mcell/s on 2 nodes vs ~100+ on one
    A100).  Scale GPUs *within* a node where possible; treat
    multi-node GPU as a correctness/capacity option until direct
    device-to-device halo exchange lands.

## Freya (MPCDF) walkthrough

!!! tip
    [Running on Freya](../freya.md) is the complete cluster guide —
    Chapel builds for single-locale, multi-locale **and** GPU
    programs, and the SLURM submission pattern for each; below is the
    GPU-only summary.

Freya has 11 nodes with **4× NVIDIA A100** each (partition
`p.gpu.ampere`, 48 cores, 380 GB) plus P100/V100 nodes (`p.gpu`).
Everything below is scripted in `tools/slurm/`:

| script | purpose |
|---|---|
| `freya-gpu-env.sh` | modules (`gcc/14 cuda/12.8 cmake/3.28 hdf5-serial/1.14.1`) + `CHPL_*` GPU environment |
| `build-chapel-gpu.slurm` | one-off Chapel build (bundled LLVM + CUDA + gasnet runtime) |
| `freya-gpu-tests.slurm` | GPU-vs-CPU comparison matrix + the full validated suite on an A100 |
| `freya-gpu-bench.slurm` | single-GPU size scan, 1/2/4-GPU strong and weak scaling |

```sh
# one-off: Chapel with GPU support (~2 h batch job)
mkdir -p ~/ptmp/Chaa/chapel-gpu && cd ~/ptmp/Chaa/chapel-gpu
wget https://github.com/chapel-lang/chapel/releases/download/2.8.0/chapel-2.8.0.tar.gz
tar xzf chapel-2.8.0.tar.gz
cd ~/ptmp/Chaa && sbatch tools/slurm/build-chapel-gpu.slurm

# Chaa GPU binary (login node is fine; no GPU needed to compile)
source tools/slurm/freya-gpu-env.sh
cmake -B build-gpu -DCHAA_GPU=ON && cmake --build build-gpu

# interactive A100 example
srun -p p.gpu.ampere --constraint=gpu --gres=gpu:a100:4 -N 1 -n 1 \
     -c 48 --mem=180G -t 0:30:00 --pty bash -l
source tools/slurm/freya-gpu-env.sh
./build-gpu/bin/chaa --problem=sedov --nx1=256 --nx2=256 --nx3=256 \
    --x1min=-0.5 --x1max=0.5 --x2min=-0.5 --x2max=0.5 \
    --x3min=-0.5 --x3max=0.5 --tstop=0.05 --outFormats=hdf5
```

A batch job requesting whole GPU nodes looks like:

```sh
#SBATCH --partition=p.gpu.ampere
#SBATCH --nodes=2            # 1 or more full GPU nodes
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:4    # 4 GPUs per node
#SBATCH --ntasks-per-node=4  # one co-locale per GPU (multi-node)
#SBATCH --cpus-per-task=12
```

(for Volta use `--partition=p.gpu` with `--gres=gpu:v100:2` and a
binary compiled with `CHPL_GPU_ARCH=sm_70`; the P100 nodes in `p.gpu`
are pre-Volta and not supported — see the limitations above).

## Validating a GPU build

`tests/run_gpu_compare.sh` runs 15 small configurations covering every
GPU-relevant code path (all reconstructions, integrators and Riemann
solvers, curvilinear geometry, viscosity, conduction, cooling, the
isothermal EOS, user BCs, OU forcing, particles, FARGO, self-gravity
and a stop/resume restart) with the GPU binary and requires agreement
with the CPU binary to round-off (and bit-identical restart
continuation):

```sh
CHAA_BIN=build/bin/chaa CHAA_GPU_BIN=build-gpu/bin/chaa \
  tests/run_gpu_compare.sh
```

On a GPU node the full validated test suite also runs against the GPU
binary (`tools/slurm/freya-gpu-tests.slurm`); on a Freya A100 the
whole comparison matrix passes and the suite scores **46 passed, 0
failed** (`cylinder-flow` is skipped by design — internal BCs are
CPU-only).

Measured GPU performance is collected in
[Benchmarks & scaling](../benchmarks.md).
