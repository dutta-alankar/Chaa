# Benchmarks & scaling

Chaa's figure of merit is the **Mcell-updates/s** printed at the end of
every run (steps × cells / wall time of the evolution loop; I/O and
setup excluded). The benchmark problem throughout is the 3D Cartesian
Sedov blast (hllc/PLM/RK2, 50 steps, no output), run by
`tools/bench.sh` (laptop) and the `tools/slurm/freya-bench-*.slurm`
scripts (cluster); the figures are produced by `tools/plot_bench.py`.

Two machines are covered: a production HPC cluster (**MPCDF Freya**)
and an Apple-silicon laptop.

## MPCDF Freya (production cluster)

**Hardware & software.** Each compute node has 2× Intel Xeon Gold 6138
(Skylake-SP, 20 cores each, 2.0 GHz base) — 40 cores and ~190 GB usable
memory per node — connected by a 100 Gb/s Intel **Omni-Path** fabric
(`hfi1`); SLES 15, SLURM. Chaa was compiled with Chapel 2.8 (C
backend, `CHPL_LLVM=none`) and gcc 14.1; the multi-locale binary uses
the GASNet-EX **ibv** conduit over Omni-Path verbs, bootstrapped by
`srun --mpi=pmix`, with one Chapel locale per node and 40 threads per
locale. HDF5 1.14.1 (serial). The complete setup and launch recipe is
in [`tools/slurm/README.md`](https://github.com/dutta-alankar/Chaa/blob/main/tools/slurm/README.md).

**Correctness first:** the full 47-case validated test suite passes on
a Freya compute node in 3 m 44 s (144 quantitative checks).

### Within one node (thread scaling, `CHPL_COMM=none`)

| threads | grid | Mcell/s | speed-up |
|---|---|---|---|
| 1 | 128³ | 1.04 | — |
| 2 | 128³ | 2.01 | 1.93× |
| 5 | 128³ | 4.47 | 4.30× |
| 10 | 128³ | 7.71 | 7.41× |
| 20 | 128³ | 10.30 | 9.90× |
| 40 | 128³ | 14.77 | 14.2× |

Weak scaling (64³ cells per 5 threads) tells the same story: 4.33 →
7.57 → 9.97 → 14.76 Mcell/s at 5/10/20/40 threads. Scaling is
near-ideal while a socket's memory bandwidth lasts and saturates
beyond ~10 cores per socket — the expected behaviour of a
bandwidth-bound stencil code on Skylake; the second socket doubles the
available bandwidth and the rate follows.

![Freya within-node scaling](assets/plots/freya-node-scaling.png)

### Across nodes (gasnet/ibv, one 40-thread locale per node)

Strong scaling, fixed 256³ box:

| locales (nodes) | Mcell/s | speed-up |
|---|---|---|
| 1 | 8.50 | — |
| 2 | 12.21 | 1.44× |
| 4 | 18.35 | 2.16× |
| 8 | 18.52 | 2.18× |

Weak scaling, 256³ cells per locale:

| locales | grid | Mcell/s | efficiency |
|---|---|---|---|
| 1 | 256³ | 8.48 | — |
| 2 | 512·256·256 | 14.90 | 88 % |
| 4 | 512·512·256 | 24.90 | 73 % |
| 8 | 512³ | 44.13 | 65 % |

![Freya across-node scaling](assets/plots/freya-multi-scaling.png)

How to read these numbers:

- **Weak scaling is the operative regime** for a distributed
  finite-volume code: grow the problem with the machine. 512³
  (134 million cells) runs at 44 Mcell/s on 8 nodes at 65 % efficiency.
- **Strong scaling saturates** once the per-locale block gets small:
  256³ split 8 ways leaves only 2 M cells per 40-core node, and the
  per-step distributed-loop synchronisations (a fixed ~tens of ms per
  step) stop shrinking. Use ≳8 M cells per node for efficient runs.
- The multi-locale binary pays a **flat ~40 % single-locale penalty**
  versus the `CHPL_COMM=none` build (8.5 vs 14.8 Mcell/s on the same
  node) — the cost of compiling every array access for a potentially
  remote address space with the C backend. Use the single-locale build
  whenever a run fits on one node.

Two multi-locale performance lessons from this campaign are now baked
into the code (they were invisible on shared memory and dominant on a
real network): restart/gather I/O uses bulk per-plane transfers
instead of element-wise remote reads, and the grid is block-split
**only along x1**, keeping every halo plane contiguous in memory (the
default 2×2×2 locale grid on 8 locales cuts the memory-fastest axis
into tens of thousands of tiny strided RDMAs and cost a third of the
8-node throughput).

## Apple-silicon laptop (8 cores: 4 performance + 4 efficiency)

Measured with `tools/bench.sh`, Chapel 2.8 (LLVM backend), macOS.

Single locale (`CHPL_COMM=none`), strong scaling at 128³:

| threads | Mcell/s | speed-up |
|---|---|---|
| 1 | 2.79 | — |
| 2 | 5.37 | 1.92× |
| 4 | 9.18 | 3.29× |
| 8 | 8.93 | 3.20× |

Weak scaling (64³ per thread): 2.73 / 5.21 / 9.00 / 9.12 Mcell/s at
1/2/4/8 threads. The 4 efficiency cores add nothing to this
bandwidth-bound kernel — 4 threads is the sweet spot.

Multi-locale (GASNet **smp** conduit, all locales on the one machine,
fixed 8 threads total — this isolates the distributed code path's
overhead, not network scaling):

| locales × threads | Mcell/s (strong, 128³) |
|---|---|
| 1 × 8 | 5.93 |
| 2 × 4 | 5.62 |
| 4 × 2 | 5.44 |

i.e. ~8 % total cost for splitting the box four ways over shared
memory. Multi-locale correctness is verified here too: 4-locale piece
output reassembles to the single-locale fields to machine precision
(2.6×10⁻¹⁵) and particle trajectories match to 10⁻¹².

## Reproducing

```sh
# laptop / any shared-memory machine:
tools/bench.sh                            # single-locale thread scaling
tools/bench.sh build-gasnet/bin/chaa 4    # + gasnet-smp locale overhead

# cluster (Freya; adapt modules/substrate elsewhere):
sbatch      tools/slurm/freya-bench-node.slurm
sbatch -N 8 tools/slurm/freya-bench-multi.slurm
python tools/plot_bench.py chaa-bench-node-*.out  --save node.png
python tools/plot_bench.py chaa-bench-multi-*.out --save multi.png
```

See [Running in parallel](user-guide/parallel.md) for building the
multi-locale runtime and [`tools/slurm/README.md`](https://github.com/dutta-alankar/Chaa/blob/main/tools/slurm/README.md)
for the full cluster recipe.
