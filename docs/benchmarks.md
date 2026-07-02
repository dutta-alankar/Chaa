# Benchmarks & scaling

Chaa's figure of merit is the **Mcell-updates/s** printed at the end of
every run (steps × cells / wall time). The numbers below were measured
with `tools/bench.sh` on the 3D Cartesian Sedov blast (hllc/linear/rk2,
no I/O, 50 steps), on an Apple-silicon laptop (8 cores: 4 performance +
4 efficiency), Chapel 2.8, `--fast`.

Reproduce them with

```sh
tools/bench.sh                            # single-locale (threads)
tools/bench.sh build-gasnet/bin/chaa 4    # + multi-locale (gasnet)
```

## Single locale — strong scaling

Fixed 128³ grid, varying `CHPL_RT_NUM_THREADS_PER_LOCALE`
(`CHPL_COMM=none` build):

| threads | Mcell/s | speed-up | efficiency |
|---|---|---|---|
| 1 | 2.50 | 1.0× | — |
| 2 | 4.95 | 1.98× | 99 % |
| 4 | 8.05 | 3.23× | 81 % |
| 8 | 7.91 | 3.17× | — |

Scaling is near-ideal across the 4 performance cores; the 4 efficiency
cores add nothing to this memory-bandwidth-bound stencil kernel (a
common result on hybrid Apple-silicon parts), so 4 threads is the sweet
spot on this machine.

## Single locale — weak scaling

64³ cells per thread (the grid grows with the thread count):

| threads | grid | Mcell/s | efficiency |
|---|---|---|---|
| 1 | 64³ | 2.60 | — |
| 2 | 128·64·64 | 5.00 | 96 % |
| 4 | 128·128·64 | 7.85 | 76 % |
| 8 | 128³ | 7.34 | — |

## Multi-locale — strong scaling

`CHPL_COMM=gasnet` (smp conduit) build, all locales on the one node,
**fixed total thread count** (8) so that only the communication /
partitioning overhead varies:

| locales | threads/locale | grid | Mcell/s | efficiency |
|---|---|---|---|---|
| 1 | 8 | 128³ | 5.56 | — |
| 2 | 4 | 128³ | 5.35 | 96 % |
| 4 | 2 | 128³ | 4.98 | 90 % |

Splitting the same box over 4 locales costs only ~10 %: the halo
exchanges (`updateFluff()`) and the distributed-array bookkeeping are
cheap relative to the hydro sweeps. (The gasnet binary is ~30 % slower
than the `CHPL_COMM=none` binary at equal resources — the standard cost
of compiling array accesses for a potentially remote address space.)

## Multi-locale — weak scaling

128·64·64 cells **per locale**, 2 threads per locale:

| locales | grid | Mcell/s | efficiency |
|---|---|---|---|
| 1 | 128·64·64 | 3.08 | — |
| 2 | 128·128·64 | 5.11 | 83 % |
| 4 | 128³ | 4.00 | 32 % |

The 4-locale point is limited by the machine, not the code: 4 locales ×
2 threads oversubscribes onto the efficiency cores (compare the
single-locale strong-scaling table, where 8 native threads are no
faster than 4).

## Reading these numbers

Everything here ran on **one laptop** — the multi-locale runs use
GASNet's shared-memory conduit, so they demonstrate correctness and
measure the *overhead* of the distributed-memory code path, not real
network scaling. On a cluster the strong-scaling efficiency per
doubling (96 %, 90 % above) is the quantity to watch, with the network
replacing shared memory in the halo exchange. Piece-wise parallel
output (per-locale HDF5/VTK + XDMF stitching) was verified in these
runs to reproduce the single-locale fields to machine precision
(max |Δρ| = 2.6×10⁻¹⁵), and distributed tracer particles match the
single-locale trajectories to 10⁻¹².

See [Running in parallel](user-guide/parallel.md) for how to build the
gasnet runtime and launch multi-locale jobs.
