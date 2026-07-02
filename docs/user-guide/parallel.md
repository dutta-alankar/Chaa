# Running in parallel

Chaa's fields live on a single block-distributed domain
(`StencilDist`) with halo ("fluff") regions. Every loop in the solver
is a data-parallel `forall`; halo refreshes happen in one call,
`updateFluff()`. There is no other communication logic anywhere in the
code — Chapel moves the data. Tracer particles are distributed the
same way (each locale owns the particles inside its grid block; see
[Tracer particles](particles.md)).

## Shared memory (default)

A standard build (`CHPL_COMM=none`) uses all cores of the node
automatically:

```sh
CHPL_RT_NUM_THREADS_PER_LOCALE=8 ./build/bin/chaa --problem=sedov ...
```

Measured thread scaling is near-ideal per physical performance core —
see [Benchmarks & scaling](../benchmarks.md) for strong- and
weak-scaling tables.

## Multi-locale (distributed)

Build Chapel with a communication layer, rebuild Chaa, and launch with
`-nl`. For a single machine (testing, development) the GASNet **smp**
conduit is enough — with a Chapel source tree (e.g. the homebrew
`libexec` directory) it is one make invocation:

```sh
export CHPL_HOME=$(chpl --print-chpl-home)
cd $CHPL_HOME
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=smp make -j8 runtime

# then rebuild Chaa against it:
cd <chaa>
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=smp cmake -B build-gasnet
CHPL_COMM=gasnet CHPL_COMM_SUBSTRATE=smp cmake --build build-gasnet

./build-gasnet/bin/chaa -nl 4 --problem=sedov --geometry=cartesian \
                        --nx1=256 --nx2=256 --nx3=256
```

(If the GASNet configure picked up an MPI it can't link against, add
`CHPL_GASNET_MORE_CFG_OPTIONS="--disable-mpi-compat"` to the runtime
make. On clusters use the substrate matching your network — `udp`,
`ibv`, … — and see the
[Chapel multilocale docs](https://chapel-lang.org/docs/usingchapel/multilocale.html).)

The grid (including its ghost padding) is block-partitioned across the
locales; stencil reads near partition edges are served from the
locale-local fluff cache, refreshed once per boundary application.

**Multi-locale correctness is verified**, not assumed: a 4-locale
Sedov/vortex run reproduces the single-locale fields to machine
precision (max |Δρ| = 2.6×10⁻¹⁵ after piece reassembly) and
single-locale particle trajectories to 10⁻¹²; strong-scaling overhead
of the distributed code path is ~10 % at 4 locales on one node
([benchmarks](../benchmarks.md)).

## Parallel I/O

On a single locale each dump is one file. On multiple locales Chaa
switches to **parallel piece output**:

- every locale concurrently writes its own block of the domain
  (`DInt.localSubdomain()`) as an independent HDF5 (and/or VTK) piece
  file — embarrassingly parallel, no inter-locale gather;
- locale 0 writes a single `.xmf` **XDMF spatial collection** that
  stitches the pieces back into one mesh, so ParaView/VisIt open the
  dump exactly as in the serial case;
- particle positions are gathered by id into the same single
  `*.particles.NNNN.txt` file regardless of locale count;
- the bundled python readers (`tools/chaa_io.py`,
  [Plotting & analysis](plotting.md)) reassemble piece files onto the
  global grid transparently.

This strategy works identically with serial and parallel (MPI) builds
of libhdf5 — independent files need no MPI-IO. True single-file
collective writes via parallel HDF5 (`H5Pset_fapl_mpio`) require an MPI
communicator, which Chapel's GASNet-based multilocale runtime does not
provide to the application; single-file MPI-IO output is therefore on
the roadmap behind Chapel-MPI interop, and the piece-plus-collection
scheme is the supported parallel path today.

Notes:

- physical ghost cells are owned elements of the padded domain, so
  boundary conditions are ordinary distributed foralls — no special
  cases at partition boundaries;
- the geometry is closed-form (no coordinate arrays), so metric
  evaluations are communication-free on any locale.
