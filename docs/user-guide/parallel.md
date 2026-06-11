# Running in parallel

Chaa's fields live on a single block-distributed domain
(`StencilDist`) with halo ("fluff") regions. Every loop in the solver
is a data-parallel `forall`; halo refreshes happen in one call,
`updateFluff()`. There is no other communication logic anywhere in the
code — Chapel moves the data.

## Shared memory (default)

A standard build (`CHPL_COMM=none`) uses all cores of the node
automatically:

```sh
CHPL_RT_NUM_THREADS_PER_LOCALE=8 ./build/bin/chaa --problem=sedov ...
```

## Multi-locale (distributed)

Build Chapel with a communication layer (typically
`CHPL_COMM=gasnet`; see the
[Chapel multilocale docs](https://chapel-lang.org/docs/usingchapel/multilocale.html)),
then rebuild Chaa and launch with `-nl`:

```sh
./build/bin/chaa -nl 8 --problem=sedov --geometry=cartesian \
                 --nx1=512 --nx2=512 --nx3=512
```

The grid (including its ghost padding) is block-partitioned across the
8 locales; stencil reads near partition edges are served from the
locale-local fluff cache, refreshed once per boundary application.

## Parallel I/O

On a single locale each dump is one file. On multiple locales Chaa
switches to **parallel piece output**:

- every locale concurrently writes its own block of the domain
  (`DInt.localSubdomain()`) as an independent HDF5 (and/or VTK) piece
  file — embarrassingly parallel, no inter-locale gather;
- locale 0 writes a single `.xmf` **XDMF spatial collection** that
  stitches the pieces back into one mesh, so ParaView/VisIt open the
  dump exactly as in the serial case.

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
