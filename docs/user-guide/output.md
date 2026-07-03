# Output & visualisation

Select formats with a comma list, and a dump cadence:

```sh
--outFormats=txt,vtk,hdf5  --outDt=0.1  --outDir=results
```

The initial and final states are always written; with `--outDt>0` you
get a time series `problem.0000`, `problem.0001`, …

## txt (1D)

Plain columns `x1 rho vx1 vx2 vx3 prs` with a `#` header carrying the
time stamp. Ideal for quick plots and the 1D validators.

## VTK

Legacy ASCII VTK, one file per dump:

- **Cartesian** runs are written as `RECTILINEAR_GRID`;
- **curvilinear** runs as `STRUCTURED_GRID` with node positions mapped
  to physical space — polar annuli, spherical wedges and shells render
  with their true shape in ParaView/VisIt.

Fields are cell data: `rho`, `vx1`, `vx2`, `vx3`, `prs`.

## HDF5 + XDMF

Each dump writes `problem.NNNN.h5` and (for 2D/3D) a matching
`problem.NNNN.xmf`. Open the `.xmf` in ParaView or VisIt — it
references the heavy data in the `.h5`.

Inside the HDF5 file:

| dataset | content |
|---|---|
| `rho`, `vx1`, `vx2`, `vx3`, `prs` | cell-centred fields, C order with x1 fastest (shape `[nx3][nx2][nx1]`, squeezed for 1D/2D) |
| `cc_x1`, `cc_x2`, `cc_x3` | cell-centre coordinates (1D arrays) |
| `node_x1`, `node_x2`, `node_x3` | face/node coordinates (1D arrays) |
| `nodes_x`, `nodes_y`[, `nodes_z`] | physically mapped node positions (curvilinear runs) |
| `time` | dump time |

```python
import h5py
with h5py.File("results/blast.0003.h5") as f:
    rho = f["rho"][:]
    t = f["time"][0]
```

!!! note
    HDF5 support is compiled in by default (`-DCHAA_HDF5=ON`) and links
    against `libhdf5`. A binary built with `-DCHAA_HDF5=OFF` reports a
    clear error if `hdf5` output is requested.

## Reading dumps back

Beyond ParaView/VisIt, the bundled python tools read every format
(including multi-locale piece files, reassembled transparently):
`tools/plot_fields.py` for quick looks, `tools/plot_compare.py` for
overlays on analytic solutions, and `tools/chaa_io.py` as an importable
reader — see [Plotting & analysis](plotting.md).

## Combining multi-locale pieces into single files

Multi-locale runs write one HDF5/VTK piece per locale
(`problem.NNNN.locL.h5` / `.locL.vtk`) plus an XDMF collection that
stitches them for ParaView/VisIt. If you want *single global files*
instead — for archiving, for tools that don't read XDMF collections,
or just for tidiness — sweep the output directory with

```sh
python tools/combine_pieces.py output/            # every snapshot
python tools/combine_pieces.py output/ --clean    # ...then delete pieces
python tools/combine_pieces.py output/ --force    # overwrite existing
```

For each snapshot this produces exactly what a single-locale run would
have written: one `problem.NNNN.h5` (with the `.xmf` rewritten as a
single grid), and one `problem.NNNN.vtk` — validated to match genuine
single-locale output to machine precision for Cartesian (rectilinear)
and curvilinear (mapped structured-grid) meshes alike. HDF5 pieces are
placed by their stored native coordinates; 1D `txt` dumps are always
written globally by Chaa itself, so there is never anything to combine
for them.
