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
