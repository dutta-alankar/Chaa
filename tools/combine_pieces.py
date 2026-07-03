#!/usr/bin/env python3
"""combine_pieces.py — merge Chaa multi-locale piece dumps into single
global files.

Multi-locale runs write one HDF5/VTK piece per locale
(``problem.NNNN.locL.h5`` / ``.locL.vtk``) plus an XDMF collection.
This script sweeps an output directory, reassembles every snapshot
onto the global grid, and writes the same files a single-locale run
would have produced:

    problem.NNNN.h5     (+ problem.NNNN.xmf rewritten as a single grid)
    problem.NNNN.vtk

    python tools/combine_pieces.py <outdir>            # combine everything
    python tools/combine_pieces.py <outdir> --clean    # ... then delete pieces
    python tools/combine_pieces.py <outdir> --force    # overwrite existing

Notes:
- HDF5 pieces are placed by their stored native coordinates (any
  locale layout); VTK pieces are concatenated along x1, matching
  Chaa's block distribution.
- 1D ``txt`` dumps are always written globally by Chaa itself, so
  there is never anything to combine for them.
- Requires numpy (and h5py for HDF5 output).
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chaa_io import Dump


# ------------------------------- HDF5 ---------------------------------
def combine_h5(outdir, base, num, force):
    target = os.path.join(outdir, f"{base}.{num:04d}.h5")
    if os.path.exists(target) and not force:
        print(f"  {os.path.basename(target)}: exists, skipping (--force)")
        return None
    import h5py
    d = Dump(outdir, num)          # reassembles the pieces
    n1, n2, n3 = len(d.x1c), len(d.x2c), len(d.x3c)
    with h5py.File(target, "w") as f:
        for name in d.fields:
            a = d[name]
            if d.ndim == 1:
                a = a.reshape(n1)
            elif d.ndim == 2:
                a = a.reshape(n2, n1)
            f.create_dataset(name, data=a)
        f.create_dataset("cc_x1", data=d.x1c)
        f.create_dataset("cc_x2", data=d.x2c)
        f.create_dataset("cc_x3", data=d.x3c)
        f.create_dataset("node_x1", data=d.x1f)
        f.create_dataset("node_x2", data=d.x2f)
        f.create_dataset("node_x3", data=d.x3f)
        if d.nodes is not None:
            comps = ("nodes_x", "nodes_y", "nodes_z")
            for c, arr in enumerate(d.nodes):
                f.create_dataset(comps[c], data=arr)
        f.create_dataset("time", data=np.array([d.time]))
    print(f"  wrote {os.path.basename(target)}")
    return d


def write_xmf(outdir, base, num, d):
    """single-grid XDMF for the combined h5 (replaces the collection)."""
    if d.ndim < 2:
        return
    n1, n2, n3 = len(d.x1c), len(d.x2c), len(d.x3c)
    h5name = f"{base}.{num:04d}.h5"
    cell = f"{n2} {n1}" if d.ndim == 2 else f"{n3} {n2} {n1}"
    node = (f"{n2+1} {n1+1}" if d.ndim == 2
            else f"{n3+1} {n2+1} {n1+1}")
    L = []
    L.append('<?xml version="1.0" ?>')
    L.append('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>')
    L.append('<Xdmf Version="2.0">')
    L.append(' <Domain>')
    L.append('  <Grid Name="chaa" GridType="Collection"'
             ' CollectionType="Spatial">')
    L.append(f'   <Time Value="{d.time}"/>')
    L.append('  <Grid Name="mesh" GridType="Uniform">')
    dataitem = ('    <DataItem Dimensions="{dims}" NumberType="Float"'
                ' Precision="8" Format="HDF">{h5}:/{ds}</DataItem>')
    if d.nodes is None:                      # Cartesian: rectilinear
        topo = "2DRectMesh" if d.ndim == 2 else "3DRectMesh"
        geo = "VXVY" if d.ndim == 2 else "VXVYVZ"
        L.append(f'   <Topology TopologyType="{topo}"'
                 f' Dimensions="{node}"/>')
        L.append(f'   <Geometry GeometryType="{geo}">')
        L.append(dataitem.format(dims=n1 + 1, h5=h5name, ds="node_x1"))
        L.append(dataitem.format(dims=n2 + 1, h5=h5name, ds="node_x2"))
        if d.ndim == 3:
            L.append(dataitem.format(dims=n3 + 1, h5=h5name,
                                     ds="node_x3"))
        L.append('   </Geometry>')
    else:                                    # curvilinear: mapped mesh
        topo = "2DSMesh" if d.ndim == 2 else "3DSMesh"
        geo = "X_Y" if d.ndim == 2 else "X_Y_Z"
        L.append(f'   <Topology TopologyType="{topo}"'
                 f' Dimensions="{node}"/>')
        L.append(f'   <Geometry GeometryType="{geo}">')
        for c in range(2 if d.ndim == 2 else 3):
            L.append(dataitem.format(dims=node, h5=h5name,
                                     ds=("nodes_x", "nodes_y",
                                         "nodes_z")[c]))
        L.append('   </Geometry>')
    for name in d.fields:
        L.append(f'   <Attribute Name="{name}" AttributeType="Scalar"'
                 ' Center="Cell">')
        L.append(dataitem.format(dims=cell, h5=h5name, ds=name))
        L.append('   </Attribute>')
    L.append('  </Grid>')
    L.append('  </Grid>')
    L.append(' </Domain>')
    L.append('</Xdmf>')
    path = os.path.join(outdir, f"{base}.{num:04d}.xmf")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"  wrote {os.path.basename(path)}")


# -------------------------------- VTK ---------------------------------
def _parse_vtk(path):
    """Parse one legacy-ASCII Chaa VTK piece."""
    with open(path) as f:
        txt = f.read().split("\n")
    out = {"header": txt[1], "fields": {}, "forder": []}
    i = 0
    while i < len(txt):
        line = txt[i]
        if line.startswith("DATASET"):
            out["type"] = line.split()[1]
        elif line.startswith("DIMENSIONS"):
            out["dims"] = [int(v) for v in line.split()[1:]]
        elif re.match(r"[XYZ]_COORDINATES", line):
            ax, n = line.split()[0][0], int(line.split()[1])
            vals, i0 = [], i + 1
            while len(vals) < n:
                vals += [float(v) for v in txt[i0].split()]
                i0 += 1
            out[ax] = np.array(vals)
            i = i0 - 1
        elif line.startswith("POINTS"):
            n = int(line.split()[1])
            vals, i0 = [], i + 1
            while len(vals) < 3 * n:
                vals += [float(v) for v in txt[i0].split()]
                i0 += 1
            out["points"] = np.array(vals).reshape(n, 3)
            i = i0 - 1
        elif line.startswith("CELL_DATA"):
            out["ncell"] = int(line.split()[1])
        elif line.startswith("SCALARS"):
            name = line.split()[1]
            n, vals, i0 = out["ncell"], [], i + 2   # skip LOOKUP_TABLE
            while len(vals) < n:
                vals += [float(v) for v in txt[i0].split()]
                i0 += 1
            out["fields"][name] = np.array(vals)
            out["forder"].append(name)
            i = i0 - 1
        i += 1
    return out


def combine_vtk(outdir, base, num, pieces, force):
    target = os.path.join(outdir, f"{base}.{num:04d}.vtk")
    if os.path.exists(target) and not force:
        print(f"  {os.path.basename(target)}: exists, skipping (--force)")
        return
    ps = [_parse_vtk(p) for p in pieces]     # pieces sorted by locale id,
    #                                          i.e. by x1 block (Chaa
    #                                          splits only along x1)
    w = open(target, "w")
    w.write("# vtk DataFile Version 3.0\n")
    w.write(ps[0]["header"] + "\n")
    w.write("ASCII\n")

    if ps[0]["type"] == "RECTILINEAR_GRID":
        # x1 (X) coordinates concatenate, dropping the duplicated
        # shared face; Y/Z are identical across pieces
        X = np.concatenate([ps[0]["X"]] + [p["X"][1:] for p in ps[1:]])
        Y, Z = ps[0]["Y"], ps[0]["Z"]
        n1, n2, n3 = len(X), len(Y), len(Z)
        w.write("DATASET RECTILINEAR_GRID\n")
        w.write(f"DIMENSIONS {n1} {n2} {n3}\n")
        for name, arr in (("X", X), ("Y", Y), ("Z", Z)):
            w.write(f"{name}_COORDINATES {len(arr)} double\n")
            w.write("\n".join(f"{v:.9e}" for v in arr) + "\n")
        c1 = [p["dims"][0] - 1 for p in ps]
        c2, c3 = max(n2 - 1, 1), max(n3 - 1, 1)
        cells = (sum(c1), c2, c3)
    else:
        # structured grid: concatenate the point block along the x1
        # (fastest-varying) index, dropping the shared node plane
        m2, m3 = ps[0]["dims"][1], ps[0]["dims"][2]
        blocks = []
        for pn, p in enumerate(ps):
            m1 = p["dims"][0]
            pts = p["points"].reshape(m3, m2, m1, 3)
            blocks.append(pts if pn == 0 else pts[:, :, 1:, :])
        P = np.concatenate(blocks, axis=2)
        n1 = P.shape[2]
        w.write("DATASET STRUCTURED_GRID\n")
        w.write(f"DIMENSIONS {n1} {m2} {m3}\n")
        w.write(f"POINTS {n1*m2*m3} double\n")
        for row in P.reshape(-1, 3):
            w.write(f"{row[0]:.9e} {row[1]:.9e} {row[2]:.9e}\n")
        c1 = [p["dims"][0] - 1 for p in ps]
        c2 = max(m2 - 1, 1)
        c3 = max(m3 - 1, 1)
        cells = (sum(c1), c2, c3)

    w.write(f"CELL_DATA {cells[0]*cells[1]*cells[2]}\n")
    for name in ps[0]["forder"]:
        w.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
        blocks = [p["fields"][name].reshape(cells[2], cells[1], c1[i])
                  for i, p in enumerate(ps)]
        A = np.concatenate(blocks, axis=2)
        w.write("\n".join(f"{v:.9e}" for v in A.ravel()) + "\n")
    w.close()
    print(f"  wrote {os.path.basename(target)}")


# -------------------------------- main --------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", help="Chaa output directory")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing combined files")
    ap.add_argument("--clean", action="store_true",
                    help="delete the piece files after combining")
    args = ap.parse_args()

    pat = re.compile(r"(.+)\.(\d{4})\.loc(\d+)\.(h5|vtk)$")
    snaps = {}
    for f in sorted(glob.glob(os.path.join(args.outdir, "*.loc*.*"))):
        m = pat.match(os.path.basename(f))
        if m:
            key = (m.group(1), int(m.group(2)), m.group(4))
            snaps.setdefault(key, []).append((int(m.group(3)), f))
    if not snaps:
        print(f"no .locN piece files found in {args.outdir} "
              "(single-locale output is already global)")
        return

    for (base, num, ext), pieces in sorted(snaps.items()):
        pieces = [f for _, f in sorted(pieces)]
        print(f"combining {base}.{num:04d}.{ext} "
              f"({len(pieces)} pieces)")
        if ext == "h5":
            d = combine_h5(args.outdir, base, num, args.force)
            if d is not None:
                write_xmf(args.outdir, base, num, d)
        else:
            combine_vtk(args.outdir, base, num, pieces, args.force)
        if args.clean:
            for f in pieces:
                os.remove(f)
            print(f"  removed {len(pieces)} pieces")


if __name__ == "__main__":
    main()
