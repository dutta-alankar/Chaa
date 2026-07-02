"""chaa_io.py — readers for Chaa simulation output.

Loads txt (1D), HDF5 (single-file or multi-locale ``.locN`` pieces,
which are reassembled onto the global grid) and tracer-particle dumps.

Typical use::

    from chaa_io import Dump, dump_ids, load_particles
    d = Dump("test-output/sod-1d-cart")            # last dump
    d0 = Dump("test-output/sod-1d-cart", num=0)    # initial dump
    rho = d["rho"]           # numpy array, shape (nx3, nx2, nx1)
    x = d.x1c                # cell centres
    print(d.time, d.ndim, d.fields)
"""
import glob
import os
import re

import numpy as np

_H5_COORDS = ("cc_x1", "cc_x2", "cc_x3", "node_x1", "node_x2", "node_x3",
              "nodes_x", "nodes_y", "nodes_z", "time")


def dump_ids(outdir):
    """Sorted list of dump numbers present in outdir."""
    ids = set()
    for f in glob.glob(os.path.join(outdir, "*.*")):
        if ".particles." in f:
            continue
        m = re.search(r"\.(\d{4})(\.loc\d+)?\.(txt|h5)$", f)
        if m:
            ids.add(int(m.group(1)))
    return sorted(ids)


def load_particles(outdir, num=None):
    """Particle positions of dump ``num`` (default: last) as an
    (n, 4) array of columns id, x, y, z — or None if the run had no
    particles."""
    fs = sorted(glob.glob(os.path.join(outdir, "*.particles.*.txt")))
    if not fs:
        return None
    if num is None:
        return np.loadtxt(fs[-1])
    for f in fs:
        if f.endswith(f".particles.{num:04d}.txt"):
            return np.loadtxt(f)
    return None


class Dump:
    """One Chaa output dump (fields + grid), any dimensionality.

    Fields are numpy arrays of shape (nx3, nx2, nx1); use
    ``d["rho"].squeeze()`` for the natural shape.  Coordinates are the
    native ones (x/R/r ...): ``x1c/x2c/x3c`` cell centres,
    ``x1f/x2f/x3f`` faces.  For curvilinear 2D/3D HDF5 dumps the
    physically mapped node positions are in ``nodes_x/nodes_y[/nodes_z]``
    (shape (nx2+1, nx1+1) in 2D) — use them for plotting.
    """

    def __init__(self, outdir, num=None):
        self.outdir = outdir
        ids = dump_ids(outdir)
        if not ids:
            raise FileNotFoundError(f"no Chaa dumps found in {outdir}")
        self.num = ids[-1] if num is None else num
        tag = f"{self.num:04d}"

        h5s = [f for f in sorted(glob.glob(os.path.join(outdir, f"*.{tag}*.h5")))
               if re.search(rf"\.{tag}(\.loc\d+)?\.h5$", f)]
        txts = [f for f in sorted(glob.glob(os.path.join(outdir, f"*.{tag}.txt")))
                if ".particles." not in f]
        if h5s:
            self._load_h5(h5s)
        elif txts:
            self._load_txt(txts[0])
        else:
            raise FileNotFoundError(f"dump {self.num} not found in {outdir}")

    # ------------------------------------------------------------------
    def _load_txt(self, path):
        with open(path) as f:
            header = f.readline()
            names = f.readline().lstrip("#").split()
        m = re.search(r"time=(\S+)", header)
        self.time = float(m.group(1)) if m else 0.0
        d = np.loadtxt(path)
        self.x1c = d[:, 0]
        self.x2c = np.zeros(1)
        self.x3c = np.zeros(1)
        self.x1f = self.x2f = self.x3f = None
        self.nodes = None
        self.ndim = 1
        self.data = {n: d[:, i + 1].reshape(1, 1, -1)
                     for i, n in enumerate(names[1:])}

    def _load_h5(self, paths):
        import h5py
        pieces = []
        for p in paths:
            with h5py.File(p, "r") as f:
                pieces.append({k: f[k][...] for k in f.keys()})
        p0 = pieces[0]
        self.time = float(p0["time"][0])

        # global native coordinate axes (pieces overlap-free by design)
        cc = [np.unique(np.concatenate([q[f"cc_x{a}"] for q in pieces]))
              for a in (1, 2, 3)]
        self.x1c, self.x2c, self.x3c = cc
        n1, n2, n3 = len(cc[0]), len(cc[1]), len(cc[2])
        self.ndim = 1 + (n2 > 1) + (n3 > 1)

        fields = [k for k in p0 if k not in _H5_COORDS]
        self.data = {k: np.empty((n3, n2, n1)) for k in fields}
        self.x1f = np.empty(n1 + 1)
        self.x2f = np.empty(n2 + 1)
        self.x3f = np.empty(n3 + 1)
        curvi = "nodes_x" in p0
        self.nodes = None
        if curvi:
            nshape = (n3 + 1, n2 + 1, n1 + 1) if self.ndim == 3 \
                else (n2 + 1, n1 + 1)
            self.nodes = [np.empty(nshape)
                          for _ in range(3 if self.ndim == 3 else 2)]

        for q in pieces:
            i0 = int(np.searchsorted(cc[0], q["cc_x1"][0]))
            j0 = int(np.searchsorted(cc[1], q["cc_x2"][0]))
            k0 = int(np.searchsorted(cc[2], q["cc_x3"][0]))
            c1, c2, c3 = len(q["cc_x1"]), len(q["cc_x2"]), len(q["cc_x3"])
            for k in fields:
                self.data[k][k0:k0 + c3, j0:j0 + c2, i0:i0 + c1] = \
                    q[k].reshape(c3, c2, c1)
            self.x1f[i0:i0 + c1 + 1] = q["node_x1"]
            self.x2f[j0:j0 + c2 + 1] = q["node_x2"]
            self.x3f[k0:k0 + c3 + 1] = q["node_x3"]
            if curvi:
                nx = q["nodes_x"]
                sl = (slice(k0, k0 + c3 + 1), slice(j0, j0 + c2 + 1),
                      slice(i0, i0 + c1 + 1)) if self.ndim == 3 else \
                     (slice(j0, j0 + c2 + 1), slice(i0, i0 + c1 + 1))
                self.nodes[0][sl] = nx
                self.nodes[1][sl] = q["nodes_y"]
                if self.ndim == 3:
                    self.nodes[2][sl] = q["nodes_z"]

    # ------------------------------------------------------------------
    @property
    def fields(self):
        return sorted(self.data.keys())

    def __getitem__(self, name):
        return self.data[name]

    def radius(self, cen=(0.0, 0.0, 0.0)):
        """Cell-centre distance from ``cen`` treating the native
        coordinates as Cartesian (correct for Cartesian runs; for a 1D
        spherical/cylindrical run x1 *is* the radius)."""
        X3, X2, X1 = np.meshgrid(self.x3c, self.x2c, self.x1c,
                                 indexing="ij")
        return np.sqrt((X1 - cen[0]) ** 2 + (X2 - cen[1]) ** 2 +
                       (X3 - cen[2]) ** 2)
