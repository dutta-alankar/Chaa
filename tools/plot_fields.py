#!/usr/bin/env python3
"""plot_fields.py — visualize the initial and final (or any) fields of a
Chaa run.

    python tools/plot_fields.py <outdir> [options]

Works on any output directory: 1D runs get line plots (initial dashed,
final solid), 2D runs get pseudocolor maps in physical coordinates
(curvilinear meshes are drawn mapped — annuli, wedges, shells), 3D runs
get slice plots (mid-plane in x3 by default; pick any axis/position
with --slice).  Tracer particles, when present, are overplotted.

Examples:
    python tools/plot_fields.py test-output/sod-1d-cart
    python tools/plot_fields.py test-output/sedov-2d-cyl --fields rho,prs --log rho
    python tools/plot_fields.py test-output/kh --dumps 0,3 --save kh.png
    python tools/plot_fields.py test-output/sedov-3d-cart --slice x1,0.5

Requires numpy, matplotlib and (for HDF5 output) h5py.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chaa_io import Dump, dump_ids, load_particles


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", help="Chaa output directory")
    ap.add_argument("--dumps", default=None,
                    help="comma-separated dump numbers (default: first,last)")
    ap.add_argument("--fields", default=None,
                    help="comma-separated field names (default: rho,prs "
                         "+ vx1 in 1D)")
    ap.add_argument("--log", default="",
                    help="comma-separated fields to plot in log10")
    ap.add_argument("--cmap", default="viridis")
    ap.add_argument("--slice", dest="slc", default="x3,0.5",
                    help="3D runs: 'x1|x2|x3[,frac]' — slice axis and "
                         "fractional position (default: x3 mid-plane)")
    ap.add_argument("--save", default=None, help="write figure to file "
                    "instead of opening a window")
    return ap.parse_args()


def field_list(args, d):
    if args.fields:
        return [f for f in args.fields.split(",") if f]
    want = ["rho", "prs", "vx1"] if d.ndim == 1 else ["rho", "prs"]
    return [f for f in want if f in d.fields]


def plot_1d(axs, dumps, fields, logs):
    for ax, f in zip(axs.flat, fields):
        for d, style in zip(dumps, ("k--", "C0-", "C1-", "C2-")):
            y = d[f].squeeze()
            y = np.log10(np.abs(y) + 1e-300) if f in logs else y
            ax.plot(d.x1c, y, style, label=f"t = {d.time:g}")
        ax.set_xlabel("x1")
        ax.set_ylabel(("log10 " if f in logs else "") + f)
        ax.legend(fontsize=8)


def slice_3d(d, f, slc):
    """Cut a 2D plane out of a 3D field: returns (z, xf, yf, xl, yl, note).

    slc is 'x1|x2|x3[,frac]' — the axis normal to the plane and the
    fractional position along it (default 0.5, the mid-plane)."""
    parts = slc.split(",")
    axis = parts[0].strip()
    frac = float(parts[1]) if len(parts) > 1 else 0.5
    n3, n2, n1 = d[f].shape
    if axis == "x1":
        i = min(int(frac * n1), n1 - 1)
        return (d[f][:, :, i], d.x2f, d.x3f, "x2", "x3",
                f"x1 = {d.x1c[i]:.3g}")
    if axis == "x2":
        j = min(int(frac * n2), n2 - 1)
        return (d[f][:, j, :], d.x1f, d.x3f, "x1", "x3",
                f"x2 = {d.x2c[j]:.3g}")
    k = min(int(frac * n3), n3 - 1)
    return (d[f][k], d.x1f, d.x2f, "x1", "x2", f"x3 = {d.x3c[k]:.3g}")


def plot_2d(fig, axs, dumps, fields, logs, cmap, outdir, slc):
    for row, f in enumerate(fields):
        for col, d in enumerate(dumps):
            ax = axs[row][col]
            note = ""
            if d.ndim == 3:
                z, xf, yf, xl, yl, note = slice_3d(d, f, slc)
                note = ", " + note
            else:
                z, xf, yf, xl, yl = d[f].squeeze(), d.x1f, d.x2f, "x1", "x2"
            if f in logs:
                z = np.log10(np.abs(z) + 1e-300)
            if d.nodes is not None and d.ndim == 2:     # curvilinear
                pc = ax.pcolormesh(d.nodes[0], d.nodes[1], z, cmap=cmap)
                ax.set_aspect("equal")
            else:
                pc = ax.pcolormesh(xf, yf, z, cmap=cmap)
                ax.set_xlabel(xl, fontsize=8)
                ax.set_ylabel(yl, fontsize=8)
            fig.colorbar(pc, ax=ax, shrink=0.85)
            p = load_particles(outdir, d.num)
            if p is not None and d.ndim == 2 and d.nodes is None:
                ax.plot(p[:, 1], p[:, 2], "w.", ms=2)
            ax.set_title(f"{('log10 ' if f in logs else '')}{f}, "
                         f"t = {d.time:g}{note}", fontsize=9)


def main():
    args = parse_args()
    if args.save:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ids = dump_ids(args.outdir)
    nums = ([int(n) for n in args.dumps.split(",")] if args.dumps
            else [ids[0], ids[-1]] if len(ids) > 1 else [ids[0]])
    dumps = [Dump(args.outdir, n) for n in nums]
    fields = field_list(args, dumps[-1])
    logs = set(args.log.split(","))

    if dumps[0].ndim == 1:
        fig, axs = plt.subplots(1, len(fields),
                                figsize=(4.2 * len(fields), 3.4),
                                squeeze=False)
        plot_1d(axs, dumps, fields, logs)
    else:
        fig, axs = plt.subplots(len(fields), len(dumps),
                                figsize=(4.8 * len(dumps), 4.0 * len(fields)),
                                squeeze=False)
        plot_2d(fig, axs, dumps, fields, logs, args.cmap, args.outdir,
                args.slc)

    fig.suptitle(os.path.basename(os.path.normpath(args.outdir)))
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"wrote {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
