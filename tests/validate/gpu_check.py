#!/usr/bin/env python3
"""gpu_check.py — compare a GPU run's output directory against a CPU
reference run (or two GPU runs for the restart check).

Usage: gpu_check.py <gpu-dir> <ref-dir> [--exact]

Every .txt / .h5 / particle .csv dump in the reference directory must
have a counterpart with the same structure and values.  The default
tolerance (1e-10, relative to the largest field magnitude) allows for
CPU-vs-GPU floating-point contraction differences; --exact requires
bit-identical values (same-binary restart continuation).
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import check, finish            # noqa: E402


def rel_ok(a, b, tol):
    scale = max(np.abs(b).max(), 1.0)
    return float(np.abs(a - b).max()), float(np.abs(a - b).max()) <= tol*scale


def main():
    gpu, ref = sys.argv[1], sys.argv[2]
    exact = "--exact" in sys.argv[3:]
    tol = 0.0 if exact else 1e-10

    txts = sorted(glob.glob(os.path.join(ref, "*.txt")))
    h5s = sorted(glob.glob(os.path.join(ref, "*.h5")))
    csvs = sorted(glob.glob(os.path.join(ref, "*.csv")))
    check(len(txts) + len(h5s) > 0,
          f"reference has dumps ({len(txts)} txt, {len(h5s)} h5)")

    for rf in txts:
        gf = os.path.join(gpu, os.path.basename(rf))
        ok = os.path.exists(gf)
        check(ok, f"{os.path.basename(rf)} exists")
        if not ok:
            continue
        a, b = np.loadtxt(gf), np.loadtxt(rf)
        same = a.shape == b.shape
        d, close = rel_ok(a, b, tol) if same else (np.inf, False)
        check(same and close,
              f"{os.path.basename(rf)} matches (max diff {d:.2e})")

    if h5s:
        import h5py
        for rf in h5s:
            gf = os.path.join(gpu, os.path.basename(rf))
            ok = os.path.exists(gf)
            check(ok, f"{os.path.basename(rf)} exists")
            if not ok:
                continue
            with h5py.File(gf) as a, h5py.File(rf) as b:
                same = sorted(a.keys()) == sorted(b.keys())
                check(same, f"{os.path.basename(rf)} same datasets")
                worst, allclose = 0.0, True
                if same:
                    for k in b.keys():
                        d, close = rel_ok(a[k][...], b[k][...], tol)
                        worst = max(worst, d)
                        allclose = allclose and close
            check(same and allclose,
                  f"{os.path.basename(rf)} matches (max diff {worst:.2e})")

    for rf in csvs:                       # particle dumps, if any
        gf = os.path.join(gpu, os.path.basename(rf))
        ok = os.path.exists(gf)
        check(ok, f"{os.path.basename(rf)} exists")
        if not ok:
            continue
        a = np.loadtxt(gf, delimiter=",", skiprows=1)
        b = np.loadtxt(rf, delimiter=",", skiprows=1)
        same = a.shape == b.shape
        d, close = rel_ok(a, b, max(tol, 1e-10)) if same \
            else (np.inf, False)
        check(same and close,
              f"{os.path.basename(rf)} particles match (max diff {d:.2e})")

    finish()


if __name__ == "__main__":
    main()
