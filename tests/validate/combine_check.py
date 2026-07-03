#!/usr/bin/env python3
"""combine_check.py — compare a combined multi-locale output directory
against a single-locale reference run.

Usage: combine_check.py <combined-dir> <reference-dir>

Requires every reference .h5/.vtk/.xmf dump to have a combined
counterpart with the same structure and values (to multi- vs
single-locale round-off), and no leftover .locN piece files.
"""
import glob
import os
import sys
import xml.dom.minidom

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "tools"))
from common import check, finish            # noqa: E402
from combine_pieces import _parse_vtk       # noqa: E402


def main():
    cmb, ref = sys.argv[1], sys.argv[2]

    h5s = sorted(glob.glob(os.path.join(ref, "*.h5")))
    check(len(h5s) > 0, f"reference has h5 dumps ({len(h5s)})")
    import h5py
    for rf in h5s:
        cf = os.path.join(cmb, os.path.basename(rf))
        ok = os.path.exists(cf)
        check(ok, f"{os.path.basename(rf)} combined file exists")
        if not ok:
            continue
        with h5py.File(cf) as a, h5py.File(rf) as b:
            same_keys = sorted(a.keys()) == sorted(b.keys())
            check(same_keys, f"{os.path.basename(rf)} same datasets")
            worst = max(float(np.abs(a[k][...] - b[k][...]).max())
                        for k in b.keys()) if same_keys else np.inf
        check(worst < 1e-11,
              f"{os.path.basename(rf)} values match (max diff {worst:.2e})")

    for rf in sorted(glob.glob(os.path.join(ref, "*.vtk"))):
        cf = os.path.join(cmb, os.path.basename(rf))
        ok = os.path.exists(cf)
        check(ok, f"{os.path.basename(rf)} combined file exists")
        if not ok:
            continue
        va, vb = _parse_vtk(cf), _parse_vtk(rf)
        check(va["dims"] == vb["dims"],
              f"{os.path.basename(rf)} same dimensions")
        geo = 0.0
        if "points" in vb:
            geo = float(np.abs(va["points"] - vb["points"]).max())
        else:
            geo = max(float(np.abs(va[a_] - vb[a_]).max())
                      for a_ in "XYZ")
        fld = max(float(np.abs(va["fields"][k] - vb["fields"][k]).max())
                  for k in vb["fields"])
        check(geo < 1e-11 and fld < 1e-11,
              f"{os.path.basename(rf)} geometry/fields match "
              f"({geo:.2e}/{fld:.2e})")

    for xf in sorted(glob.glob(os.path.join(cmb, "*.xmf"))):
        try:
            xml.dom.minidom.parse(xf)
            check(True, f"{os.path.basename(xf)} well-formed")
        except Exception as e:
            check(False, f"{os.path.basename(xf)} well-formed ({e})")

    left = glob.glob(os.path.join(cmb, "*.loc*"))
    check(len(left) == 0, "no piece files left after --clean")
    finish()


if __name__ == "__main__":
    main()
