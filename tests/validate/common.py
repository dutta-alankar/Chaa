"""Shared helpers for Chaa test validation."""
import sys
import numpy as np

_failures = []


def load_txt(path):
    """1D txt dump -> columns x1, rho, vx1, vx2, vx3, prs."""
    return np.loadtxt(path)


def load_h5(path):
    import h5py
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            out[k] = f[k][...]
    return out


def check(cond, msg):
    tag = "PASS" if cond else "FAIL"
    print(f"  [{tag}] {msg}")
    if not cond:
        _failures.append(msg)


def finish():
    if _failures:
        print(f"{len(_failures)} check(s) failed")
        sys.exit(1)
    print("all checks passed")
    sys.exit(0)


def sedov_radius(t, gamma=1.4, E0=1.0, rho0=1.0):
    """Spherical Sedov-Taylor shock radius R = (E t^2 / (alpha rho))^(1/5).

    alpha values from Kamm & Timmes (2007), spherical symmetry.
    """
    alpha = {1.4: 0.851072, 1.666667: 0.493610}[round(gamma, 6)]
    return (E0 * t * t / (alpha * rho0)) ** 0.2
