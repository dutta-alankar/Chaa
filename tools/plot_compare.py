#!/usr/bin/env python3
"""plot_compare.py — plot a Chaa run against the analytic estimate.

    python tools/plot_compare.py <kind> <outdir> [options]

kinds and the analytic reference used:
    sod             exact Riemann solution (Toro sampling)
    sod-iso         exact isothermal Riemann solution      [--cs]
    sedov           Sedov-Taylor similarity shock radius   [--gamma --e0]
    taylor-couette  analytic steady Couette profile        [--omega-in ...]
    thermal-wave    conduction decay  exp(-kappa (g-1)/g k^2 t)  [--kappa]
    cooling         exact Townsend power-law cooling       [--lambda0 --alpha]
    linear-wave     the initial eigenmode (returns after each period)
    vortex          the initial vortex (exact solution after one period)
    epicycle        epicyclic oscillation <vx> = -A cos(kappa_ep t)

Examples:
    python tools/plot_compare.py sod test-output/sod-1d-cart
    python tools/plot_compare.py sedov test-output/sedov-1d-sph
    python tools/plot_compare.py epicycle out-epicycle --amp 0.01

Requires numpy, matplotlib (and h5py for HDF5 output).
"""
import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "tests", "validate"))
from chaa_io import Dump, dump_ids  # noqa: E402
import exact_riemann                # noqa: E402
from common import sedov_radius     # noqa: E402


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["sod", "sod-iso", "sedov",
                                     "taylor-couette", "thermal-wave",
                                     "cooling", "linear-wave", "vortex",
                                     "epicycle"])
    ap.add_argument("outdir")
    ap.add_argument("--save", default=None)
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--e0", type=float, default=1.0)
    ap.add_argument("--cen", default="0,0,0", help="blast centre (sedov)")
    ap.add_argument("--cs", type=float, default=1.0)
    ap.add_argument("--kappa", type=float, default=0.02)
    ap.add_argument("--lambda0", type=float, default=0.1)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--omega-in", type=float, default=1.0)
    ap.add_argument("--omega-out", type=float, default=0.0)
    ap.add_argument("--r-in", type=float, default=1.0)
    ap.add_argument("--r-out", type=float, default=2.0)
    ap.add_argument("--omega", type=float, default=1.0, help="epicycle Omega")
    ap.add_argument("--q", type=float, default=1.5, help="epicycle shear q")
    ap.add_argument("--amp", type=float, default=0.01, help="epicycle amp")
    return ap.parse_args()


def l1(a, b):
    return np.abs(np.asarray(a) - np.asarray(b)).mean()


def cmp_sod(args, d, ax):
    x = d.x1c
    re_, ue, pe = exact_riemann.solution(x, t=d.time, gamma=args.gamma)
    ax.plot(x, d["rho"].squeeze(), "C0.", ms=3, label="Chaa rho")
    ax.plot(x, re_, "k-", lw=1, label="exact")
    ax.plot(x, d["prs"].squeeze(), "C1.", ms=3, label="Chaa p")
    ax.plot(x, pe, "k--", lw=1)
    ax.set_xlabel("x")
    ax.set_title(f"Sod, t={d.time:g}: L1(rho) = "
                 f"{l1(d['rho'].squeeze(), re_):.2e}")
    ax.legend()


def cmp_sod_iso(args, d, ax):
    x = d.x1c
    re_, ue = exact_riemann.iso_solution(x, t=d.time, cs=args.cs)
    ax.plot(x, d["rho"].squeeze(), "C0.", ms=3, label="Chaa rho")
    ax.plot(x, re_, "k-", lw=1, label="exact (isothermal)")
    ax.set_xlabel("x")
    ax.set_title(f"isothermal Sod, t={d.time:g}: L1(rho) = "
                 f"{l1(d['rho'].squeeze(), re_):.2e}")
    ax.legend()


def cmp_sedov(args, d, ax):
    cen = [float(c) for c in args.cen.split(",")]
    r = d.radius(cen).ravel() if d.ndim > 1 else d.x1c
    rho = d["rho"].ravel() if d.ndim > 1 else d["rho"].squeeze()
    rs = sedov_radius(d.time, gamma=round(args.gamma, 6), E0=args.e0)
    ax.plot(r, rho, "C0.", ms=1.5, alpha=0.6, label="Chaa cells")
    ax.axvline(rs, color="k", ls="--",
               label=f"similarity shock radius {rs:.3f}")
    peak = r[np.argmax(rho)]
    ax.set_xlabel("radius")
    ax.set_ylabel("rho")
    ax.set_title(f"Sedov, t={d.time:g}: peak at r={peak:.3f} "
                 f"({100 * abs(peak - rs) / rs:.1f}% off)")
    ax.legend()


def cmp_taylor_couette(args, d, ax):
    R = d.x1c
    R1, R2 = args.r_in, args.r_out
    a = (args.omega_out * R2**2 - args.omega_in * R1**2) / (R2**2 - R1**2)
    b = (args.omega_in - args.omega_out) * R1**2 * R2**2 / (R2**2 - R1**2)
    ana = a * R + b / R
    vphi = d["vx3"].squeeze()
    ax.plot(R, vphi, "C0.", ms=4, label="Chaa v_phi")
    ax.plot(R, ana, "k-", lw=1, label="analytic Couette")
    ax.set_xlabel("R")
    ax.set_title(f"Taylor-Couette: rel. L1 = "
                 f"{l1(vphi, ana) / np.abs(ana).mean():.2e}")
    ax.legend()


def cmp_thermal_wave(args, d0, d, ax):
    g, kw = args.gamma, 2 * np.pi / (d.x1f[-1] - d.x1f[0]
                                     if d.x1f is not None else 1.0)
    a0 = np.abs(d0["rho"].squeeze() - 1.0).max()
    dec = np.exp(-args.kappa * (g - 1) / g * kw**2 * (d.time - d0.time))
    ax.plot(d0.x1c, d0["rho"].squeeze() - 1.0, "k--", label=f"t={d0.time:g}")
    ax.plot(d.x1c, d["rho"].squeeze() - 1.0, "C0-", label=f"t={d.time:g}")
    ax.plot(d0.x1c, (d0["rho"].squeeze() - 1.0) * dec, "C3:",
            label="analytic decay")
    a1 = np.abs(d["rho"].squeeze() - 1.0).max()
    ax.set_xlabel("x")
    ax.set_ylabel("drho")
    ax.set_title(f"conduction decay {a1/a0:.4f} vs analytic {dec:.4f}")
    ax.legend()


def cmp_cooling(args, outdir, ax):
    ts, Ts = [], []
    for n in dump_ids(outdir):
        d = Dump(outdir, n)
        ts.append(d.time)
        Ts.append((d["prs"] / d["rho"]).mean())
    ts, Ts = np.array(ts), np.array(Ts)
    a, g = args.alpha, args.gamma
    ana = (1.0 - (1 - a) * (g - 1) * args.lambda0 * ts) ** (1 / (1 - a))
    ax.plot(ts, Ts, "C0o", label="Chaa <T>")
    tt = np.linspace(ts[0], ts[-1], 200)
    ax.plot(tt, (1.0 - (1 - a) * (g - 1) * args.lambda0 * tt) ** (1 / (1 - a)),
            "k-", lw=1, label="exact Townsend")
    ax.set_xlabel("t")
    ax.set_ylabel("T")
    ax.set_title(f"power-law cooling: max rel. err "
                 f"{np.abs(Ts / ana - 1).max():.2e}")
    ax.legend()


def cmp_periodic_return(args, d0, d, ax, label):
    ax.plot(d0.x1c, d0["rho"].squeeze() if d.ndim == 1
            else d0["rho"].mean(axis=(0, 1)), "k--", label=f"t={d0.time:g}")
    ax.plot(d.x1c, d["rho"].squeeze() if d.ndim == 1
            else d["rho"].mean(axis=(0, 1)), "C0-", label=f"t={d.time:g}")
    ax.set_xlabel("x")
    ax.set_ylabel("rho" + ("" if d.ndim == 1 else " (y-mean)"))
    ax.set_title(f"{label}: L1 vs initial = "
                 f"{l1(d['rho'], d0['rho']):.2e}")
    ax.legend()


def cmp_epicycle(args, outdir, ax):
    ts, vxs = [], []
    for n in dump_ids(outdir):
        d = Dump(outdir, n)
        ts.append(d.time)
        vxs.append(d["vx1"].mean())
    kap = np.sqrt(2 * (2 - args.q)) * args.omega
    tt = np.linspace(min(ts), max(ts), 300)
    ax.plot(ts, vxs, "C0o", label="Chaa <vx>")
    ax.plot(tt, args.amp * np.cos(kap * tt), "k-", lw=1,
            label=f"A cos(kappa t), kappa={kap:g}")
    ax.set_xlabel("t")
    ax.set_ylabel("<vx>")
    ax.set_title("epicyclic oscillation (shearing box)")
    ax.legend()


def main():
    args = parse_args()
    if args.save:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.6))
    ids = dump_ids(args.outdir)
    if args.kind == "sod":
        cmp_sod(args, Dump(args.outdir), ax)
    elif args.kind == "sod-iso":
        cmp_sod_iso(args, Dump(args.outdir), ax)
    elif args.kind == "sedov":
        cmp_sedov(args, Dump(args.outdir), ax)
    elif args.kind == "taylor-couette":
        cmp_taylor_couette(args, Dump(args.outdir), ax)
    elif args.kind == "thermal-wave":
        cmp_thermal_wave(args, Dump(args.outdir, ids[0]), Dump(args.outdir), ax)
    elif args.kind == "cooling":
        cmp_cooling(args, args.outdir, ax)
    elif args.kind == "linear-wave":
        cmp_periodic_return(args, Dump(args.outdir, ids[0]),
                            Dump(args.outdir), ax, "linear wave")
    elif args.kind == "vortex":
        cmp_periodic_return(args, Dump(args.outdir, ids[0]),
                            Dump(args.outdir), ax, "isentropic vortex")
    elif args.kind == "epicycle":
        cmp_epicycle(args, args.outdir, ax)

    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"wrote {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
