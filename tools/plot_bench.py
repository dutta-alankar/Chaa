#!/usr/bin/env python3
"""plot_bench.py — plot Chaa scaling-benchmark results.

Parses the result lines emitted by tools/slurm/freya-bench-node.slurm
and tools/slurm/freya-bench-multi.slurm (or anything printing the same
format):

    NODE   8 threads   128x 128x 128  4.95 Mcell/s
    MULTI  4 locales   256x 256x 256  61.2 Mcell/s

under "== strong scaling" / "== weak scaling" section headers, and
produces a two-panel figure per input file: measured throughput with
the ideal-scaling line (strong), and weak-scaling efficiency.

    python tools/plot_bench.py chaa-bench-node-*.out --save node.png
    python tools/plot_bench.py chaa-bench-multi-*.out --save multi.png
"""
import argparse
import re
import sys

import numpy as np


def parse(path):
    """-> dict with 'strong' and 'weak' lists of (units, rate)."""
    out = {"strong": [], "weak": []}
    sec = None
    pat = re.compile(r"^(NODE|MULTI)\s+(\d+)\s+(?:threads|locales)\s+"
                     r"(\d+)x\s*(\d+)x\s*(\d+)\s+([0-9.]+)\s+Mcell/s")
    kind = None
    for line in open(path):
        if "strong scaling" in line:
            sec = "strong"
        elif "weak scaling" in line:
            sec = "weak"
        m = pat.match(line.strip())
        if m and sec:
            kind = "threads" if m.group(1) == "NODE" else "locales"
            out[sec].append((int(m.group(2)), float(m.group(6))))
    return out, kind


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="benchmark output file")
    ap.add_argument("--save", default=None)
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    if args.save:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    data, kind = parse(args.results)
    if not data["strong"] and not data["weak"]:
        sys.exit(f"no benchmark lines found in {args.results}")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    if data["strong"]:
        n = np.array([u for u, _ in data["strong"]], dtype=float)
        r = np.array([v for _, v in data["strong"]])
        a1.loglog(n, r, "o-", label="measured")
        a1.loglog(n, r[0] * n / n[0], "k--", lw=1, label="ideal")
        a1.set_xlabel(kind)
        a1.set_ylabel("Mcell-updates/s")
        a1.set_title("strong scaling (fixed problem)", fontsize=10)
        a1.set_xticks(n)
        a1.set_xticklabels([f"{int(x)}" for x in n])
        a1.xaxis.set_minor_formatter(NullFormatter())
        a1.legend(fontsize=9)

    if data["weak"]:
        n = np.array([u for u, _ in data["weak"]], dtype=float)
        r = np.array([v for _, v in data["weak"]])
        eff = 100.0 * (r / n) / (r[0] / n[0])
        a2.semilogx(n, eff, "s-", color="C1", label="measured")
        a2.axhline(100.0, color="k", ls="--", lw=1, label="ideal")
        a2.set_xlabel(kind)
        a2.set_ylabel("weak-scaling efficiency [%]")
        a2.set_title("weak scaling (fixed work per unit)", fontsize=10)
        a2.set_ylim(0, 115)
        a2.set_xticks(n)
        a2.set_xticklabels([f"{int(x)}" for x in n])
        a2.xaxis.set_minor_formatter(NullFormatter())
        a2.legend(fontsize=9)

    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"wrote {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
