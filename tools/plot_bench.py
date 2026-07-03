#!/usr/bin/env python3
"""plot_bench.py — plot Chaa scaling-benchmark results.

Parses the result lines emitted by tools/slurm/freya-bench-node.slurm
and tools/slurm/freya-bench-multi.slurm:

    NODE   8 threads   128x 128x 128  4.95 Mcell/s
    MULTI  4 locales   256x 256x 256  61.2 Mcell/s

as well as the markdown tables printed by tools/bench.sh:

    | 1 | 8 | 128x 128x 128 | 5.93 |     (locales | threads | grid | rate)

under "strong scaling" / "weak scaling" section headers, and produces
a two-panel figure: measured throughput with the ideal-scaling line
(strong), and weak-scaling efficiency. For bench.sh output, --select
single|multi picks which pair of sections to plot; the fixed-total-
threads multi-locale protocol gets a flat ideal line.

    python tools/plot_bench.py chaa-bench-node-*.out --save node.png
    python tools/plot_bench.py bench.out --select single --save laptop.png
"""
import argparse
import re
import sys

import numpy as np


def parse(path, select=None):
    """-> (dict with 'strong'/'weak' lists of (units, rate), kind,
    flat_ideal) — flat_ideal marks a fixed-total-resources strong
    protocol whose ideal line is horizontal."""
    out = {"strong": [], "weak": []}
    sec = None
    keep = True
    flat = False
    pat = re.compile(r"^(NODE|MULTI)\s+(\d+)\s+(?:threads|locales)\s+"
                     r"(\d+)x\s*(\d+)x\s*(\d+)\s+([0-9.]+)\s+Mcell/s")
    tab = re.compile(r"^\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*"
                     r"(\d+)x\s*(\d+)x\s*(\d+)\s*\|\s*([0-9.]+)\s*\|")
    kind = None
    for line in open(path):
        if "strong scaling" in line or "weak scaling" in line:
            sec = "strong" if "strong" in line else "weak"
            keep = select is None or select in line
            if sec == "strong" and keep:
                flat = "fixed total" in line
                kind = "locales" if "multi" in line.lower() else kind
        if not keep:
            continue
        m = pat.match(line.strip())
        if m and sec:
            kind = "threads" if m.group(1) == "NODE" else "locales"
            out[sec].append((int(m.group(2)), float(m.group(6))))
            continue
        m = tab.match(line.strip())
        if m and sec:
            locales, threads = int(m.group(1)), int(m.group(2))
            if kind is None:
                kind = "locales" if locales > 1 else "threads"
            out[sec].append((locales if kind == "locales" else threads,
                             float(m.group(6))))
    return out, kind, flat


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="benchmark output file")
    ap.add_argument("--save", default=None)
    ap.add_argument("--title", default=None)
    ap.add_argument("--select", default=None,
                    help="only sections whose header contains this "
                         "(e.g. 'single-locale' / 'multi-locale' for "
                         "bench.sh output)")
    args = ap.parse_args()
    if args.save:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    data, kind, flat = parse(args.results, args.select)
    if not data["strong"] and not data["weak"]:
        sys.exit(f"no benchmark lines found in {args.results}")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    if data["strong"]:
        n = np.array([u for u, _ in data["strong"]], dtype=float)
        r = np.array([v for _, v in data["strong"]])
        a1.loglog(n, r, "o-", label="measured")
        ideal = np.full_like(n, r[0]) if flat else r[0] * n / n[0]
        a1.loglog(n, ideal, "k--", lw=1, label="ideal")
        a1.set_xlabel(kind)
        a1.set_ylabel("Mcell-updates/s")
        a1.set_title("strong scaling (fixed problem"
                     + (", fixed total threads)" if flat else ")"),
                     fontsize=10)
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
