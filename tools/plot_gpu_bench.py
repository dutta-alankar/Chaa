#!/usr/bin/env python3
"""plot_gpu_bench.py — plot the GPU benchmark results produced by
tools/slurm/freya-gpu-bench.slurm (and multi-node GPU logs).

Parses `GPUBENCH <section> n=<N> gpus=<G> rate=<Mcell/s>` lines:
  section 'size'   -> single-GPU throughput vs problem size
  section 'strong' -> fixed total size, 1..G GPUs
  section 'weak'   -> fixed size per GPU, 1..G GPUs

Usage:
  python tools/plot_gpu_bench.py chaa-gpu-bench-*.out --save gpu.png \
      [--cpu-node <Mcell/s>]      # optional CPU-node reference line
"""
import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LINE = re.compile(
    r"GPUBENCH\s+(\w+)\s+n=(\d+)\s+gpus=(\d+)\s+rate=([\d.]+)")


def parse(paths):
    data = {"size": [], "strong": [], "weak": []}
    for path in paths:
        with open(path) as f:
            for line in f:
                m = LINE.search(line)
                if m:
                    sec, n, g, r = (m.group(1), int(m.group(2)),
                                    int(m.group(3)), float(m.group(4)))
                    if sec in data:
                        data[sec].append((n, g, r))
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--save", default="gpu-bench.png")
    ap.add_argument("--cpu-node", type=float, default=None,
                    help="CPU-node reference rate (Mcell/s)")
    args = ap.parse_args()

    d = parse(args.logs)
    npanels = sum(1 for k in ("size", "strong", "weak") if d[k])
    if npanels == 0:
        sys.exit("no GPUBENCH lines found")

    fig, axs = plt.subplots(1, npanels, figsize=(4.6*npanels, 3.8))
    if npanels == 1:
        axs = [axs]
    ip = 0

    if d["size"]:
        ax = axs[ip]; ip += 1
        pts = sorted(d["size"])
        ax.plot([p[0] for p in pts], [p[2] for p in pts], "o-",
                color="tab:green", label="1 GPU")
        if args.cpu_node:
            ax.axhline(args.cpu_node, ls="--", color="tab:gray",
                       label="full CPU node")
        ax.set_xlabel(r"$n$  (problem size $n^3$)")
        ax.set_ylabel("Mcell-updates / s")
        ax.set_title("single-GPU throughput")
        ax.legend()

    if d["strong"]:
        ax = axs[ip]; ip += 1
        pts = sorted(d["strong"], key=lambda p: p[1])
        gs = [p[1] for p in pts]
        rs = [p[2] for p in pts]
        ax.plot(gs, rs, "o-", color="tab:green", label="measured")
        ax.plot(gs, [rs[0]*g/gs[0] for g in gs], ":", color="tab:gray",
                label="ideal")
        ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
        ax.set_xticks(gs); ax.set_xticklabels([str(g) for g in gs])
        ax.set_xlabel("GPUs")
        ax.set_ylabel("Mcell-updates / s")
        ax.set_title(f"strong scaling (${pts[0][0]}^3$ total)")
        ax.legend()

    if d["weak"]:
        ax = axs[ip]; ip += 1
        pts = sorted(d["weak"], key=lambda p: p[1])
        gs = [p[1] for p in pts]
        rs = [p[2] for p in pts]
        ax.plot(gs, rs, "o-", color="tab:green", label="measured")
        ax.plot(gs, [rs[0]*g/gs[0] for g in gs], ":", color="tab:gray",
                label="ideal")
        ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
        ax.set_xticks(gs); ax.set_xticklabels([str(g) for g in gs])
        ax.set_xlabel("GPUs")
        ax.set_ylabel("Mcell-updates / s")
        ax.set_title(r"weak scaling ($256^3$ per GPU)")
        ax.legend()

    fig.tight_layout()
    fig.savefig(args.save, dpi=150)
    print("wrote", args.save)


if __name__ == "__main__":
    main()
