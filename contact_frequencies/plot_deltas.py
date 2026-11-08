#!/usr/bin/env python3
"""Heatmap per system of the contact frequency differences between two chains (output of frequency_deltas.py)."""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

CMAP = LinearSegmentedColormap.from_list("delta", ["#2a78d6", "#f5f4ef", "#e34948"])   # less / same / more


def between_chains(df, rows, cols):
    c1, r1 = df.position_1.str.split(":", expand=True).T.values
    c2, r2 = df.position_2.str.split(":", expand=True).T.values
    fwd, rev = (c1 == rows) & (c2 == cols), (c1 == cols) & (c2 == rows)
    out = pd.DataFrame({"system": df.system, "row": np.where(fwd, r1, r2), "col": np.where(fwd, r2, r1),
                        "delta": df.delta})[fwd | rev]
    return out.astype({"row": int, "col": int})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("deltas")
    p.add_argument("output", help=".svg, .png or .pdf")
    p.add_argument("--rows", required=True, help="chain on the rows")
    p.add_argument("--cols", required=True, help="chain on the columns")
    p.add_argument("--group", default="all", help="interaction group (default: all)")
    p.add_argument("--row-offset", type=int, default=0, help="added to the row residue numbers in the labels")
    p.add_argument("--ncols", type=int, default=4, help="panels per row")
    a = p.parse_args()

    df = between_chains(pd.read_csv(a.deltas, sep="\t").query("group == @a.group"), a.rows, a.cols)
    if df.empty:
        raise SystemExit(f"no {a.rows}-{a.cols} pairs in group {a.group!r}")
    row_ids = np.arange(df.row.min(), df.row.max() + 1)
    col_ids = np.sort(df.col.unique())          # only residues in contact somewhere, to keep the panels readable
    systems = sorted(df.system.unique())
    vmax = max(df.delta.abs().max(), 1e-6)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)   # same scale in every panel

    ncols = min(a.ncols, len(systems))
    nrows = -(-len(systems) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 2.8 * nrows), squeeze=False)
    axes = axes.ravel()
    for ax, system in zip(axes, systems):
        d = df[df.system == system]
        m = np.zeros((len(row_ids), len(col_ids)))
        m[d.row - row_ids[0], np.searchsorted(col_ids, d.col)] = d.delta
        im = ax.imshow(m, cmap=CMAP, norm=norm, aspect="auto", interpolation="nearest")
        ax.set_title(system, fontsize=10)
        xt = np.arange(0, len(col_ids), max(1, len(col_ids) // 8))
        ax.set_xticks(xt, col_ids[xt], fontsize=6, rotation=90)
        yt = [i for i, r in enumerate(row_ids) if (r + a.row_offset) % 5 == 0]
        ax.set_yticks(yt, [row_ids[i] + a.row_offset for i in yt], fontsize=6)
        ax.set_xlabel(f"chain {a.cols}", fontsize=7)
        ax.set_ylabel(f"chain {a.rows}", fontsize=7)
    for ax in axes[len(systems):]:
        ax.axis("off")
    fig.colorbar(im, ax=axes[:len(systems)].tolist(), fraction=0.02, pad=0.02).set_label(f"Δ frequency ({a.group})")
    fig.savefig(a.output, bbox_inches="tight", dpi=200)
    print(f"{len(systems)} panel(s) -> {a.output}")


if __name__ == "__main__":
    main()
