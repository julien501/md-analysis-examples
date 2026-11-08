#!/usr/bin/env python3
"""Lipid-protein contact frequencies per chain and lipid type, averaged over replicates,
as tables, heatmaps and a structure with the frequency in the B-factor column."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

CMAP = LinearSegmentedColormap.from_list("white_green", ["#ffffff", "#006400"])
AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
      "HIS": "H", "HSD": "H", "HSE": "H", "HSP": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
      "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def read_frequencies(path):
    """GetContacts residue frequency file: 'chain:resname:resid  chain:resname:resid  frequency'."""
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, usecols=[0, 1, 2], names=["res1", "res2", "freq"])
    return df.set_index(["res1", "res2"]).freq


def replicate_mean(paths):
    # a pair that a replicate never formed is not in its file, so it counts as 0 there
    table = pd.concat([read_frequencies(p) for p in paths], axis=1).fillna(0.0)
    return table.mean(axis=1).rename("freq").reset_index()


def per_chain(df, lipids, chains=None, exclude_resids=()):
    """One row per (chain, residue, lipid type): the highest frequency with any single lipid molecule of that type."""
    rows = []
    for r1, r2, f in df.itertuples(index=False):
        a, b = r1.split(":"), r2.split(":")
        if a[1] in lipids:                       # make sure the protein residue comes first
            a, b = b, a
        if b[1] not in lipids or a[1] in lipids:
            continue
        chain, resname, resid = a[0], a[1], int(a[2])
        if (chains and chain not in chains) or resid in exclude_resids:
            continue
        rows.append((chain, resid, resname, b[1], f))
    out = pd.DataFrame(rows, columns=["chain", "resid", "resname", "lipid", "freq"])
    return out.groupby(["chain", "resid", "resname", "lipid"], as_index=False).freq.max()


def chain_mean(pc, chains, lipids):
    """Mean over chains per residue and lipid type. A chain without contact at a residue counts as 0."""
    wide = pc.pivot_table(index=["resid", "resname"], columns=["lipid", "chain"], values="freq", fill_value=0.0)
    return pd.DataFrame({lip: wide[lip].reindex(columns=chains, fill_value=0.0).mean(axis=1)
                         if lip in wide.columns.get_level_values(0) else 0.0 for lip in lipids})


def label(resid, resname):
    return f"{AA.get(resname[:3], 'X')}{resid}"


def heatmap(ax, table, title):
    im = ax.imshow(table.to_numpy(), cmap=CMAP, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(table.shape[1]), table.columns, rotation=90, fontsize=7)
    ax.set_title(title, fontsize=9)
    return im


def plot_chains(pc, chains, lipids, out):
    # all chains share one residue axis, so a row means the same residue in every panel
    residues = pc[["resid", "resname"]].drop_duplicates("resid").sort_values("resid")
    names = [label(r, n) for r, n in residues.itertuples(index=False)]
    fig, axes = plt.subplots(1, len(chains), figsize=(1.0 + 0.9 * len(chains), 1.0 + 0.16 * len(residues)),
                             sharey=True, squeeze=False)
    for ax, chain in zip(axes[0], chains):
        t = (pc[pc.chain == chain].pivot_table(index="resid", columns="lipid", values="freq")
             .reindex(index=residues.resid, columns=lipids).fillna(0.0))     # columns by name, not by position
        im = heatmap(ax, t, f"chain {chain}")
    axes[0, 0].set_yticks(range(len(names)), names, fontsize=6)
    fig.colorbar(im, ax=axes[0].tolist(), fraction=0.02, pad=0.01, label="contact frequency")
    fig.savefig(out, bbox_inches="tight", dpi=200)
    plt.close(fig)


def plot_mean(cm, out, top=None):
    t = cm
    if top:                       # strongest residues, still shown in sequence order
        t = cm.loc[cm.max(axis=1).nlargest(top).index].sort_index()
    t.index = [label(r, n) for r, n in t.index]
    fig, ax = plt.subplots(figsize=(1.0 + 0.5 * t.shape[1], 1.0 + 0.16 * len(t)))
    im = heatmap(ax, t, "mean over chains")
    ax.set_yticks(range(len(t)), t.index, fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.08, pad=0.04, label="contact frequency")
    fig.savefig(out, bbox_inches="tight", dpi=200)
    plt.close(fig)


def write_bfactor_pdb(structure, values, out):
    """Copy a PDB file with B-factor = 100 x frequency for the residues in values {(chain, resid): freq}, 0 elsewhere."""
    with open(structure) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                b = 100.0 * values.get((line[21], int(line[22:26])), 0.0)
                line = f"{line[:60]}{b:6.2f}{line[66:]}"     # B-factor = PDB columns 61-66
            fout.write(line)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--freq", action="append", required=True,
                   help="GetContacts residue frequency file of one replicate (repeat for each replicate)")
    p.add_argument("--lipids", required=True, help="lipid residue names, comma-separated, e.g. POPE,POPG")
    p.add_argument("--chains", default="", help="protein chains to analyse, comma-separated (default: all)")
    p.add_argument("--exclude-resids", default="",
                   help="residue numbers to leave out, e.g. a lipid-anchored residue that is always in contact")
    p.add_argument("--structure", help="PDB file to write B-factor structures from (optional)")
    p.add_argument("--top", type=int, default=25, help="residues in the short mean heatmap (default 25)")
    p.add_argument("--outdir", default="lipid_contacts_out")
    a = p.parse_args()

    lipids = a.lipids.split(",")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    excl = {int(x) for x in a.exclude_resids.split(",") if x}
    pc = per_chain(replicate_mean(a.freq), lipids, set(a.chains.split(",")) - {""}, excl)
    chains = a.chains.split(",") if a.chains else sorted(pc.chain.unique())
    cm = chain_mean(pc, chains, lipids)

    pc.to_csv(out / "per_chain.csv", index=False, float_format="%.4f")
    cm.to_csv(out / "chain_mean.csv", float_format="%.4f")
    plot_chains(pc, chains, lipids, out / "per_chain_heatmap.png")
    plot_mean(cm.copy(), out / "chain_mean_heatmap.png")
    plot_mean(cm.copy(), out / f"chain_mean_heatmap_top{a.top}.png", top=a.top)

    if a.structure:
        wide = pc.pivot_table(index=["chain", "resid"], columns="lipid", values="freq", fill_value=0.0)
        for lip in lipids + ["all"]:
            col = wide.max(axis=1) if lip == "all" else wide.get(lip)
            if col is not None:
                write_bfactor_pdb(a.structure, col.to_dict(), out / f"bfactor_{lip}.pdb")
    print(f"{len(a.freq)} replicate(s), {len(chains)} chain(s), {cm.shape[0]} residues in contact -> {out}")


if __name__ == "__main__":
    main()
