#!/usr/bin/env python3
"""How different are the replicas of one system? Contacts and RMSF per run, replica noise per system,
and the similarity of two replicas of one system vs two runs of different systems."""
import argparse
import os
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

CUTOFF = 4.5          # heavy-atom contact distance (A)
MIN_SEQ = 4           # receptor-receptor contacts only between residues more than 4 apart
LOW, HIGH = 0.05, 0.95
FAMILIES = {"pp": "ligand-receptor contacts", "rr": "receptor-receptor contacts", "rmsf": "RMSF"}


def kabsch(mob, ref):
    mc, rc = mob.mean(axis=0), ref.mean(axis=0)
    u, _, vt = np.linalg.svd((mob - mc).T @ (ref - rc))
    d = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1, 1, d]) @ vt, mc, rc


def one_run(job):
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance, self_capped_distance

    run, cfg = job
    lo, hi = cfg["receptor_range"]
    nr = hi - lo + 1
    u = mda.Universe(run["topology"], run["trajectory"])
    sel = f"({run['receptor_selection']}) and resid {lo}:{hi}"
    rec, ca = u.select_atoms(f"{sel} and not name H*"), u.select_atoms(f"{sel} and name CA")
    if len(ca) != nr:
        raise RuntimeError(f"{run['system']} {run['replica']}: {len(ca)} C-alphas, expected {nr}")
    fit = np.isin(ca.resids, cfg["fit_resids"]) if cfg["fit_resids"] is not None else np.ones(nr, bool)
    ri = rec.resids - lo
    lig = None
    if isinstance(run.get("ligand_selection"), str) and run["ligand_selection"].strip():
        llo, lhi = cfg["ligand_range"]
        off = int(run.get("ligand_offset") or 0)     # puts the ligands of all data sets on one numbering
        lig = u.select_atoms(f"({run['ligand_selection']}) and resid {llo + off}:{lhi + off} and not name H*")
        li, nl = lig.resids - off - llo, lhi - llo + 1

    start, end = cfg["window"]
    stride = max(1, round(cfg["step_ns"] / (u.trajectory.dt / 1000.0)))
    pp = np.zeros((nl, nr)) if lig is not None else None
    rr, xs, ref = np.zeros((nr, nr)), [], None
    for ts in u.trajectory[::stride]:
        if not start <= ts.time / 1000.0 <= end + 1e-6:
            continue
        # contact matrices as flat indices, so each residue pair counts once per frame
        pairs = self_capped_distance(rec.positions, CUTOFF, box=None, return_distances=False)
        a, b = ri[pairs[:, 0]], ri[pairs[:, 1]]
        far = np.abs(a - b) > MIN_SEQ
        rr.ravel()[np.unique(np.minimum(a[far], b[far]) * nr + np.maximum(a[far], b[far]))] += 1
        if lig is not None:
            pairs = capped_distance(lig.positions, rec.positions, CUTOFF, box=None, return_distances=False)
            pp.ravel()[np.unique(li[pairs[:, 0]] * nr + ri[pairs[:, 1]])] += 1
        x = ca.positions.astype(float)
        ref = x[fit].copy() if ref is None else ref
        rot, mc, rc = kabsch(x[fit], ref)
        xs.append((x - mc) @ rot + rc)
    n = len(xs)
    if n < 10:
        raise RuntimeError(f"{run['system']} {run['replica']}: only {n} frames in {start}-{end} ns")
    xs = np.array(xs)
    out = {"rr": rr / n, "rmsf": np.sqrt(((xs - xs.mean(axis=0)) ** 2).sum(axis=2).mean(axis=0))}
    if pp is not None:
        out["pp"] = pp / n
    print(f"{run['system']} {run['replica']}: {n} frames", flush=True)
    return (run["system"], run["replica"]), out


def vec(o, fam, nr):
    return o["rr"][np.triu_indices(nr, MIN_SEQ + 1)] if fam == "rr" else o[fam].ravel()


def variable(X, fam):
    # contacts that are always on or always off say nothing about noise
    m = X.mean(axis=0)
    return np.ones(X.shape[1], bool) if fam == "rmsf" else (m >= LOW) & (m <= HIGH)


def system_noise(X, fam):
    keep = variable(X, fam)
    sd = X.std(axis=0, ddof=1) / (X.mean(axis=0) if fam == "rmsf" else 1.0)    # RMSF: relative SD
    spread = X.max(axis=0) - X.min(axis=0)
    r = np.mean([np.corrcoef(X[i, keep], X[j, keep])[0, 1] for i, j in combinations(range(len(X)), 2)])
    return {"features": int(keep.sum()), "replica_sd_median": float(np.median(sd[keep])),
            "share_range_gt_0.3": float((spread[keep] > 0.3).mean()) if fam != "rmsf" else np.nan,
            "replica_r": float(r)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run_table", help="CSV: system, replica, topology, trajectory, receptor_selection, "
                                     "ligand_selection (may be empty), ligand_offset")
    p.add_argument("outdir")
    p.add_argument("--receptor-range", nargs=2, type=int, required=True, help="receptor residues in every run")
    p.add_argument("--ligand-range", nargs=2, type=int, default=(0, -1), help="ligand residues in every run")
    p.add_argument("--fit-resids", default="", help="residues for the RMSF fit, e.g. 10-30,40-60 (default: all)")
    p.add_argument("--window", nargs=2, type=float, default=(0.0, 1e9), help="start and end in ns")
    p.add_argument("--step-ns", type=float, default=0.2, help="time between analysed frames")
    a = p.parse_args()

    fit = [r for part in a.fit_resids.split(",") if part for r in
           range(int(part.split("-")[0]), int(part.split("-")[-1]) + 1)]
    cfg = {"receptor_range": a.receptor_range, "ligand_range": a.ligand_range,
           "fit_resids": np.array(fit) if fit else None, "window": a.window, "step_ns": a.step_ns}
    nr = a.receptor_range[1] - a.receptor_range[0] + 1
    runs = pd.read_csv(a.run_table).to_dict("records")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    with Pool(max(1, min(len(runs), len(os.sched_getaffinity(0))))) as pool:
        res = dict(pool.map(one_run, [(r, cfg) for r in runs], chunksize=1))

    rows = []
    for system in dict.fromkeys(k[0] for k in res):
        reps = [o for k, o in res.items() if k[0] == system]
        for fam, name in FAMILIES.items():
            if len(reps) > 1 and all(fam in o for o in reps):
                X = np.array([vec(o, fam, nr) for o in reps])
                rows.append({"system": system, "replicas": len(reps), "descriptor": name, **system_noise(X, fam)})
    noise = pd.DataFrame(rows)
    noise.to_csv(out / "noise_per_system.csv", index=False, float_format="%.4g")

    sim = []
    for fam, name in FAMILIES.items():
        keys = [k for k in res if fam in res[k]]
        X = np.array([vec(res[k], fam, nr) for k in keys])
        C = np.corrcoef(X[:, variable(X, fam)])
        iu = np.triu_indices(len(keys), 1)
        same = np.array([[i[0] == j[0] for j in keys] for i in keys])[iu]
        if same.any() and (~same).any():
            sim.append({"descriptor": name, "r_same_system": float(np.median(C[iu][same])),
                        "r_different_systems": float(np.median(C[iu][~same]))})
    pd.DataFrame(sim).to_csv(out / "similarity.csv", index=False, float_format="%.4g")
    print(noise.round(3).to_string(index=False))
    print(pd.DataFrame(sim).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
