#!/usr/bin/env python3
"""Contact frequency of each system averaged over its replicas, and the difference to a reference system."""
import argparse

import pandas as pd

KEY = ["group", "position_1", "position_2"]


def replica_means(runs):
    stats = []
    for system, rows in runs.groupby("system"):
        tables = [pd.read_csv(path, sep="\t", comment="#").set_index(KEY) for path in rows.path]
        freq = pd.concat([t.frequency for t in tables], axis=1)
        names = pd.concat([t[["resname_1", "resname_2"]] for t in tables])
        # GetContacts lists only pairs that occur, so a pair missing in a replica is 0 there
        freq = freq.fillna(0.0)
        s = pd.DataFrame({"n_replicas": freq.shape[1], "mean": freq.mean(axis=1),
                          "sd": freq.std(axis=1, ddof=1) if freq.shape[1] > 1 else float("nan")})
        s = s.join(names[~names.index.duplicated()])
        stats.append(s.assign(system=system).reset_index())
    return pd.concat(stats, ignore_index=True)


def deltas(stats, reference):
    ref = stats[stats.system == reference].set_index(KEY)
    out = []
    for system, s in stats[stats.system != reference].groupby("system"):
        m = s.set_index(KEY).join(ref.add_prefix("ref_"), how="outer")
        # pair formed on one side only: frequency 0 and replica SD 0 on the other side
        m[["mean", "ref_mean"]] = m[["mean", "ref_mean"]].fillna(0.0)
        if s.n_replicas.iloc[0] > 1:
            m["sd"] = m["sd"].fillna(0.0)
        if len(ref) and ref.n_replicas.iloc[0] > 1:
            m["ref_sd"] = m["ref_sd"].fillna(0.0)
        for col in ("resname_1", "resname_2"):
            m[col] = m[col].fillna(m[f"ref_{col}"])
        m["system"], m["n_replicas"] = system, s.n_replicas.iloc[0]
        m["delta"] = m["mean"] - m["ref_mean"]
        out.append(m.reset_index())
    cols = ["system"] + KEY + ["resname_1", "resname_2", "n_replicas", "mean", "sd", "ref_mean", "ref_sd", "delta"]
    return pd.concat(out, ignore_index=True)[cols]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run_table", help="CSV with the columns system, replica, path (output of contact_frequencies.py)")
    p.add_argument("output", help="output TSV")
    p.add_argument("--reference", required=True, help="system the others are compared with")
    a = p.parse_args()

    runs = pd.read_csv(a.run_table)
    if a.reference not in set(runs.system):
        raise SystemExit(f"reference {a.reference!r} is not in {a.run_table}")
    out = deltas(replica_means(runs), a.reference)
    out.sort_values(["system"] + KEY).to_csv(a.output, sep="\t", index=False, float_format="%.4f")
    print(f"{out.system.nunique()} system(s) vs {a.reference}: {len(out)} rows -> {a.output}")


if __name__ == "__main__":
    main()
