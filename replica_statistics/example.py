#!/usr/bin/env python3
"""replica_stats on simulated numbers: six systems with three replicas each, only V3 differs from WT."""
import numpy as np

import replica_stats as rs

rng = np.random.default_rng(1)
systems = ["WT", "V1", "V2", "V3", "V4", "V5"]
true_mean = {"WT": 20.0, "V1": 20.2, "V2": 19.9, "V3": 22.5, "V4": 20.1, "V5": 20.3}
groups = np.repeat(systems, 3)
values = np.array([rng.normal(true_mean[s], 0.8) for s in groups])     # one value per replica
assay = {"WT": 100, "V1": 95, "V2": 104, "V3": 60, "V4": 98, "V5": 90}

for s in systems:
    d = rs.describe(values[groups == s])
    print(f"{s:3s} {d['mean']:6.2f} +- {d['sd']:.2f}")

anova, dunnett = rs.variant_test(values, groups, control="WT")
print(f"\nsystems differ? ICC {anova['icc']:.2f}, p {rs.fmt_p(anova['p'])}")
for r in dunnett.itertuples():
    print(f"  {r.variant} vs WT: diff {r.diff:+.2f}, d {r.d:+.2f}, p {rs.fmt_p(r.p)}")

c = rs.spearman([values[groups == s].mean() for s in systems], [assay[s] for s in systems])
print(f"\nMD vs assay: Spearman r {c['r']:.2f} ({c['strength']}), p {rs.fmt_p(c['p'])}")
print(f"replicas needed for a 0.5 difference: {rs.replicas_needed(anova['sd_replica'], 0.5)}")
