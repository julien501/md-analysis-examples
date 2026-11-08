# Replica statistics

Frames of one simulation are correlated in time, so here every test uses the replica
(one run) as the unit.

`replica_stats.py`

| Function | Question |
|---|---|
| `variant_anova` | do the systems differ at all? (ICC and p of a one-way ANOVA) |
| `variant_test`, `dunnett_vs_wt` | which systems differ from the reference? (Dunnett, only after a significant ANOVA) |
| `spearman`, `strength` | does an MD descriptor follow an experimental readout? |
| `free_energy` | free-energy profile from a distribution |
| `replicas_needed` | how many replicas would detect a given difference? |
| `bh` | Benjamini–Hochberg q values for many descriptors |

`replica_noise.py` measures, with the same code for every data set, how much the replicas of
one system differ in contacts and RMSF, and compares two replicas of one system with two runs
of different systems.

```bash
python example.py
python -m pytest tests
python replica_noise.py runs.csv out/ --receptor-range 1 300 --ligand-range 1 30 --window 100 500
```

References: Grossfield et al. 2018, LiveCoMS 1:5067 (uncertainty in simulations);
Schober et al. 2018, Anesth Analg 126:1763 (correlation strength).
