# Contact frequencies

Which residue contacts between two chains (e.g. ligand and receptor) does a point mutation
break or create, over several replicas per system?

```
run_getcontacts.sh       trajectory -> per-frame contacts (GetContacts) -> frequencies
contact_frequencies.py   per-frame contacts -> frequency per residue pair and interaction group
frequency_deltas.py      frequencies of all runs -> mean and SD over replicas, difference to the reference
plot_deltas.py           differences -> one heatmap per system
```

Interaction groups: salt bridges, π-cation, π-stacking, T-stacking, van der Waals,
hydrophobic, direct and water-mediated hydrogen bonds, and all together.

```bash
export GETCONTACTS_DIR=/path/to/getcontacts
sbatch run_getcontacts.sh wt_1.pdb wt_1.xtc out/wt_1 "chain A" "chain B" 1000
python frequency_deltas.py runs.csv deltas.tsv --reference WT
python plot_deltas.py deltas.tsv hb_direct.svg --group hb_direct --rows B --cols A
```

`runs.csv` lists the frequency file of every run:

```
system,replica,path
WT,1,out/wt_1_frequencies.tsv
WT,2,out/wt_2_frequencies.tsv
M1,1,out/m1_1_frequencies.tsv
```
