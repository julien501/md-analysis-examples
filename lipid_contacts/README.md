# Lipid contacts

Which residues of a membrane protein are in contact with which lipid type, per chain and
averaged over chains and replicates. Input: the GetContacts residue frequency file of each
replicate.

```bash
python lipid_contacts.py --freq rep1.tsv --freq rep2.tsv --freq rep3.tsv \
    --lipids POPE,POPG --chains A,B,C --structure system.pdb --outdir out
chimerax --offscreen --script "color_bfactor_chimerax.py out/bfactor_all.pdb"
```

Output: `per_chain.csv`, `chain_mean.csv`, heatmaps per chain and for the chain mean, and one
PDB per lipid type (plus `all`) with B-factor = 100 × contact frequency, coloured by the
ChimeraX script.
