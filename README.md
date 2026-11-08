# MD analysis examples

Analysis scripts from my molecular dynamics work on membrane proteins (GPCR–ligand complexes
and a bacterial outer membrane complex). Code only, no simulation data.

| Folder | What it does |
|---|---|
| [`contact_frequencies/`](contact_frequencies/) | contact and hydrogen bond frequencies from GetContacts, averaged over replicas, and their change vs a reference system |
| [`replica_statistics/`](replica_statistics/) | statistics with the replica as the unit, and how much the replicas of one system differ |
| [`lipid_contacts/`](lipid_contacts/) | lipid–protein contacts per chain and lipid type, as tables, heatmaps and B-factor structures |

```bash
pip install -r requirements.txt
cd replica_statistics && python example.py && python -m pytest tests
```

Julien Reißmann
