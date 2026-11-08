#!/usr/bin/env python3
"""Residue-pair contact frequencies per interaction group from a GetContacts per-frame contact file."""
import argparse
import re
from collections import defaultdict

# GetContacts interaction codes -> groups
GROUPS = {
    "sb": ["sb"], "pc": ["pc"], "ps": ["ps"], "ts": ["ts"], "vdw": ["vdw"], "hp": ["hp"],
    "hb_direct": ["hb", "hbbb", "hbsb", "hbss", "hblb", "hbls"],
    "hb_water": ["wb", "wb2", "lwb", "lwb2"],
}
RAW_TO_GROUP = {code: group for group, codes in GROUPS.items() for code in codes}


def residue(atom):
    # 'A:ARG:120:NH1' -> key 'A:120' and name 'ARG'. Keyed by position so that a point mutant
    # can later be compared with the wild type at the same position.
    chain, resname, resid = atom.split(":")[:3]
    return f"{chain}:{resid}", resname


def count_contacts(path, skip_frames=0):
    counts = defaultdict(int)
    last_frame = {}
    resnames = {}
    total_frames = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                m = re.search(r"total_frames:(\d+)", line)
                if m:
                    total_frames = int(m.group(1))
                continue
            parts = line.rstrip("\n").split("\t")
            frame = int(parts[0])
            if frame < skip_frames:            # equilibration
                continue
            pos1, res1 = residue(parts[2])
            pos2, res2 = residue(parts[3])
            if pos1 == pos2:                   # contact within one residue
                continue
            resnames[pos1], resnames[pos2] = res1, res2
            pos1, pos2 = sorted((pos1, pos2))
            for group in ["all", RAW_TO_GROUP.get(parts[1])]:
                key = (group, pos1, pos2)
                # several atom contacts between the same two residues count once per frame.
                # Counting here, and not with GetContacts' get_contact_frequencies.py, because that
                # script missed the residue pair written first in each frame.
                if group and last_frame.get(key) != frame:
                    counts[key] += 1
                    last_frame[key] = frame
    if total_frames is None:
        raise ValueError(f"{path}: no total_frames in the header")
    return counts, resnames, total_frames - skip_frames


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("contacts", help="output of get_dynamic_contacts.py")
    p.add_argument("output", help="output TSV")
    p.add_argument("--skip-frames", type=int, default=0, help="leave out the first N frames")
    a = p.parse_args()

    counts, resnames, n_frames = count_contacts(a.contacts, a.skip_frames)
    with open(a.output, "w") as fh:
        fh.write(f"# analysed_frames:{n_frames}\n")
        fh.write("group\tposition_1\tposition_2\tresname_1\tresname_2\tframes\tfrequency\n")
        for (group, pos1, pos2), n in sorted(counts.items()):
            fh.write(f"{group}\t{pos1}\t{pos2}\t{resnames[pos1]}\t{resnames[pos2]}\t{n}\t{n / n_frames:.4f}\n")
    print(f"{a.contacts}: {n_frames} frames, {sum(g == 'all' for g, _, _ in counts)} residue pairs")


if __name__ == "__main__":
    main()
