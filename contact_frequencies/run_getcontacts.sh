#!/bin/bash
#SBATCH --job-name=getcontacts
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
# Contacts between two selections for one run, then residue-pair frequencies.
# Usage: sbatch run_getcontacts.sh TOPOLOGY TRAJECTORY OUT_PREFIX SELECTION_1 SELECTION_2 [SKIP_FRAMES]
set -euo pipefail

TOPOLOGY=$1 TRAJECTORY=$2 OUT_PREFIX=$3 SELE1=$4 SELE2=$5 SKIP_FRAMES=${6:-0}
GETCONTACTS_DIR=${GETCONTACTS_DIR:?set GETCONTACTS_DIR to your GetContacts folder}
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$0")}      # sbatch runs a copy of this script elsewhere

mkdir -p "$(dirname "$OUT_PREFIX")"
python "$GETCONTACTS_DIR/get_dynamic_contacts.py" --topology "$TOPOLOGY" --trajectory "$TRAJECTORY" \
    --sele "$SELE1" --sele2 "$SELE2" --itypes all --cores "${SLURM_CPUS_PER_TASK:-1}" \
    --output "${OUT_PREFIX}_contacts.tsv"
python "$HERE/contact_frequencies.py" "${OUT_PREFIX}_contacts.tsv" "${OUT_PREFIX}_frequencies.tsv" \
    --skip-frames "$SKIP_FRAMES"
