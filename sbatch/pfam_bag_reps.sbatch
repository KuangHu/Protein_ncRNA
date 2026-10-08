#!/bin/bash
#SBATCH --job-name=pn_pfam
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=08:00:00
#SBATCH --output=logs/pn_pfam_%j.out
#SBATCH --error=logs/pn_pfam_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# POST-HOC annotation of the Arm 2 anchor proteins. This is interpretation, not
# discovery: the blind pipeline anchored on `cds_cluster_id` and never saw a
# Pfam label. Nothing written here is read by any script in BLIND_CHAIN, and a
# Pfam hit must never decide which bag is carried forward -- that is the whole
# reason the discovery run came first.
#
# Pfam gathering thresholds (--cut_ga) rather than an e-value: they are the
# curators' own per-family cutoffs, so the hits are comparable across 20k
# models of wildly different length and information content.

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
HMMSCAN=/global/home/users/kh36969/.conda/envs/opfi/bin/hmmscan
PFAM=/global/scratch/users/kh36969/pfam_db/Pfam-A.hmm
E=/global/scratch/users/kh36969/protein_ncrna/expand
OUT=$E/pfam
mkdir -p "$OUT"

echo "=== hmmscan: $(grep -c '^>' "$E/bag_reps.faa") bag representatives vs Pfam-A"
$HMMSCAN --cpu "${SLURM_CPUS_PER_TASK}" --cut_ga \
    --domtblout "$OUT/bag_reps.domtbl" \
    -o /dev/null \
    "$PFAM" "$E/bag_reps.faa"

echo "=== per-bag summary"
$PY -u scripts/pfam_annotate.py \
    --domtbl "$OUT/bag_reps.domtbl" \
    --nr100 "$E/bag_nr100.tsv" \
    --blocks /global/scratch/users/kh36969/protein_ncrna/expand/pilot/blocks/blocks.tsv \
    --focus "$E/pilot_focus.txt" \
    --out "$OUT"

echo "=== outputs: $OUT"
