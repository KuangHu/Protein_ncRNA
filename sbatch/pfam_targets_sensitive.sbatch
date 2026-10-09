#!/bin/bash
#SBATCH --job-name=pn_pfsens
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --time=04:00:00
#SBATCH --output=logs/pn_pfsens_%j.out
#SBATCH --error=logs/pn_pfsens_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Sensitive Pfam read-out for the 18 hard-null targets. The production scan used
# --cut_ga, the curated per-family threshold, which is the right bar for calling
# a domain present. Five of these bags cleared it with nothing at all, and for an
# orphan like CDS07978 "no domain above the curated threshold" and "no homology
# to any Pfam family" are very different statements. So this is a second pass at
# -E 1 with no gathering cutoff: everything it finds is BELOW the curated bar and
# is reported as a lead, never as an annotation.
#
# Still outside BLIND_CHAIN -- this reads the discovery out, it does not steer it.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
HMMSCAN=/global/home/users/kh36969/.conda/envs/opfi/bin/hmmscan
PFAM=/global/scratch/users/kh36969/pfam_db/Pfam-A.hmm
O=$E/pfam_sensitive; mkdir -p "$O"
echo "=== hmmscan -E 1 (no --cut_ga): $(grep -c '^>' "$E/hpool/targets.faa") reps of 18 target bags"
$HMMSCAN --cpu "${SLURM_CPUS_PER_TASK}" -E 1 --domE 1 \
    --domtblout "$O/targets_sensitive.domtbl" -o /dev/null \
    "$PFAM" "$E/hpool/targets.faa"
echo "=== also phmmer vs Pfam-A seed-free: reps vs all 495 bag reps already done"
wc -l "$O/targets_sensitive.domtbl"
