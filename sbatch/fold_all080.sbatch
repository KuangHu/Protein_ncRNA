#!/bin/bash
#SBATCH --job-name=pn_fold080
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --time=10:00:00
#SBATCH --output=logs/pn_fold080_%j.out
#SBATCH --error=logs/pn_fold080_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# BLIND. Soft-padded export -> MAFFT -> RNAalifold -> R-scape on each block's own
# members, plus a locus report for the two GroupII_RT blocks. No curated RNA
# family, no covariance model from any database, is read anywhere here.
__conda_root=${__conda_root:-/global/home/users/kh36969/miniconda3}
. "$__conda_root/etc/profile.d/conda.sh"
conda activate rnafold
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CAND=${CAND:-/global/scratch/users/kh36969/protein_ncrna/step4_candidates}
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
K3=${K3:-/global/scratch/users/kh36969/protein_ncrna/refsweep_26723281/k3}
FA=${FA:-$CAND/fasta_all080}
OUT=${OUT:-$CAND/fold_all080}
LOC=${LOC:-$CAND/locus_groupII}
export TMPDIR=/tmp/pnfold080_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# All non-Cas3 candidates at the frozen 0.80 score floor, soft ORF padding.
$PY -u scripts/export_blocks.py --candidates "$CAND/step4_candidates.tsv" \
    --windows "$WIN" --out "$FA" --soft-orf-padding-bp "${PAD:-75}" \
    --min-score "${MINSCORE:-0.80}" --exclude-benchmark-only --tmp "$TMPDIR"

$PY -u scripts/fold_blocks.py --fasta "$FA" --out "$OUT" \
    --which padded --threads "${SLURM_CPUS_PER_TASK}"

# Locus report for the two GroupII_RT blocks of the CAUXVC010000002.1 clade.
C='GroupII_RT__thr0.35__GCA_963568795.1_CAUXVC010000002.1_22007-23260__'
$PY -u scripts/locus_report.py --blocks "$K3" --windows "$WIN" --out "$LOC" \
    --fold "$OUT" --tmp "$TMPDIR" \
    --block-id "${C}__B000" --block-id "${C}__B001"
