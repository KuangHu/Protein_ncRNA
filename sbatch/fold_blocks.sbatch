#!/bin/bash
#SBATCH --job-name=pn_fold
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --time=08:00:00
#SBATCH --output=logs/pn_fold_%j.out
#SBATCH --error=logs/pn_fold_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# BLIND. MAFFT -> RNAalifold -> R-scape on the block's own members. No curated
# RNA family, no covariance model from any database, is read here.
__conda_root=${__conda_root:-/global/home/users/kh36969/miniconda3}
. "$__conda_root/etc/profile.d/conda.sh"
conda activate rnafold          # R-scape / Infernal / RNAalifold / MAFFT
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CAND=${CAND:-/global/scratch/users/kh36969/protein_ncrna/step4_candidates}
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
FA=${FA:-$CAND/fasta_hc}
OUT=${OUT:-$CAND/fold_hc}
export TMPDIR=/tmp/pnfold_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# High-confidence, non-Cas3 only.
$PY -u scripts/export_blocks.py --candidates "$CAND/step4_candidates.tsv" \
    --windows "$WIN" --out "$FA" --pad "${PAD:-100}" \
    --min-score "${MINSCORE:-0.90}" --exclude-benchmark-only --tmp "$TMPDIR"

$PY -u scripts/fold_blocks.py --fasta "$FA" --out "$OUT" \
    --which padded --threads "${SLURM_CPUS_PER_TASK}"
