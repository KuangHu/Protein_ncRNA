#!/bin/bash
#SBATCH --job-name=pn_triage
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=01:00:00
#SBATCH --output=logs/pn_triage_%j.out
#SBATCH --error=logs/pn_triage_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# BLIND. Coordinates, ORF geometry and each block's own fold. Benchmark columns
# are annotation only and are joined after every other column is computed.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CAND=${CAND:-/global/scratch/users/kh36969/protein_ncrna/step4_candidates}
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
K3=${K3:-/global/scratch/users/kh36969/protein_ncrna/refsweep_26723281/k3}
OUT=${OUT:-$CAND/triage}

$PY -u scripts/triage_report.py --candidates "$CAND/step4_candidates.tsv" \
    --blocks "$K3" --windows "$WIN" \
    --fold "$CAND/fold_all080" --padding "$CAND/fasta_all080/block_padding.tsv" \
    --out "$OUT" --merge-window "${MERGE:-500}" --exclude-benchmark-only
