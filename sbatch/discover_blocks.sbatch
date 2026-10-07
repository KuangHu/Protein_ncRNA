#!/bin/bash
#SBATCH --job-name=pn_blocks
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --output=logs/pn_blocks_%j.out
#SBATCH --error=logs/pn_blocks_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

# BLIND. No RetronDB, no Rfam, no covariance model is read by anything here.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
OUT=${OUT:-/global/scratch/users/kh36969/protein_ncrna/blocks_${SLURM_JOB_ID}}
export TMPDIR=/tmp/pnblocks_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$OUT"
trap 'rm -rf "$TMPDIR"' EXIT

# Decoys are a permutation of the same windows, so they cost no extraction.
[ -d "$WIN/decoys" ] || $PY -u scripts/permutation_null.py \
    --windows "$WIN" --out "$WIN/decoys" --replicates "${REPLICATES:-3}" --seed 1

$PY -u scripts/discover_blocks.py \
    --windows "$WIN" --decoys "$WIN/decoys" --out "$OUT" \
    --workers 8 --threads 4 --tmp "$TMPDIR" --proteins "$WIN/proteins"
