#!/bin/bash
#SBATCH --job-name=pn_blocks2
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --output=logs/pn_blocks2_%j.out
#SBATCH --error=logs/pn_blocks2_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# BLIND. Reference choice uses anchor-protein diversity only -- no RNA, no
# known-RNA coordinates. min-recurrence stays at the 0.5 operating point.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
WIN=${WIN:-/global/scratch/users/kh36969/protein_ncrna/windows_26714415}
OUT=${OUT:-/global/scratch/users/kh36969/protein_ncrna/blocks_v2_${SLURM_JOB_ID}}
export TMPDIR=/tmp/pnblocks2_${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$OUT"
trap 'rm -rf "$TMPDIR"' EXIT

$PY -u scripts/discover_blocks.py \
    --windows "$WIN" --decoys "$WIN/decoys" --proteins "$WIN/proteins" \
    --out "$OUT" --references "${REFS:-12}" --min-recurrence 0.5 \
    --workers 8 --threads 4 --tmp "$TMPDIR"
