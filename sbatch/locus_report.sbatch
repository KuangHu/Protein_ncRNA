#!/bin/bash
#SBATCH --job-name=pn_locus
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=02:00:00
#SBATCH --output=logs/pn_locus_%j.out
#SBATCH --error=logs/pn_locus_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# BLIND. Per-member ORF geometry around one block. No curated RNA model enters
# here; whether the locus *is* a known element is asked separately, in the
# unblind branch, after this table exists.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
S=/global/scratch/users/kh36969/protein_ncrna
RUN=${RUN:-final_26739540}
BLK=${BLK:-$S/blocks_$RUN}
WIN=${WIN:-$S/windows_$RUN}
CAND=${CAND:-$S/step4_candidates_$RUN}
OUT=${OUT:-$CAND/loci}
BLOCK_IDS=${BLOCK_IDS:?set BLOCK_IDS to a space-separated list of block ids}

ARGS=()
for b in $BLOCK_IDS; do ARGS+=(--block-id "$b"); done

$PY -u scripts/locus_report.py --blocks "$BLK" --windows "$WIN" \
    --fold "$CAND/fold_all080" --out "$OUT" "${ARGS[@]}"
