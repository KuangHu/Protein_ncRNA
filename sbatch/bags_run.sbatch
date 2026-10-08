#!/bin/bash
#SBATCH --job-name=pn_bags
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=24:00:00
#SBATCH --output=logs/pn_bags_%j.out
#SBATCH --error=logs/pn_bags_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# The FROZEN step-4 chain (tag v1.0-frozen) applied to the MGE insertion
# catalogue from fna_based_mgefinder_project. Nothing in the caller changes;
# only the corpus does. The published score floors are NOT applied as a filter
# here -- inserts are ~1 kb against the 11 kb genomic windows those floors came
# from, so they are re-derived from this corpus's own decoys at the end.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
S=/global/scratch/users/kh36969/protein_ncrna
F=/global/scratch/users/kh36969/fna_ins_discovery
TAG=${TAG:-bags_$SLURM_JOB_ID}
WIN=$S/windows_$TAG
BLK=$S/blocks_$TAG
CAND=$S/step4_candidates_$TAG
NT=${SLURM_CPUS_PER_TASK}
export TMPDIR=/tmp/pnb_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== bags -> windows in the anchor frame"
$PY -u scripts/bags_to_windows.py --bags "$F/bags_v4" --cds "$F/cds_v3" \
    --inserts "$F/database_v3/insertions_inserts.fna" --out "$WIN" \
    --max-per-clade "${CAP:-500}" --min-sites "${MINSITES:-8}"

echo "=== post-cap diagnostics"
$PY -u scripts/cap_diagnostics.py --windows "$WIN" --census "$WIN" \
    --threads "$NT" --out "$WIN/cap_diagnostics.tsv" || echo "(non-fatal)"

echo "=== within-corpus permutation null"
$PY -u scripts/permutation_null.py --windows "$WIN" --out "$WIN/decoys" \
    --replicates "${REPS:-3}" --seed 1 --min-source-clades 4 \
    --max-source-fraction 0.5

echo "=== step 4 (frozen parameters)"
$PY -u scripts/discover_blocks.py --windows "$WIN" --decoys "$WIN/decoys" \
    --proteins "$WIN/proteins" --out "$BLK" --references 3 \
    --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR"

echo "=== candidate table"
$PY -u scripts/candidate_table.py --blocks "$BLK" --out "$CAND" \
    --decoys "$WIN/decoys"

echo "=== export + fold"
$PY -u scripts/export_blocks.py --candidates "$CAND/step4_candidates.tsv" \
    --windows "$WIN" --out "$CAND/fasta_all" --soft-orf-padding-bp 75 \
    --min-score "${FOLD_FLOOR:-0.80}" --tmp "$TMPDIR"
$PY -u scripts/fold_blocks.py --fasta "$CAND/fasta_all" \
    --out "$CAND/fold_all" --which padded --threads "$NT"

echo "=== triage"
$PY -u scripts/triage_report.py --candidates "$CAND/step4_candidates.tsv" \
    --blocks "$BLK" --windows "$WIN" --fold "$CAND/fold_all" \
    --padding "$CAND/fasta_all/block_padding.tsv" \
    --out "$CAND/triage" --merge-window 500
$PY -u scripts/discovery_shortlist.py --triage "$CAND/triage" \
    --exclude-families NONE

echo "=== outputs: $WIN $BLK $CAND"
