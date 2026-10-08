#!/bin/bash
#SBATCH --job-name=pn_strat
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=12:00:00
#SBATCH --output=logs/pn_strat_%j.out
#SBATCH --error=logs/pn_strat_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Full re-run from step 3 after replacing every order-dependent subset with a
# deterministic one: diversity-spread clade cap, hash-sampled identity estimate,
# hash-tie-broken reference picking. Fresh directories throughout; the previous
# results stay on disk for comparison.
__conda_root=${__conda_root:-/global/home/users/kh36969/miniconda3}
. "$__conda_root/etc/profile.d/conda.sh"
conda activate rnafold
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
S=/global/scratch/users/kh36969/protein_ncrna
CEN=$S/census_26708202/corrected
NBR=$S/nbr_26710278/neighborhood.tsv
J=${SLURM_JOB_ID}
CL=${CL:-$S/clades_resample_26738463}
WIN=$S/windows_resample_$J
BLK=$S/blocks_resample_$J
CAND=$S/step4_candidates_resample_$J
BENCH=$S/known_retron_benchmark_resample_$J
export TMPDIR=/tmp/pnrs_$J; mkdir -p "$TMPDIR" "$WIN" "$BLK" "$CAND"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== reusing step 3: $CL (hash-sampled identity, 37 clades)"

echo "=== step 3b: windows with the diversity-spread cap"
$PY -u scripts/extract_windows.py --census "$CEN" \
    --clades "$CL/selected_clades.tsv" --members "$CL/cluster_members.tsv" \
    --out "$WIN" --workers "${SLURM_CPUS_PER_TASK}" --proteins \
    --max-per-clade "${MAX_PER_CLADE:-500}" --cap-mode "${CAP_MODE:-stratified}" --seed 1

echo "=== null"
$PY -u scripts/permutation_null.py --windows "$WIN" --out "$WIN/decoys" \
    --replicates 3 --seed 1 --min-source-clades 4 --max-source-fraction 0.5

echo "=== step 4"
$PY -u scripts/discover_blocks.py --windows "$WIN" --decoys "$WIN/decoys" \
    --proteins "$WIN/proteins" --out "$BLK" --references 3 \
    --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR"

echo "=== benchmark (annotation only)"
mkdir -p "$BENCH"; cp -r "$S"/known_retron_benchmark_26714577/* "$BENCH"/ 2>/dev/null || true
$PY -u scripts/benchmark_clade_overlap.py --bench "$BENCH" --clades "$CL" \
    --out "$BENCH" --blocks "$BLK" --windows "$WIN"

echo "=== candidate table"
$PY -u scripts/candidate_table.py --blocks "$BLK" --out "$CAND" \
    --benchmark "$BENCH" --decoys "$WIN/decoys"

echo "=== step 5 + 6"
$PY -u scripts/export_blocks.py --candidates "$CAND/step4_candidates.tsv" \
    --windows "$WIN" --out "$CAND/fasta_all080" --soft-orf-padding-bp 75 \
    --min-score 0.80 --exclude-benchmark-only --tmp "$TMPDIR"
$PY -u scripts/fold_blocks.py --fasta "$CAND/fasta_all080" \
    --out "$CAND/fold_all080" --which padded --threads "${SLURM_CPUS_PER_TASK}"

echo "=== triage"
$PY -u scripts/triage_report.py --candidates "$CAND/step4_candidates.tsv" \
    --blocks "$BLK" --windows "$WIN" --fold "$CAND/fold_all080" \
    --padding "$CAND/fasta_all080/block_padding.tsv" \
    --out "$CAND/triage" --merge-window 500 --exclude-benchmark-only

echo "=== outputs: $CL $WIN $BLK $CAND $BENCH"
