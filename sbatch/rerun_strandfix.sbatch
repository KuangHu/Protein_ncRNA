#!/bin/bash
#SBATCH --job-name=pn_strandfix
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=10:00:00
#SBATCH --output=logs/pn_strandfix_%j.out
#SBATCH --error=logs/pn_strandfix_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Full re-run of steps 4-6 after the minus-strand member-trim fix in
# trim_member(). Writes to fresh directories so the pre-fix results stay on
# disk for comparison; nothing is overwritten in place.
__conda_root=${__conda_root:-/global/home/users/kh36969/miniconda3}
. "$__conda_root/etc/profile.d/conda.sh"
conda activate rnafold
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
S=/global/scratch/users/kh36969/protein_ncrna
WIN=$S/windows_26714415
BENCH=$S/known_retron_benchmark_26714577
BLK=$S/blocks_strandfix_${SLURM_JOB_ID}
CAND=$S/step4_candidates_strandfix_${SLURM_JOB_ID}
export TMPDIR=/tmp/pnsf_${SLURM_JOB_ID}; mkdir -p "$TMPDIR" "$BLK" "$CAND"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== step 4: block calling, real + 129 decoys, frozen k=3 / r=0.5"
$PY -u scripts/discover_blocks.py --windows "$WIN" --decoys "$WIN/decoys" \
    --proteins "$WIN/proteins" --out "$BLK" --references 3 \
    --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR"

echo "=== benchmark overlap (annotation only, written under \$BENCH2)"
BENCH2=$S/known_retron_benchmark_strandfix_${SLURM_JOB_ID}
mkdir -p "$BENCH2"; cp -r "$BENCH"/* "$BENCH2"/ 2>/dev/null || true
$PY -u scripts/benchmark_clade_overlap.py --bench "$BENCH2" \
    --clades "$S/clades_v1" --out "$BENCH2" --blocks "$BLK" --windows "$WIN"

echo "=== candidate table"
$PY -u scripts/candidate_table.py --blocks "$BLK" --out "$CAND" --benchmark "$BENCH2"

echo "=== step 5: soft-padded export, all non-Cas3 >= 0.80"
$PY -u scripts/export_blocks.py --candidates "$CAND/step4_candidates.tsv" \
    --windows "$WIN" --out "$CAND/fasta_all080" --soft-orf-padding-bp 75 \
    --min-score 0.80 --exclude-benchmark-only --tmp "$TMPDIR"

echo "=== step 6: fold"
$PY -u scripts/fold_blocks.py --fasta "$CAND/fasta_all080" \
    --out "$CAND/fold_all080" --which padded --threads "${SLURM_CPUS_PER_TASK}"

echo "=== triage"
$PY -u scripts/triage_report.py --candidates "$CAND/step4_candidates.tsv" \
    --blocks "$BLK" --windows "$WIN" --fold "$CAND/fold_all080" \
    --padding "$CAND/fasta_all080/block_padding.tsv" \
    --out "$CAND/triage" --merge-window 500 --exclude-benchmark-only

echo "=== outputs: $BLK  $CAND  $BENCH2"
