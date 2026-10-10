#!/bin/bash
#SBATCH --job-name=pn_fs4
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --array=0-15%16
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --time=48:00:00
#SBATCH --output=logs/pn_fs4_%A_%a.out
#SBATCH --error=logs/pn_fs4_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Step 4 over the whole corpus, sharded by clade. Every clade is scored
# independently -- references come from its own windows and the anchor-protein
# pool handed to each task is the same one -- so splitting the clade files
# across tasks gives exactly what one job would have produced. The 117-subject
# run took 5h24m single-job; this is 3.4x the subjects.
#
# Settings are frozen and identical to that run: --references 3
# --min-recurrence 0.5. Widening the corpus is the change; the caller is not
# being retuned, or nothing here would be comparable to what came before.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
E=/global/scratch/users/kh36969/protein_ncrna/expand
P=$E/full_run
S=$P/shard_$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")
export TMPDIR=/tmp/pnfs4_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
rm -rf "$S"; mkdir -p "$S/real" "$S/decoys"

# Deterministic round-robin over the sorted clade list, real and decoy split
# separately so every task carries a comparable mix of both.
link_shard () {
  local SRC=$1 DST=$2 I=0
  for F in "$SRC"/*.fasta; do
    [ -e "$F" ] || continue
    if [ $((I % SLURM_ARRAY_TASK_COUNT)) -eq "${SLURM_ARRAY_TASK_ID}" ]; then
      B=$(basename "$F" .fasta)
      ln -sf "$SRC/$B.fasta" "$DST/$B.fasta"
      [ -f "$SRC/$B.jsonl" ] && ln -sf "$SRC/$B.jsonl" "$DST/$B.jsonl"
    fi
    I=$((I + 1))
  done
}
link_shard "$P/real" "$S/real"
link_shard "$P/decoys" "$S/decoys"
echo "shard ${SLURM_ARRAY_TASK_ID}: $(ls "$S"/real/*.fasta 2>/dev/null | wc -l) real, $(ls "$S"/decoys/*.fasta 2>/dev/null | wc -l) decoy clades"

$PY -u scripts/discover_blocks.py --windows "$S/real" --decoys "$S/decoys" \
    --proteins "$P/proteins" --out "$S/blocks" --references 3 \
    --min-recurrence 0.5 --workers 8 --threads 4 --tmp "$TMPDIR"
