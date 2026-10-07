#!/bin/bash
#SBATCH --job-name=pn_discovery
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --output=logs/pn_anchors_%j.out
#SBATCH --error=logs/pn_anchors_%j.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
OUT=/global/scratch/users/kh36969/protein_ncrna/${SLURM_JOB_NAME}_${SLURM_JOB_ID}
# prodigal and hmmsearch write many small temp files; keep them node-local.
SCRATCH=/tmp/pn_${SLURM_JOB_ID}
mkdir -p "$OUT" "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

# Discovery mode: carves and stores every window (~10 kb/anchor). Run the
# census first -- it decides which families are worth this.
GENOMES=${GENOMES:-/global/scratch/users/kh36969/protein_ncrna/corpus/real_genomes.tsv}
if [ ! -f "$GENOMES" ]; then
    echo "missing $GENOMES -- run sbatch/genomes_db_report.sbatch first" >&2
    exit 1
fi

# prodigal dominates runtime (~10 s/genome, single-threaded), so one worker per
# core and chunksize 1 keeps the pool balanced across uneven genome sizes.
$PY -u scripts/collect_anchors.py \
    --genomes "$GENOMES" \
    --out "$OUT" --windows all \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --scratch "$SCRATCH"

echo "windows: $(wc -l < "$OUT/windows.jsonl")"
cat "$OUT/summary.json"
