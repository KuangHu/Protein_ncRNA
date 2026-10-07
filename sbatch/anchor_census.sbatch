#!/bin/bash
#SBATCH --job-name=pn_census
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/pn_census_%A_%a.out
#SBATCH --error=logs/pn_census_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CORPUS=/global/scratch/users/kh36969/protein_ncrna/corpus
OUT=/global/scratch/users/kh36969/protein_ncrna/census_${SLURM_ARRAY_JOB_ID}
SCRATCH=/tmp/pn_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$OUT" "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

GENOMES=${GENOMES:-$CORPUS/real_genomes.tsv}
if [ ! -f "$GENOMES" ]; then
    echo "missing $GENOMES -- run sbatch/genomes_db_report.sbatch first" >&2
    exit 1
fi

NSHARDS=${SLURM_ARRAY_TASK_COUNT:-16}
echo "shard ${SLURM_ARRAY_TASK_ID}/${NSHARDS} of $(($(wc -l < "$GENOMES") - 1)) real genomes"

# --windows none: the census needs flank *feasibility*, not 10 kb per anchor.
# prodigal dominates (~19 s/genome/core), so one worker per core.
$PY -u scripts/collect_anchors.py \
    --genomes "$GENOMES" \
    --out "$OUT/shard_${SLURM_ARRAY_TASK_ID}" \
    --windows none \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${NSHARDS}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --scratch "$SCRATCH"

echo "shard ${SLURM_ARRAY_TASK_ID} done"
