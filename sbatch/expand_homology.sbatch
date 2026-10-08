#!/bin/bash
#SBATCH --job-name=pn_expand
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --time=24:00:00
#SBATCH --array=0-31
#SBATCH --output=logs/pn_expand_%A_%a.out
#SBATCH --error=logs/pn_expand_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Arm 2: recover copies of each bag's protein that the empty-allele filter never
# saw. Prodigal every genome, MMseqs every ORF against the bag representatives,
# and keep what matches -- with no requirement that an empty allele exist.
# Blind: sequence search only, no HMM, no IS library, no RNA model.
PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
S=/global/scratch/users/kh36969/protein_ncrna
TAG=${TAG:-expand}
Q=${Q:-$S/$TAG/bag_reps.faa}
OUT=$S/$TAG/anchors
export TMPDIR=/tmp/pnx_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT

$PY -u scripts/expand_by_homology.py \
    --genomes out/real_genomes.tsv --queries "$Q" --out "$OUT" \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}" \
    --workers "${SLURM_CPUS_PER_TASK}" --scratch "$TMPDIR" \
    --min-identity "${MINID:-0.30}" --min-coverage "${MINCOV:-0.70}" \
    --evalue "${EVAL:-1e-10}"
