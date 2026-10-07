#!/bin/bash
#SBATCH --job-name=pn_nbr
#SBATCH --partition=lr6
#SBATCH --account=pc_rubinlab
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --array=0-7
#SBATCH --output=logs/pn_nbr_%A_%a.out
#SBATCH --error=logs/pn_nbr_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA

PY=/global/home/users/kh36969/.conda/envs/opfi/bin/python
CEN=${CEN:-/global/scratch/users/kh36969/protein_ncrna/census_26708202/corrected}
OUT=/global/scratch/users/kh36969/protein_ncrna/nbr_${SLURM_ARRAY_JOB_ID}
SCRATCH=/tmp/pnnbr_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$OUT" "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

# Cas3 is provisional until cascade genes are confirmed beside it; RT needs
# effector neighbours to separate defence RTs from bare ones.
$PY -u scripts/neighborhood_scan.py \
    --census "$CEN" \
    --families Cas3 RT \
    --out "$OUT/shard_${SLURM_ARRAY_TASK_ID}" \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT:-8}" \
    --workers "${SLURM_CPUS_PER_TASK}" --scratch "$SCRATCH"
