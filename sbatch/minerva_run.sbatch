#!/bin/bash
#SBATCH --job-name=pn_mrun
#SBATCH --partition=es1
#SBATCH --account=pc_rubinlab
#SBATCH --qos=es_normal
#SBATCH --array=0-7%8
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --gres=gpu:A40:1 --time=12:00:00
#SBATCH --output=logs/pn_mrun_%A_%a.out
#SBATCH --error=logs/pn_mrun_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# The 175 bags the recurrence method cannot serve: scored and blockless, or
# below the 5 distinct anchors a block needs. For these, recurrence is not a
# weak result, it is undefined. Minerva asks a per-locus question instead.
#
# Proposer only. Nothing here re-enters the block caller, the decoys or the
# floor, and every column is prefixed minerva_.
PY=/global/home/users/kh36969/.conda/envs/minerva/bin/python
export PYTHONNOUSERSITE=1 HF_HOME=/global/scratch/users/kh36969/hf_cache
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -I -u scripts/minerva_score.py \
    --windows "$E/strat_windows" "$E/hpool_windows" "$E/pilot_windows" "$E/full_windows" \
    --bags "$E/minerva/datapoor_bags.txt" \
    --blocks "$E/full_run/blocks/blocks.tsv" \
    --out "$E/minerva/datapoor" --max-cases 10 \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
