#!/bin/bash
#SBATCH --job-name=pn_mbench
#SBATCH --partition=es1
#SBATCH --account=pc_rubinlab
#SBATCH --qos=es_normal
#SBATCH --array=0-7%8
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --gres=gpu:A40:1 --time=12:00:00
#SBATCH --output=logs/pn_mbench_%A_%a.out
#SBATCH --error=logs/pn_mbench_%A_%a.err
set -euo pipefail
cd /global/home/users/kh36969/tools/Protein_ncRNA
# Does Minerva separate the three families that carry a known ncRNA -- IS110
# (bridge RNA), group II intron (ribozyme), IS200/605 (omegaRNA) -- from IS
# families that carry none: IS1, IS4/IS5/IS10, IS3?
#
# 42 RNA+ bags vs 137 RNA- bags, scored on ELEMENT-TRIMMED windows so the
# regions compared are in-element non-coding on both sides. The family labels
# are used to READ the answer, never to produce it: Minerva sees only sequence,
# and nothing here re-enters the block caller, the decoys or the floor.
PY=/global/home/users/kh36969/.conda/envs/minerva/bin/python
export PYTHONNOUSERSITE=1 HF_HOME=/global/scratch/users/kh36969/hf_cache
E=/global/scratch/users/kh36969/protein_ncrna/expand
$PY -I -u scripts/minerva_score.py \
    --windows "$E/elem_windows" \
    --bags data/minerva_bench_bags.txt \
    --blocks "$E/elem_run/blocks/blocks.tsv" \
    --out "$E/minerva/bench" --max-cases 10 \
    --shard "${SLURM_ARRAY_TASK_ID}" --nshards "${SLURM_ARRAY_TASK_COUNT}"
